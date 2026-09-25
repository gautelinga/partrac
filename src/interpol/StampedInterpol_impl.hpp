#ifndef __STAMPEDINTERPOL_IMPL_HPP
#define __STAMPEDINTERPOL_IMPL_HPP

// StampedInterpol's loading and update, included by each format's translation
// unit and instantiated there, two cells a unit; the evaluation is in
// StampedInterpol_eval.hpp

#include "Error.hpp"
#include "StampedInterpol.hpp"
#include "mesh_tables.hpp"
#include "p12_eval.hpp"
#include "phase_timing.hpp"
#include <array>
#include <algorithm>
#include <cassert>
#include <iostream>

template<typename Cell, typename Format>
StampedInterpol<Cell, Format>::StampedInterpol(const std::string& infilename)
  : MeshCore<Cell>(infilename)
{
  partrac::phase_begin("load");
  dolfin_params = partrac::parse_file_or_exit(Format::schema(D), infilename);

  const std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  read_mesh_params();
  include_phi = dolfin_params.template get<bool>("include_phi");
  // P2 near walls: edge or none
  wall_p2_ = dolfin_params.template get<std::string>("wall_p2") == "edge" ? WallP2::Edge : WallP2::None;

  // The mesh, its tables, and the first stamp where the format reads it here
  fmt_.load(*this, infilename);

  if (wall_p2_ == WallP2::Edge && ncoeffs_u != Uint(Cell::n_verts))
    partrac::fail(infilename, ": wall_p2 = edge is for a P1 velocity, and this one has ",
                  ncoeffs_u, " dofs a cell");

  // Identify edge cells
  near_wall::label_cell_type(cell_type_, facet_neigh_, Cell::n_verts);
  std::cout << "Built neighbour list" << std::endl;

  if (wall_p2_ == WallP2::Edge){
    if (vclass_.empty())
      vclass_ = mesh_tables::match_periodic_vertices(coords_, nverts_, dim, periodic,
                                                     x_min, x_max, periodic_tol);
    near_wall::build_wall_edges<Cell>(topo_, coords_, ncells_, nverts_, dim, facet_neigh_,
                                      vclass_, wall_index_, wall_cells_, verbose);
    if (!stamps_.a().u.empty())
      stamps_.a().rest_tol = rest_tol(stamps_.a().u);
    partrac::phase("wall edges");
  }

  // The per-cell geometry and the tree, from the topology
  mesh_tables::build_cells_and_tree<Cell>(topo_, coords_, ncells_, nverts_, dim, cells_, tree_);
  partrac::phase_total();

  std::cout << "Setting max threads: " << omp_get_max_threads() << std::endl;
}

template<typename Cell, typename Format>
void StampedInterpol<Cell, Format>::update(const double t)
{
  const auto sp = fmt_.ts.get(t);

  if (partrac::stamp_reload(is_initialized, t, fmt_.ts.get_t_max(),
                            {t_prev, t_next}, {sp.prev.t, sp.next.t})){
    partrac::phase_begin("update");
    const auto fill = stamps_.load(fmt_.key(*this, sp.prev), fmt_.key(*this, sp.next),
                                   [&](const typename Format::Key& key, Stamp& s){
      fmt_.read(*this, key, s);
      if (wall_p2_ == WallP2::Edge)
        s.rest_tol = rest_tol(s.u);
    });
    std::cout << "Prev: Timestep = " << sp.prev.t << ", "
              << partrac::stamp_note(fill.first, fmt_.name(sp.prev)) << std::endl;
    std::cout << "Next: Timestep = " << sp.next.t << ", "
              << partrac::stamp_note(fill.second, fmt_.name(sp.next)) << std::endl;

    u_prev_ = stamps_.prev().u.data();
    u_next_ = stamps_.next().u.data();
    p_prev_ = stamps_.prev().p.data();
    p_next_ = stamps_.next().p.data();
    phi_prev_ = stamps_.prev().phi.data();
    phi_next_ = stamps_.next().phi.data();
    rest_tol_prev_ = stamps_.prev().rest_tol;
    rest_tol_next_ = stamps_.next().rest_tol;

    partrac::phase("vector read");
    partrac::phase_total();

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  t_update = t;
}

template<typename Cell, typename Format>
void StampedInterpol<Cell, Format>::freeze(const double t)
{
  const double t_min = fmt_.ts.get_t_min(), t_max = fmt_.ts.get_t_max();
  if (t < t_min || t > t_max)
    partrac::fail(this->infilename, ": the fields cannot be frozen at t = ", t, ", outside the stamps' ",
                  t_min, " to ", t_max);
  update(t);
  const double a = stamp_weight(t, t_prev, t_next);
  if (a == 1.){
    u_prev_ = u_next_;
    p_prev_ = p_next_;
    phi_prev_ = phi_next_;
    rest_tol_prev_ = rest_tol_next_;
  }
  else if (a != 0.){
    // Blend the two stamps into one
    const auto blend = [a](const std::vector<double>& prev, const std::vector<double>& next,
                           std::vector<double>& out){
      out.resize(prev.size());
      for (std::size_t i = 0; i < prev.size(); ++i) out[i] = a*next[i] + (1 - a)*prev[i];
    };
    blend(stamps_.prev().u, stamps_.next().u, frozen_.u);
    blend(stamps_.prev().p, stamps_.next().p, frozen_.p);
    blend(stamps_.prev().phi, stamps_.next().phi, frozen_.phi);
    frozen_.rest_tol = wall_p2_ == WallP2::Edge ? rest_tol(frozen_.u) : 0.;
    u_prev_ = frozen_.u.data();
    p_prev_ = frozen_.p.data();
    phi_prev_ = frozen_.phi.data();
    rest_tol_prev_ = frozen_.rest_tol;
  }
  u_next_ = u_prev_;
  p_next_ = p_prev_;
  phi_next_ = phi_prev_;
  rest_tol_next_ = rest_tol_prev_;
  // A single stamp: the weight is zero and the rate too at every t
  t_prev = t;
  t_next = t;
  std::cout << "Fields frozen at t = " << t << std::endl;
}

template<typename Cell, typename Format>
double StampedInterpol<Cell, Format>::rest_tol(const std::vector<double>& u_data) const
{
  // Round-off of the largest velocity a cell can read: an image vertex's own
  // value is never gathered, its master's is
  double u_max = 0.;
  if constexpr (Format::vertex_fields){
    for (std::size_t v = 0; v < nverts_; ++v){
      if (vclass_[v] != v) continue;
      for (Uint c = 0; c < dim; ++c) u_max = std::max(u_max, std::abs(u_data[v*dim + c]));
    }
  }
  else {
    for (const double u : u_data) u_max = std::max(u_max, std::abs(u));
  }
  return 1e-12*u_max;
}

#endif
