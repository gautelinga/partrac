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
  // Near walls: the velocity's order raised or not
  wall_upgrade_ = dolfin_params.template get<bool>("wall_upgrade");

  // The mesh, its tables, and the first stamp where the format reads it here
  fmt_.load(*this, infilename);

  upgrade_p2_ = wall_upgrade_ && ncoeffs_u == Uint(Cell::n_dofs_max);
  if (wall_upgrade_ && !upgrade_p2_ && ncoeffs_u != Uint(Cell::n_verts))
    partrac::fail(infilename, ": wall_upgrade = true is for a P1 or P2 velocity, and this one has ",
                  ncoeffs_u, " dofs a cell");

  // Identify edge cells
  near_wall::label_cell_type(cell_type_, facet_neigh_, Cell::n_verts);
  std::cout << "Built neighbour list" << std::endl;

  if (wall_upgrade_){
    if (vclass_.empty())
      vclass_ = mesh_tables::match_periodic_vertices(coords_, nverts_, dim, periodic,
                                                     x_min, x_max, periodic_tol);
    if (upgrade_p2_){
      near_wall::build_wall_terms<Cell>(topo_, coords_, ncells_, nverts_, dim, facet_neigh_,
                                        vclass_, wall_index_, upgrade_, verbose);
    }
    else
      near_wall::build_wall_edges<Cell>(topo_, coords_, ncells_, nverts_, dim, facet_neigh_,
                                        vclass_, wall_index_, wall_cells_, verbose);
    if (!stamps_.a().u.empty()) prepare(stamps_.a());
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

  if (partrac::stamp_reload(is_initialized, {t_prev, t_next}, {sp.prev.t, sp.next.t})){
    partrac::phase_begin("update");
    const auto fill = stamps_.load(fmt_.key(*this, sp.prev), fmt_.key(*this, sp.next),
                                   [&](const typename Format::Key& key, Stamp& s){
      fmt_.read(*this, key, s);
      prepare(s);
    });
    std::cout << "Prev: Timestep = " << sp.prev.t << ", "
              << partrac::stamp_note(fill.first, fmt_.name(sp.prev)) << std::endl;
    std::cout << "Next: Timestep = " << sp.next.t << ", "
              << partrac::stamp_note(fill.second, fmt_.name(sp.next)) << std::endl;

    point(stamps_.prev(), stamps_.next());

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
  update(t);
  const double a = stamp_weight(t, t_prev, t_next);
  if (a != 0.){
    // Blend the two stamps into one
    partrac::blend(a, stamps_.prev().u, stamps_.next().u, frozen_.u);
    partrac::blend(a, stamps_.prev().p, stamps_.next().p, frozen_.p);
    partrac::blend(a, stamps_.prev().phi, stamps_.next().phi, frozen_.phi);
    prepare(frozen_);
    point(frozen_, frozen_);
  }
  else
    point(stamps_.prev(), stamps_.prev());
  // A single stamp: the weight is zero and the rate too at every t
  t_prev = t;
  t_next = t;
  std::cout << "Fields frozen at t = " << t << std::endl;
}

template<typename Cell, typename Format>
void StampedInterpol<Cell, Format>::prepare(Stamp& s) const
{
  if (!wall_upgrade_) return;
  s.rest_tol = rest_tol(s.u);
  if constexpr (!Format::vertex_fields)
    if (upgrade_p2_) stamp_terms(s);
}

template<typename Cell, typename Format>
void StampedInterpol<Cell, Format>::point(const Stamp& prev, const Stamp& next)
{
  u_prev_ = prev.u.data();
  u_next_ = next.u.data();
  p_prev_ = prev.p.data();
  p_next_ = next.p.data();
  phi_prev_ = prev.phi.data();
  phi_next_ = next.phi.data();
  rest_tol_prev_ = prev.rest_tol;
  rest_tol_next_ = next.rest_tol;
  coef_prev_ = prev.coef.data();
  coef_next_ = next.coef.data();
}

template<typename Cell, typename Format>
void StampedInterpol<Cell, Format>::stamp_terms(Stamp& s) const
{
  constexpr int nt = near_wall::n_terms<Cell>;
  const std::size_t ncw = upgrade_.cells.size();
  s.coef.resize(upgrade_.n_coef);
#pragma omp parallel for schedule(dynamic, 1024)
  for (std::size_t w = 0; w < ncw; ++w){
    const near_wall::UpgradeCell& wc = upgrade_.cells[w];
    const near_wall::TermSpan sp = upgrade_.spans[w];
    std::array<double, Cell::n_dofs_max*3> u;
    gather_cell_nodes<Cell::n_verts, Cell::n_dofs_max, D>(u_dofs_[wc.cell], u_dofs_.stride(), s.u.data(), u.data());
    std::array<double, D*nt> c;
    const bool rest = near_wall::upgrade_coefficients<Cell>(u.data(), c.data(), wc, sp.mask,
                                                            upgrade_.normals.data() + wc.normals, s.rest_tol);
    // The cell's terms only, in their order; none at rest: zero
    double* out = s.coef.data() + sp.offset;
    for (std::uint32_t m = sp.mask; m; m &= m - 1, out += D)
      for (int d = 0; d < D; ++d) out[d] = rest ? c[std::size_t(D*__builtin_ctz(m) + d)] : 0.;
  }
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
