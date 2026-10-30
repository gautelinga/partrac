#ifndef __STAMPEDINTERPOL_IMPL_HPP
#define __STAMPEDINTERPOL_IMPL_HPP

// StampedInterpol's loading and update, included by each format's translation
// unit and instantiated there, two cells a unit; the evaluation is in
// StampedInterpol_eval.hpp

#include "Error.hpp"
#include "StampedInterpol.hpp"
#include "mesh_tables.hpp"
#include "near_wall.hpp"
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

  upgrade_p2_ = wall_upgrade_ && near_wall::velocity_p2<Cell>(ncoeffs_u);

  // Identify edge cells
  near_wall::label_cell_type(cell_type_, facet_neigh_, Cell::n_verts);
  std::cout << "Built neighbour list" << std::endl;

  if (wall_upgrade_){
    near_wall::build_tables<Cell>(topo_, coords_, ncells_, nverts_, dim, periodic, x_min, x_max,
                                  periodic_tol(), facet_neigh_, upgrade_p2_, vclass_, wall_index_,
                                  wall_cells_, upgrade_, verbose);
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
  partrac::update_bracket(stamps_, is_initialized, t_prev, t_next, sp.prev, sp.next,
                          [&](const auto& end){ return fmt_.key(*this, end); },
                          [&](const auto& end){ return fmt_.name(end); },
                          [&](const typename Format::Key& key, Stamp& s){
                            fmt_.read(*this, key, s);
                            prepare(s);
                          },
                          [&](const Stamp& prev, const Stamp& next){ point(prev, next); });
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
  if constexpr (Format::vertex_fields)
    s.rest_tol = near_wall::rest_tolerance_by_vertex(s.u, vclass_, nverts_, dim);
  else {
    s.rest_tol = near_wall::rest_tolerance(s.u);
    if (upgrade_p2_) near_wall::field_terms<Cell>(upgrade_, u_dofs_, s.u.data(), s.rest_tol, s.coef);
  }
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

#endif
