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
#include "wall_facets.hpp"
#include <array>
#include <algorithm>
#include <atomic>
#include <cassert>
#include <iostream>
#include <type_traits>

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
    // The no-slip walls: wall_facets', else the exterior facets off the periodic sides
    std::vector<std::uint8_t> walls;
    if constexpr (Format::named_walls){
      walls_given_ = dolfin_params.has("wall_facets");
      if (walls_given_){
        walls = wall_facets::load(get_folder(), dolfin_params.template get<std::string>("wall_facets"),
                                  dolfin_params.template get<std::string>("wall_markers"),
                                  topo_, coords_, ncells_, dim, facet_neigh_, verbose);
        near_wall::label_cell_type(cell_type_, facet_neigh_, Cell::n_verts, walls);
      }
    }
    near_wall::build_tables<Cell>(topo_, coords_, ncells_, nverts_, dim, periodic, x_min, x_max,
                                  periodic_tol(), facet_neigh_, walls, upgrade_p2_, vclass_, wall_index_,
                                  wall_cells_, split_.cell, verbose);
    if constexpr (Format::named_walls)
      if (!walls_given_) near_wall::warn_velocity_walls(wall_index_);
    // wall_facets' walls are at rest
    if (walls_given_){
      if (upgrade_p2_)
        near_wall::build_wall_terms<Cell>(topo_, coords_, ncells_, nverts_, dim, facet_neigh_, walls, vclass_,
                                          wall_index_, split_, verbose);
      rest_fixed_ = true;
    }
    wall_eval_ = rest_fixed_ || rest_by_stamp;
    if (!stamps_.a().u.empty()){
      prepare(stamps_.a());
      keep_rest(stamps_.a(), "the first stamp");
      check_walls(stamps_.a());
    }
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
                            if constexpr (std::is_convertible_v<typename Format::Key, std::string>) keep_rest(s, key);
                            else keep_rest(s, fmt_.name(key));
                            check_walls(s);
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
  if (!wall_upgrade_ || walls_given_) return;
  if constexpr (Format::vertex_fields)
    s.rest_tol = near_wall::rest_tolerance_by_vertex(s.u, vclass_, nverts_, dim);
  else
    s.rest_tol = near_wall::rest_tolerance(s.u);
}

template<typename Cell, typename Format>
void StampedInterpol<Cell, Format>::keep_rest(const Stamp& s, const std::string& name)
{
  using S = wall_split::Shape<Cell::n_verts>;
  constexpr unsigned facets = (1u << Cell::n_verts) - 1;
  constexpr unsigned others = S::pattern_mask & ~facets;
  if (!wall_upgrade_ || walls_given_ || rest_by_stamp) return;
  const bool first = !rest_fixed_;
  // No flow: a later stamp fixes the walls
  if (first && s.rest_tol == 0.) return;
  std::size_t moved_facets = 0, moved_others = 0, moved_vertices = 0;
  if (upgrade_p2_){
    if constexpr (!Format::vertex_fields){
#pragma omp parallel for schedule(static) reduction(+:moved_facets, moved_others, moved_vertices)
      for (std::size_t i = 0; i < ncells_; ++i){
        const std::int32_t w = wall_index_[i];
        if (w < 0) continue;
        std::array<double, Cell::n_dofs_max*3> u;
        gather_node_values_n<Cell::n_dofs_max, D>(u_dofs_[i], s.u.data(), u.data());
        const unsigned held = split_.cell[std::size_t(w)].rest;
        const unsigned rest = wall_split::rest_entities<Cell>(u.data(), held, s.rest_tol);
        if (first) split_.cell[std::size_t(w)].rest = std::uint16_t(rest);
        else if (rest != held){
          if ((rest ^ held) & facets) ++moved_facets;
          else if ((rest ^ held) & others) ++moved_others;
          else ++moved_vertices;
        }
      }
      if (first)
        near_wall::build_wall_terms<Cell>(topo_, coords_, ncells_, nverts_, dim, facet_neigh_, {}, vclass_,
                                          wall_index_, split_, verbose);
    }
  }
  else {
    // P1: the wall vertices at rest; a cell with none reads as without the rule
#pragma omp parallel for schedule(static) reduction(+:moved_vertices)
    for (std::size_t i = 0; i < ncells_; ++i){
      const std::int32_t w = wall_index_[i];
      if (w < 0) continue;
      std::array<double, Cell::n_verts*3> u;
      gather_node_values_n<Cell::n_verts, D>(u_dofs_[i], s.u.data(), u.data());
      std::uint8_t& held = wall_cells_[std::size_t(w)].wall;
      const unsigned rest = near_wall::rest_vertices<Cell>(u.data(), Uint(Cell::n_verts), held, s.rest_tol);
      if (first){
        held = std::uint8_t(rest);
        if (!rest) wall_index_[i] = -1;
      }
      else if (rest != held) ++moved_vertices;
    }
  }
  rest_fixed_ = true;
  wall_eval_ = true;
  if (moved_facets || moved_others || moved_vertices){
    const char* others_name = D == 3 ? "wall edges" : "wall vertices";
    std::string what;
    const auto add = [&](const std::size_t n, const std::string& kind){
      if (!n) return;
      what += std::string(what.empty() ? "" : "; ") + kind + " at rest in the first stamp with flow move in "
              "this one, in " + std::to_string(n) + " cells";
    };
    add(moved_facets, "wall facets");
    add(moved_others, others_name);
    add(moved_vertices, "wall vertices");
    partrac::fail(name, ": ", what, "; wall_upgrade = true needs the walls at rest in the first stamp with flow "
                  "at rest in every stamp, or wall_facets naming the no-slip walls");
  }
}

template<typename Cell, typename Format>
void StampedInterpol<Cell, Format>::check_walls(const Stamp& s)
{
  if (!walls_given_ || walls_checked_) return;
  // The largest |u| by vertex where the fields are
  static const std::vector<std::uint32_t> by_node;
  walls_checked_ = near_wall::check_moving_walls<Cell>(u_dofs_, wall_index_, upgrade_p2_, split_.cell, wall_cells_,
                                                       {&s.u}, Format::vertex_fields ? vclass_ : by_node,
                                                       "the first stamp with flow");
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
  static std::atomic<std::uint64_t> count{0};
  stamps_in_play_ = ++count;
}

#endif
