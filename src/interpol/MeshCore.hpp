#ifndef __MESHCORE_HPP
#define __MESHCORE_HPP

// What an interpolator on a simplex mesh does once its tables exist: hold the
// cells and the facet table, locate a point, and walk a move off the walls.
// None of it depends on where the tables came from, so it compiles without
// dolfin; a loader fills them, DolfInterpol from a dolfin mesh and the rest
// from the file's own arrays.

#include <cassert>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include "Interpol.hpp"
#include "Params.hpp"
#include "cell_tree.hpp"
#include "cell_walk.hpp"
#include "geometry.hpp"
#include "simplex_load.hpp"

template<typename Cell>
class MeshCore : public Interpol {
public:
  MeshCore(const std::string& infilename) : Interpol(infilename) {}
  bool locate(const Vector3d &x, const double t, CellPos& pos){
    assert(t <= t_next && t >= t_prev);
    const Vector3d xx = _modx(x);
    return walk_to_cell(cells_, facet_neigh_, xx, pos) || locate_tree(xx, pos);
  }
  bool reflect(const Vector3d& x, Vector3d& dx, CellPos& pos);
  void enable_reflection();
  // Outward unit normal of the cell's wall facets, their mean at an edge; zero off the wall
  Vector3d get_boundary_normal(const Vector3d &x, int& cell_id);
  double hmin() const { return hmin_; }
  using Interpol::locate;
  using Interpol::evaluate;
protected:
  // Periodicity and the pressure flag from the parameter file
  void read_mesh_params();
  // What build_tables produced, into the names the base owns; the phase field,
  // the mesh arrays and the cell counts are the loader's own
  void adopt_tables(simplex_load::Tables& t){
    dim = t.mesh.gdim;
    x_min = t.mesh.x_min;
    x_max = t.mesh.x_max;
    set_period();
    ncoeffs_u = t.ncoeffs_u;
    ncoeffs_p = t.ncoeffs_p;
    u_dofs_ = std::move(t.u_dofs);
    p_dofs_ = std::move(t.p_dofs);
    facet_neigh_ = std::move(t.facets);
    hmin_ = t.hmin;
  }
  // The box lengths along the periodic axes, once x_min and x_max are final
  void set_period();
  // Out of the step loops: the cell tree, from nothing known
  bool locate_tree(const Vector3d& xx, CellPos& pos);
  // Into the box along the periodic axes
  Vector3d _modx(const Vector3d& x) const {
    Vector3d x_loc = x;
    for (Uint i = 0; i < dim; ++i)
      if (period_[i] > 0. && (x[i] < x_min[i] || x[i] >= x_max[i]))
        x_loc[i] = x_min[i] + modulox(x[i] - x_min[i], period_[i]);
    return x_loc;
  }

  partrac::Params dolfin_params;
  double t_prev = 0.;   // the stamps the fields are between
  double t_next = 0.;
  std::vector<bool> periodic = {false, false, false};
  bool include_pressure = true;
  double periodic_tol = 1e-12;   // heuristic

  Uint dim;
  double hmin_ = 0.;    // the shortest cell diameter
  Uint ncoeffs_u;
  Uint ncoeffs_p = 0;   // stays 0 when pressure is ignored

  std::vector<Cell> cells_;
  // The cells from nothing known; it addresses the loader's own arrays
  std::unique_ptr<partrac::CellTree> tree_;
  std::vector<std::int32_t> facet_neigh_;     // walk_to_cell, reflect_in_cells
  Vector3d period_ = Vector3d::Zero();
  CellDofs u_dofs_;
  CellDofs p_dofs_;
};

#endif
