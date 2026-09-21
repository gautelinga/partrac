#ifdef USE_DOLFIN
#ifndef __DOLFINTERPOL_HPP
#define __DOLFINTERPOL_HPP

// Velocity and pressure stamps as dolfin HDF5 checkpoints, in any of the
// compiled Lagrange elements (P1-P3), evaluated by dolfin's own basis: the
// reference the other loaders are compared with. The only loader left that
// builds its tables from a dolfin mesh and a dofmap; what a step reads once
// they exist is MeshCore's.

#include <memory>
#include <string>
#include <vector>
#include <dolfin.h>
#include "MeshCore.hpp"
#include "Timestamps.hpp"
#include "Triangle.hpp"
#include "Tet.hpp"

template<typename Cell>
class DolfInterpol final
  : public MeshCore<Cell>
{
public:
  DolfInterpol(const std::string& infilename);
  void update(const double t);
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& ptvals);
  double get_t_min() { return ts.get_t_min(); };
  double get_t_max() { return ts.get_t_max(); };
protected:
  static constexpr int D = Cell::n_verts - 1;
  using Base = MeshCore<Cell>;
  // Names the base owns
  using Base::dolfin_params; using Base::periodic; using Base::periodic_tol;
  using Base::include_pressure; using Base::dim; using Base::hmin_;
  using Base::cells_; using Base::facet_neigh_; using Base::period_;
  using Base::u_dofs_; using Base::p_dofs_;
  using Base::t_prev; using Base::t_next; using Base::x_min; using Base::x_max;
  using Base::is_initialized; using Base::t_update;
  using Base::get_folder; using Base::set_folder; using Base::wants_gradient;
  using Base::read_mesh_params; using Base::set_period;
  using Base::_modx;

  // Dimension, mesh tables, bounding box tree and the domain bounds
  void init_mesh_geometry();
  // Cells in dof order, with their facet table
  void build_cells(const dolfin::GenericDofMap& dofmap);
  // The neighbour across each facet: a cell, a wall or a periodic image
  void build_facet_table();
  // Out of the step loops: the tree, when the walk does not find the cell
  bool locate_tree(const Vector3d& xx, CellPos& pos) override;

  std::shared_ptr<dolfin::Mesh> mesh;
  std::shared_ptr<dolfin::FunctionSpace> u_space_;
  std::shared_ptr<dolfin::FunctionSpace> p_space_;
  std::vector<dolfin::Cell> dolfin_cells_;
  std::vector<std::uint32_t> dolfin2local_;   // empty: dolfin's cell order

  Timestamps ts;

  std::shared_ptr<dolfin::Function> u_prev_;
  std::shared_ptr<dolfin::Function> u_next_;
  std::shared_ptr<dolfin::Function> p_prev_;
  std::shared_ptr<dolfin::Function> p_next_;
  // Per cell, once: what evaluate used to rebuild on every call
  std::vector<int> cell_orientations_;   // empty when the mesh carries none
  std::vector<double> coordinate_dofs_;  // ncoords_ of them per cell, flat
  Uint ncoords_ = 0;
  std::shared_ptr<const dolfin::FiniteElement> u_element_, p_element_;
  Uint u_dim_ = 0, p_dim_ = 0;           // p_dim_ stays 0 when pressure is ignored
  // Read out whole at each load: dolfin's vector is not safe to read in parallel
  std::vector<double> u_prev_data_, u_next_data_, p_prev_data_, p_next_data_;
};

// Geometric dimension of the mesh a fenics parameter file names; 0 if unknown
Uint dolfin_mesh_dim(const std::string& infilename);

#endif
#endif
