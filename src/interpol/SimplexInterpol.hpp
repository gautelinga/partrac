#ifndef __SIMPLEXINTERPOL_HPP
#define __SIMPLEXINTERPOL_HPP

// Velocity and pressure stamps written as dolfin HDF5 checkpoints, on
// triangles or tets: a P1 or P2 Taylor-Hood pair read from the stored dofmap
// and the mesh arrays, and evaluated through this code's own cell walk, cell
// tree and basis. Nothing here needs dolfin.

#include <memory>
#include <string>
#include <vector>

#include "MeshCore.hpp"
#include "cell_tree.hpp"
#include "mesh_tables.hpp"
#include "strings.hpp"
#include "Timestamps.hpp"
#include "Triangle.hpp"
#include "Tet.hpp"
#include <omp.h>

template<typename Cell>
class SimplexInterpol final
  : public MeshCore<Cell>
{
public:
  SimplexInterpol(const std::string& infilename);
  void update(const double t);
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& );
  // What a step reads: velocity, acceleration and their gradients; P, Phi, cell_type stay zero
  void evaluate_motion(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  double get_t_min() { return ts.get_t_min(); };
  double get_t_max() { return ts.get_t_max(); };
  // The two stamps are the same values, not a copy of them
  bool stamps_aliased() const { return u_prev_ == u_next_; }
protected:
  template<bool Scalars>
  void evaluate_impl(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  static constexpr int D = Cell::n_verts - 1;
  static constexpr const char* mode = D == 2 ? "triangle" : "tet";
  using Base = MeshCore<Cell>;
  // Names the base owns
  using Base::dolfin_params; using Base::periodic; using Base::include_pressure;
  using Base::periodic_tol; using Base::dim; using Base::hmin_;
  using Base::ncoeffs_u; using Base::ncoeffs_p;
  using Base::cells_; using Base::facet_neigh_; using Base::period_;
  using Base::u_dofs_; using Base::p_dofs_; using Base::t_prev; using Base::t_next;
  using Base::x_min; using Base::x_max; using Base::is_initialized; using Base::t_update;
  using Base::get_folder; using Base::set_folder; using Base::wants_gradient;
  using Base::read_mesh_params; using Base::set_period;
  // Out of the step loops: the cell tree, from nothing known
  bool locate_tree(const Vector3d& xx, CellPos& pos) override;
  // One stamp's vectors into the buffers, through the retained dof mappings
  void read_stamp(const std::string& filename, std::vector<double>& u_buf, std::vector<double>& p_buf);

  Timestamps ts;

  // The arrays the tree addresses; the cells hold the geometry
  std::vector<std::uint32_t> topo_;
  std::vector<double> coords_;
  std::size_t ncells_ = 0, nverts_ = 0;
  std::unique_ptr<partrac::CellTree> tree_;

  // A stamp's values by node, D (or one) doubles a node; a single-stamp field
  // aliases both pointers onto the same buffer
  std::vector<double> u_a_, u_b_, p_a_, p_b_;
  std::string file_a_, file_b_;
  const double* u_prev_ = nullptr;
  const double* u_next_ = nullptr;
  const double* p_prev_ = nullptr;
  const double* p_next_ = nullptr;
  // Stored dof -> every node slot it feeds, for the later stamps
  mesh_tables::DofNodes u_map_, p_map_;
  std::string u_field_, p_field_;
};

#endif
