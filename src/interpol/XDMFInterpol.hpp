#ifndef __XDMFINTERPOL_HPP
#define __XDMFINTERPOL_HPP

// P1 velocity, pressure and phase fields written as XDMF, on triangles or
// tets. The dofs of a P1 field are the mesh vertices, so the values are read
// by vertex and every table comes from the mesh arrays the XDMF names; nothing
// here needs dolfin. The near-wall rule (near_wall.hpp) differs by dimension,
// so wall_block is specialised, and it lives here, where it inlines into the
// evaluation.

#include <memory>
#include <string>
#include <vector>

#include "MeshCore.hpp"
#include "cell_tree.hpp"
#include "near_wall.hpp"
#include "strings.hpp"
#include "Timestamps.hpp"
#include "Triangle.hpp"
#include "Tet.hpp"
#include <omp.h>

template<typename Cell>
class XDMFInterpol final
  : public MeshCore<Cell>
{
public:
  XDMFInterpol(const std::string& infilename);
  void update(const double t);
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& );
  // What a step reads: velocity, acceleration and their gradients; P, Phi, cell_type stay zero
  void evaluate_motion(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  double get_t_min() { return ts.get_t_min(); };
  double get_t_max() { return ts.get_t_max(); };
protected:
  template<bool Scalars>
  void evaluate_impl(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  static constexpr int D = Cell::n_verts - 1;
  using Base = MeshCore<Cell>;
  // Names the base owns
  using Base::dolfin_params; using Base::periodic; using Base::include_pressure;
  using Base::periodic_tol; using Base::dim; using Base::hmin_;
  using Base::ncoeffs_u; using Base::ncoeffs_p;
  using Base::cells_; using Base::facet_neigh_; using Base::period_;
  using Base::u_dofs_; using Base::t_prev; using Base::t_next;
  using Base::x_min; using Base::x_max; using Base::is_initialized; using Base::t_update;
  using Base::get_folder; using Base::set_folder; using Base::wants_gradient;
  using Base::verbose;
  using Base::read_mesh_params; using Base::set_period;
  // Out of the step loops: the cell tree, from nothing known
  bool locate_tree(const Vector3d& xx, CellPos& pos) override;
  // One stamp of one field by vertex, the first ncols of each stored row
  void read_stamp(const std::vector<std::string>& path, std::vector<double>& data, const int ncols);

  MultiTimestamps ts;

  bool include_phi = true;

  // A field's values by vertex; the pressure and the phase field share the
  // velocity's node table, since every P1 dof is a vertex
  std::vector<double> u_prev_data_;
  std::vector<double> u_next_data_;
  std::vector<double> p_prev_data_;
  std::vector<double> p_next_data_;
  std::vector<double> phi_prev_data_;
  std::vector<double> phi_next_data_;

  std::vector<int> cell_type_;

  enum class WallP2 { Edge, None };
  WallP2 wall_p2_ = WallP2::Edge;

  using WallEdges = near_wall::WallEdges<Cell>;
  std::vector<std::int32_t> wall_index_;   // -1: no wall vertex
  std::vector<WallEdges> wall_cells_;

  // The near-wall rule, per dimension
  bool wall_block(const double* u, double* u2, const WallEdges& w, const double tol) const;

  double rest_tol_prev_ = 0.;   // |u| at rest on a wall
  double rest_tol_next_ = 0.;
  double rest_tol(const std::vector<double>& u_data) const;

  // The arrays the tree addresses; the cells hold the geometry
  std::vector<std::uint32_t> topo_;
  std::vector<double> coords_;
  std::size_t ncells_ = 0, nverts_ = 0;
  // A vertex's master across the periodic faces, its own id where there is none
  std::vector<std::uint32_t> vclass_;
  std::unique_ptr<partrac::CellTree> tree_;
};

#endif
