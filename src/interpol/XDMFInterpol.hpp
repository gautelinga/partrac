#ifdef USE_DOLFIN
#ifndef __XDMFINTERPOL_HPP
#define __XDMFINTERPOL_HPP

// P1 velocity, pressure and phase fields written as XDMF, on triangles or
// tets. The near-wall rule differs by dimension: the midpoint of an edge from
// a wall vertex is a fixed linear map of the velocity at its fluid end, a
// matrix in 2D and u_v/2 + (u_v.n) q in 3D, so build_wall_edges and
// wall_block are specialised.

#include "MeshInterpol.hpp"
#include "strings.hpp"
#include "Timestamps.hpp"
#include "Triangle.hpp"
#include "Tet.hpp"
#include "cell_locate.hpp"
#include <omp.h>

// What the near-wall rule stores per edge end, and the mesh's name in XDMF
template<typename Cell> struct XDMFCell;

template<> struct XDMFCell<Triangle> {
  static constexpr const char* name = "triangle";
  // Midpoint of a wall-vertex edge: M u_v
  struct WallEnd { double mxx, mxy, myx, myy; };
  // Ends by edge (01, 02, 12), then (first, second)
  static constexpr int n_ends = 6;
};

template<> struct XDMFCell<Tet> {
  static constexpr const char* name = "tetrahedron";
  // Midpoint of a wall-vertex edge: u_v/2 + (u_v.n) q
  struct WallEnd { double qx, qy, qz, nx, ny, nz; };
  // Ends by edge (01, 02, 03, 12, 13, 23), then (first, second)
  static constexpr int n_ends = 12;
};

template<typename Cell>
class XDMFInterpol final
  : public MeshInterpol<Cell>
{
public:
  XDMFInterpol(const std::string& infilename);
  ~XDMFInterpol() { std::cout << "Destructing XDMFInterpol (" << XDMFCell<Cell>::name << ")." << std::endl; };
  void update(const double t);
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& );
  double get_t_min() { return ts.get_t_min(); };
  double get_t_max() { return ts.get_t_max(); };
  // Only the tet loader fills the normals
  Vector3d get_boundary_normal(const Vector3d &x, int & cell_id)
  {
    if constexpr (D == 3) return cell_normal_[cell_id];
    else                  return Vector3d::Zero();
  }
protected:
  static constexpr int D = Cell::n_verts - 1;
  using Base = MeshInterpol<Cell>;
  // Names the base owns
  using Base::dolfin_params; using Base::periodic; using Base::include_pressure;
  using Base::periodic_tol; using Base::mesh; using Base::dim;
  using Base::u_space_; using Base::p_space_; using Base::ncoeffs_u; using Base::ncoeffs_p;
  using Base::cells_; using Base::dolfin_cells_; using Base::cell2cells_;
  using Base::u_dofs_; using Base::p_dofs_; using Base::t_prev; using Base::t_next;
  using Base::x_min; using Base::x_max; using Base::is_initialized; using Base::t_update;
  using Base::get_folder; using Base::set_folder; using Base::wants_gradient;
  using Base::read_mesh_params; using Base::init_mesh_geometry; using Base::build_cells;

  MultiTimestamps ts;

  bool include_phi = true;

  // Assuming phi_space_ == p_space_
  std::vector<double> u_prev_data_;
  std::vector<double> u_next_data_;
  std::vector<double> p_prev_data_;
  std::vector<double> p_next_data_;
  std::vector<double> phi_prev_data_;
  std::vector<double> phi_next_data_;

  std::string h5filename_u;
  std::string h5filename_p;
  std::string h5filename_phi;

  std::vector<Uint> i2j;
  std::vector<Uint> j2i;

  std::vector<int> cell_type_;
  std::vector<Vector3d> cell_normal_;
  std::vector<Vector3d> cell_facet_midpoint_;

  enum class WallP2 { Edge, None };
  WallP2 wall_p2_ = WallP2::Edge;

  using WallEnd = typename XDMFCell<Cell>::WallEnd;
  struct WallEdges {
    std::array<WallEnd, XDMFCell<Cell>::n_ends> ends;
    std::uint8_t wall;   // bit k: vertex k lies on a wall
  };
  std::vector<std::int32_t> wall_index_;   // -1: no wall vertex
  std::vector<WallEdges> wall_cells_;

  // The near-wall rule, per dimension
  void build_wall_edges(const double tol);
  bool wall_block(const double* u, double* u2, const WallEdges& w, const double tol) const;

  double rest_tol_prev_ = 0.;   // |u| at rest on a wall
  double rest_tol_next_ = 0.;
  double rest_tol(const std::vector<double>& u_data) const;
  // The element spaces, per dimension
  void make_spaces(std::shared_ptr<const dolfin::SubDomain> constrained_domain);
};

#endif
#endif
