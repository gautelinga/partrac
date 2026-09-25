#ifndef __STAMPEDINTERPOL_HPP
#define __STAMPEDINTERPOL_HPP

// Velocity, pressure and phase fields as stamps in time on triangles or tets,
// blended between the two stamps around t. Format reads the files -- the mesh,
// its node tables and one stamp into a buffer (DolfinH5Format in
// SimplexInterpol.hpp, XDMFFormat in XDMFInterpol.hpp, OpenFoamFormat in
// OpenFoamInterpol.hpp); the evaluation, the near-wall rule (near_wall.hpp)
// and the stamp buffer are the same for every format. Nothing here needs
// dolfin.
//
// A Format is a friend that fills the interpolator's tables: it gives Key (what
// names a stamp), vertex_fields (every field on the velocity's vertex table),
// schema(D), the timestamps ts, load (the mesh, its tables and any first stamp),
// read (one stamp into a buffer), key and name (a timestamp's Key and its name
// in the log). Its schema declares include_phi and wall_p2, which are read here.

#include <cstdint>
#include <string>
#include <vector>

#include "MeshCore.hpp"
#include "near_wall.hpp"
#include "stamp_buffer.hpp"
#include "Triangle.hpp"
#include "Tet.hpp"
#include <omp.h>

template<typename Cell, typename Format>
class StampedInterpol final
  : public MeshCore<Cell>
{
public:
  StampedInterpol(const std::string& infilename);
  void update(const double t);
  // One stamp for every t: the blend at t
  void freeze(const double t);
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  // What a step reads: velocity, acceleration and their gradients; P, Phi, cell_type stay zero
  void evaluate_motion(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  double get_t_min() { return fmt_.ts.get_t_min(); };
  double get_t_max() { return fmt_.ts.get_t_max(); };
  // The two stamps are the same values, not a copy of them
  bool stamps_aliased() const { return stamps_.aliased(); }
  bool has_phase_field() const override { return include_phi; }
  bool has_phase_gradient() const override { return phase_gradient_; }
  void evaluate_phase_gradient(const Vector3d &x, const double t, const CellPos& pos, Vector3d& g) override;
protected:
  friend Format;
  template<bool Scalars>
  void evaluate_impl(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  static constexpr int D = Cell::n_verts - 1;
  using Base = MeshCore<Cell>;
  // Names the base owns
  using Base::dolfin_params; using Base::periodic; using Base::include_pressure;
  using Base::periodic_tol; using Base::dim; using Base::hmin_;
  using Base::ncoeffs_u; using Base::ncoeffs_p;
  using Base::cells_; using Base::facet_neigh_; using Base::period_;
  using Base::u_dofs_; using Base::p_dofs_; using Base::t_prev; using Base::t_next;
  using Base::x_min; using Base::x_max; using Base::is_initialized; using Base::t_update;
  using Base::get_folder; using Base::set_folder; using Base::wants_gradient;
  using Base::verbose;
  using Base::read_mesh_params; using Base::set_period; using Base::tree_;

  // A stamp's fields by node, and the near-wall rule's tolerance computed from its velocity
  struct Stamp {
    std::vector<double> u, p, phi;   // phi: nverts_ values, then with phase_gradient_ D a vertex
    double rest_tol = 0.;   // |u| at rest on a wall
  };

  // The node tables of the scalar fields; a format that reads every field by
  // vertex has the velocity's alone
  const CellDofs& p_table() const {
    if constexpr (Format::vertex_fields) return u_dofs_; else return p_dofs_;
  }
  const CellDofs& phi_table() const {
    if constexpr (Format::vertex_fields) return u_dofs_; else return phi_dofs_;
  }

  Format fmt_;

  bool include_phi = false;
  Uint ncoeffs_phi = 0;
  CellDofs phi_dofs_;

  // The stamps the fields are blended between, keyed by the format's name for a stamp
  partrac::StampBuffer<Stamp, typename Format::Key> stamps_;
  const double* u_prev_ = nullptr;
  const double* u_next_ = nullptr;
  const double* p_prev_ = nullptr;
  const double* p_next_ = nullptr;
  const double* phi_prev_ = nullptr;
  const double* phi_next_ = nullptr;

  std::vector<int> cell_type_;

  enum class WallP2 { Edge, None };
  WallP2 wall_p2_ = WallP2::None;

  using WallEdges = near_wall::WallEdges<Cell>;
  std::vector<std::int32_t> wall_index_;   // -1: no wall vertex
  std::vector<WallEdges> wall_cells_;

  // Velocity, acceleration and their gradients in a wall cell, in a unit of
  // its own: P2 for a stamp whose wall vertices are at rest, else P1
  void wall_motion(const int id, const CellPos& pos, const double alpha_t, PointValues& fields) const;

  double rest_tol_prev_ = 0.;   // the tolerance of the stamps in play
  double rest_tol_next_ = 0.;
  double rest_tol(const std::vector<double>& u_data) const;

  // The arrays the tree addresses; the cells hold the geometry
  std::vector<std::uint32_t> topo_;
  std::vector<double> coords_;
  std::size_t ncells_ = 0, nverts_ = 0;
  // A vertex's master across the periodic faces, its own id where there is none
  std::vector<std::uint32_t> vclass_;

  Stamp frozen_;   // a blend frozen between two stamps

  // The phase field's own gradient after its values, on the velocity's vertex table
  bool phase_gradient_ = false;
};

#endif
