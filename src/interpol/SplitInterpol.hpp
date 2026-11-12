#ifndef __SPLITINTERPOL_HPP
#define __SPLITINTERPOL_HPP

// Velocity, pressure and phase stamps written as dolfin HDF5 checkpoints, read
// as SimplexInterpol reads them, but with the velocity evaluated on each cell's
// barycentric (Alfeld) split, where it is pointwise divergence-free
// (split_eval.hpp). The file's every cell must already have zero net flux,
// which python/divfree/divfree_clean.py prepares; the interior values are recomputed
// from the P2 boundary data at every stamp and never read. Nothing here needs
// dolfin.

#include <cstdint>
#include <string>
#include <vector>

#include "MeshCore.hpp"
#include "mesh_tables.hpp"
#include "split_eval.hpp"
#include "stamp_buffer.hpp"
#include "SubCellRegions.hpp"
#include "Timestamps.hpp"
#include "Triangle.hpp"
#include "Tet.hpp"

template<typename Cell>
class SplitInterpol final
  : public SubCellRegions<Cell, SplitInterpol<Cell>, true>
{
public:
  SplitInterpol(const std::string& infilename);
  void update(const double t);
  // The two stamps blended at t, held for every later time
  void freeze(const double t) override;
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  // What a step reads: velocity, acceleration and their gradients; P, Phi, cell_type stay zero
  void evaluate_motion(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  double get_t_min() const override { return ts_.get_t_min(); };
  double get_t_max() const override { return ts_.get_t_max(); };
  double next_stamp_after(const double t) const override { return ts_.next_after(t); }
  // The two stamps are the same values, not a copy of them
  bool stamps_aliased() const { return stamps_.aliased(); }
  bool has_phase_field() const override { return include_phi; }
  // A sub-cell's velocity nodes at both stamps, gathered once for the evaluations in it
  struct Held {
    std::array<double, split_eval::Split<Cell>::n_sub*3> prev, next;
    int sub;
  };
  void hold(const SubRegion& R, Held& h) const;
  // evaluate_motion at the levels lev of the held sub-cell of cell id
  void held_motion(const int id, const std::array<double, 4>& lev, const double t, const Held& h,
                   PointValues& fields);
  // The same without the gradients
  void held_velocity(const int id, const std::array<double, 4>& lev, const double t, const Held& h,
                     PointValues& fields);
protected:
  template<bool Scalars>
  void evaluate_impl(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  // Velocity and its rate in sub-cell i of cell id at mu, from the sub-cell's gathered nodes; Velocity: no gradients
  template<bool Velocity = false>
  void sub_at(const int id, const int i, const double* mu, const double alpha_t, const double* prev,
              const double* next, PointValues& fields) const;
  template<bool Velocity = false>
  void held_at(const int id, const std::array<double, 4>& lev, const double t, const Held& h,
               PointValues& fields) const;
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
  using Base::read_mesh_params; using Base::set_period; using Base::tree_;

  using Split = split_eval::Split<Cell>;

  // A stamp's fields by node, and the interior values of every cell's split
  struct Stamp {
    std::vector<double> u, p, phi;
    std::vector<double> u_int;   // ncells blocks of n_int nodes, D doubles a node
  };

  // One stamp's vectors, then its interior values; the reference the failures name
  void read_stamp(const std::string& file, Stamp& s);
  // Point the evaluation at two stamps
  void point(const Stamp& prev, const Stamp& next);
  // Every vertex of every exterior facet at rest in the stamp
  bool walls_at_rest(const Stamp& s) const;
  // Every cell's net flux against its scale, floored at the stamp's, and the
  // interior values through the reference matrix; cells in parallel and independent
  void build_interior(const std::string& file, Stamp& s);

  Timestamps ts_;
  std::string u_field, p_field, phi_field;
  // Stored dof -> every node slot it feeds, for the later stamps
  mesh_tables::DofNodes u_map_, p_map_, phi_map_;

  bool include_phi = false;
  Uint ncoeffs_phi = 0;
  CellDofs phi_dofs_;

  typename Split::RefMatrix R_;

  // The stamps the fields are blended between, keyed by the stamp's file
  partrac::StampBuffer<Stamp, std::string> stamps_;
  Stamp frozen_;   // a blend frozen between two stamps
  const double* u_prev_ = nullptr;
  const double* u_next_ = nullptr;
  const double* int_prev_ = nullptr;
  const double* int_next_ = nullptr;
  const double* p_prev_ = nullptr;
  const double* p_next_ = nullptr;
  const double* phi_prev_ = nullptr;
  const double* phi_next_ = nullptr;

  std::vector<int> cell_type_;

  // The arrays the tree addresses; the cells hold the geometry
  std::vector<std::uint32_t> topo_;
  std::vector<double> coords_;
  std::size_t ncells_ = 0, nverts_ = 0;
};

#endif
