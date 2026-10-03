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
#include "Timestamps.hpp"
#include "Triangle.hpp"
#include "Tet.hpp"

// A sub-cell of the split: the macro cell's region and the macro vertex it leaves out
struct SubRegion : Region {
  int sub = 0;
};

template<typename Cell>
class SplitInterpol final
  : public MeshCore<Cell>
{
public:
  SplitInterpol(const std::string& infilename);
  void update(const double t);
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  // What a step reads: velocity, acceleration and their gradients; P, Phi, cell_type stay zero
  void evaluate_motion(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  double get_t_min() { return ts_.get_t_min(); };
  double get_t_max() { return ts_.get_t_max(); };
  double next_stamp_after(const double t) const override { return ts_.next_after(t); }
  // The two stamps are the same values, not a copy of them
  bool stamps_aliased() const { return stamps_.aliased(); }
  bool has_phase_field() const override { return include_phi; }
  // Regions are sub-cells; levels the sub-cell's barycentrics, the last one the macro facet's
  using region_type = SubRegion;
  using typename MeshCore<Cell>::levels_type;
  bool region_of(const int id, const Vector3d& x, const double band, SubRegion& R, levels_type& lev){
    levels_type b;
    if (!Base::region_of(id, x, band, R, b))
      return false;
    R.sub = split_eval::sub_cell<Cell::n_verts>(b, lev.data());
    return true;
  }
  // The sub-cell of the smallest barycentric
  SubRegion region_in(const int id, const Vector3d& x) const {
    SubRegion R;
    static_cast<Region&>(R) = Base::region_in(id, x);
    levels_type b, mu;
    Base::region_point(R, x, b);
    R.sub = split_eval::sub_cell<Cell::n_verts>(b, mu.data());
    return R;
  }
  Vector3d region_point(const SubRegion& R, const Vector3d& x, levels_type& lev) const {
    levels_type b;
    const Vector3d xe = Base::region_point(R, x, b);
    split_eval::sub_levels<Cell::n_verts>(R.sub, b, lev.data());
    return xe;
  }
  // An internal plane to the sub-cell of the same cell, the macro facet to the neighbour's sub-cell on it
  int across(const SubRegion& R, const int k, const Vector3d& x, SubRegion& next) const {
    constexpr int nv = Cell::n_verts;
    if (k < nv - 1){
      next = R;
      next.sub = split_eval::sub_beyond(R.sub, k);
      return split_eval::sub_plane(next.sub, R.sub);
    }
    Region macro;
    const int e = Base::across(R, R.sub, x, macro);
    if (e < 0)
      return -1;
    static_cast<Region&>(next) = macro;
    next.sub = e;
    return nv - 1;
  }
  void level_rates(const SubRegion& R, const Vector3d& v, levels_type& rate) const {
    constexpr int nv = Cell::n_verts;
    levels_type d;
    Base::level_rates(R, v, d);
    for (int m = 0; m < nv - 1; ++m) rate[m] = d[split_eval::sub_beyond(R.sub, m)] - d[R.sub];
    rate[nv - 1] = double(nv)*d[R.sub];
  }
  Vector3d level_grad(const SubRegion& R, const int k) const {
    constexpr int nv = Cell::n_verts;
    const Vector3d gi = Base::level_grad(R, R.sub);
    return k < nv - 1 ? Vector3d(Base::level_grad(R, split_eval::sub_beyond(R.sub, k)) - gi) : Vector3d(double(nv)*gi);
  }
  bool is_wall(const SubRegion& R, const int k) const {
    return k == Cell::n_verts - 1 && Base::is_wall(R, R.sub);
  }
  // A third of the cell's: a path crosses its sub-cells' planes too
  double region_size(const int id) const { return Base::cell_size(id)/3.; }
  // A sub-cell's velocity nodes at both stamps, gathered once for the evaluations in it
  struct Held {
    std::array<double, split_eval::Split<Cell>::n_sub*3> prev, next;
    int sub;
  };
  void hold(const SubRegion& R, Held& h) const;
  // evaluate_motion at the barycentrics lev of the held sub-cell of cell id
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
