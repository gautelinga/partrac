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
// in the log). Its schema declares include_phi, wall_upgrade, wall_facets and
// wall_markers, which are read here; named_walls: accepts wall_facets.

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "MeshCore.hpp"
#include "near_wall.hpp"
#include "wall_split.hpp"
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
  double get_t_min() const override { return fmt_.ts.get_t_min(); };
  double get_t_max() const override { return fmt_.ts.get_t_max(); };
  double next_stamp_after(const double t) const override { return fmt_.ts.next_after(t); }
  // The two stamps are the same values, not a copy of them
  bool stamps_aliased() const { return stamps_.aliased(); }
  bool has_phase_field() const override { return include_phi; }
  bool has_phase_gradient() const override { return phase_gradient_; }
  void evaluate_phase_gradient(const Vector3d &x, const double t, const CellPos& pos, Vector3d& g) override;
  // A cell's velocity nodes at both stamps, gathered once for the evaluations in it
  struct Held {
    std::array<double, Cell::n_dofs_max*3> prev, next;
    std::array<double, (Cell::n_verts - 1)*Cell::n_dofs_max> prev2, next2;   // a wall cell's P2 blocks
    wall_split::Ready<Cell> ready;   // a P2 wall cell's, of both stamps
    bool wall, at_rest_prev, at_rest_next;
  };
  void hold(const Region& R, Held& h) const;
  // evaluate_motion at the barycentrics of cell id from its held nodes
  void held_motion(const int id, const std::array<double, 4>& lev, const double t, const Held& h,
                   PointValues& fields);
  // The same without the gradients
  void held_velocity(const int id, const std::array<double, 4>& lev, const double t, const Held& h,
                     PointValues& fields);
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

  // A stamp's fields by node, and the tolerance its walls are found at rest by
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

  // Near walls: P1 to P2 by the edge rule, or P2 plus the wall correction
  bool wall_upgrade_ = false;
  bool upgrade_p2_ = false;   // never read where every field is by vertex
  bool wall_eval_ = false;    // the wall cells read as such: once their rest is fixed
  bool walls_given_ = false;  // wall_facets: its walls, all at rest
  bool walls_checked_ = false;   // against the first stamp with flow

  using WallEdges = near_wall::WallEdges<Cell>;
  std::vector<std::int32_t> wall_index_;   // -1: no wall vertex or entity at rest
  // P1: the edge ends; wall: the wall vertices, then those at rest
  std::vector<WallEdges> wall_cells_;
  // P2: the wall patterns at rest and their tables
  wall_split::Tables<Cell::n_verts> split_;
  bool rest_fixed_ = false;   // by wall_facets or the first stamp with flow
  // Without wall_facets to name the walls (OpenFOAM, whose patches may start
  // moving): a P1 stamp's wall vertices at rest are its own
  static constexpr bool rest_by_stamp = !Format::named_walls;
  double rest_tol_prev_ = 0.;   // rest_by_stamp: the tolerance of the stamps in play
  double rest_tol_next_ = 0.;
  // A thread's last located P2 wall cell, ready while the stamps in play stay
  using WallMemo = wall_split::Memo<Cell, std::uint64_t>;
  static WallMemo& wall_memo(){ static thread_local WallMemo m; return m; }
  std::uint64_t stamps_in_play_ = 0;   // new with every point
  // What a stamp read needs computed: the rest tolerance, without wall_facets
  void prepare(Stamp& s) const;
  // Without wall_facets: the wall vertices (P1) or entities (P2) at rest in the
  // first stamp with flow, which every later one must keep
  void keep_rest(const Stamp& s, const std::string& name);
  // With wall_facets: warns once of its walls' velocity in the first stamp with flow
  void check_walls(const Stamp& s);
  // The stamps in play
  void point(const Stamp& prev, const Stamp& next);

  // A wall cell, out of line: gathers, then wall_at
  __attribute__((noinline))
  void wall_motion(const int id, const CellPos& pos, const double alpha_t, PointValues& fields) const;
  // The per-point part of the evaluations, from gathered nodes; Velocity: no gradients
  template<bool Velocity = false>
  void plain_at(const int id, const std::array<double, 4>& bary, const double alpha_t,
                const double* prev, const double* next, PointValues& fields) const;
  // A P1 wall cell's P2 blocks, for the stamps whose wall vertices are at rest
  std::array<bool, 2> wall_blocks(const int id, const double* prev, const double* next, double* prev2,
                                  double* next2) const;
  template<bool Velocity = false>
  void wall_at(const int id, const std::array<double, 4>& bary, const double alpha_t, const double* prev,
               const double* next, const double* prev2, const double* next2, const bool at_rest_prev,
               const bool at_rest_next, PointValues& fields) const;
  // A P2 wall cell, out of line: the P2 field plus the wall correction at each
  // stamp; Hold: from what h holds, else from a thread's last cell made ready
  template<bool Velocity = false, bool Hold = false>
  __attribute__((noinline))
  void upgrade_at(const int id, const std::array<double, 4>& bary, const double alpha_t, const Held* h,
                  PointValues& fields) const;
  // What a P2 wall cell's held evaluation reads, from h's blocks, out of line
  __attribute__((noinline))
  void wall_hold(const int id, Held& h) const;
  template<bool Velocity = false>
  void held_at(const int id, const std::array<double, 4>& lev, const double t, const Held& h,
               PointValues& fields) const;

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
