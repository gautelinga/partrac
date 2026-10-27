#ifndef __SIMPLEXFREQINTERPOL_HPP
#define __SIMPLEXFREQINTERPOL_HPP

// A time series given as frequency components: one steady Taylor-Hood pair per
// component, each a dolfin HDF5 checkpoint on the same triangle or tet mesh,
// summed with a cosine of the base frequency. The components are read the way
// SimplexInterpol reads a stamp and held by node, so an evaluation gathers the
// cell's nodes out of each component in turn. With wall_upgrade, a wall cell's
// components are summed and the sum raised by the near-wall rule, at the wall
// vertices where every component is at rest.

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "MeshCore.hpp"
#include "FreqStamps.hpp"
#include "near_wall.hpp"
#include "Triangle.hpp"
#include "Tet.hpp"
#include <omp.h>

// Every component's weight at one time: the cosine, and its rate
struct FreqWeights {
  std::uint64_t owner;
  double t;
  const double* w;
  const double* wt;
};

template<typename Cell>
class SimplexFreqInterpol final
  : public MeshCore<Cell>
{
public:
  SimplexFreqInterpol(const std::string& infilename);
  void update(const double t);
  // The components summed at t into one, held for every later time
  void freeze(const double t) override;
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  // What a step reads: velocity, acceleration and their gradients; P, Phi, cell_type stay zero
  void evaluate_motion(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  double get_t_min() { return dolfin_params.template get<double>("t_min"); };
  double get_t_max() { return dolfin_params.template get<double>("t_max"); };
  // A cell's velocity nodes in every component, gathered once for the evaluations in it
  struct Held {
    std::vector<double> blocks;   // a component's block every n_dofs_max*3
    bool wall = false;
  };
  void hold(const Region& R, Held& h) const;
  // evaluate_motion at the barycentrics of cell id from its held nodes
  void held_motion(const int id, const std::array<double, 4>& lev, const double t, const Held& h,
                   PointValues& fields);
protected:
  template<bool Scalars>
  void evaluate_impl(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  static constexpr int D = Cell::n_verts - 1;
  static constexpr const char* mode = D == 2 ? "trianglefreq" : "tetfreq";
  using Base = MeshCore<Cell>;
  // Names the base owns
  using Base::dolfin_params; using Base::periodic; using Base::include_pressure;
  using Base::periodic_tol; using Base::dim; using Base::hmin_;
  using Base::ncoeffs_u; using Base::ncoeffs_p;
  using Base::cells_; using Base::facet_neigh_; using Base::period_;
  using Base::u_dofs_; using Base::p_dofs_;
  using Base::x_min; using Base::x_max; using Base::is_initialized; using Base::t_update;
  using Base::get_folder; using Base::set_folder; using Base::wants_gradient;
  using Base::read_mesh_params; using Base::set_period; using Base::tree_;
  using Base::verbose;

  // They depend on the time alone; kept per thread for the times last asked
  const FreqWeights& weights(const double t);
  __attribute__((noinline)) const FreqWeights& fill_weights(const double t);

  std::uint64_t id_ = 0;   // tells this loader from one later at the same address
  void new_id();
  FreqStamps fs;   // frequencies holder
  double omega0 = 0.;

  // Per component, the field by node: D doubles a node for the velocity, one
  // for the pressure
  std::vector<std::vector<double>> u_nodes_;
  std::vector<std::vector<double>> p_nodes_;

  // The arrays the tree addresses; the cells hold the geometry
  std::vector<std::uint32_t> topo_;
  std::vector<double> coords_;
  std::size_t ncells_ = 0, nverts_ = 0;

  // Near walls: P1 to P2 by the edge rule, or P2 plus the quartic terms
  bool wall_upgrade_ = false;
  bool upgrade_p2_ = false;
  std::vector<std::uint32_t> vclass_;
  std::vector<std::int32_t> wall_index_;   // -1: no wall vertex at rest
  std::vector<near_wall::WallEdges<Cell>> wall_cells_;
  near_wall::UpgradeTable<Cell> upgrade_;
  std::vector<std::uint8_t> rest_;   // by wall cell: the vertices at rest in every component
  // The wall vertices at rest, from every component's values
  void rest_masks();
  // A wall cell at bary, out of line: the components summed with their weights,
  // then raised; component k's block at blocks + k*stride, gathered here if
  // blocks is null
  __attribute__((noinline))
  void wall_at(const int id, const std::array<double, 4>& bary, const FreqWeights& wf, const double* blocks,
               const std::size_t stride, PointValues& fields) const;
  // A P2 wall cell, out of line: P2 plus the terms of the summed blocks B (the
  // velocity's) and Bt (its rate's)
  __attribute__((noinline))
  void upgrade_at(const int id, const std::array<double, 4>& bary, const double* B, const double* Bt,
                  PointValues& fields) const;
  // Without gradients: the velocity and its rate at bary in cell id, out of line
  __attribute__((noinline))
  void velocity_at(const int id, const std::array<double, 4>& bary, const FreqWeights& wf,
                   PointValues& fields) const;
  // The pressure at bary in cell id
  double pressure_at(const int id, const std::array<double, 4>& bary, const FreqWeights& wf) const;
};

#endif
