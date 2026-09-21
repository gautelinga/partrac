#ifndef __SIMPLEXFREQINTERPOL_HPP
#define __SIMPLEXFREQINTERPOL_HPP

// A time series given as frequency components: one steady Taylor-Hood pair per
// component, each a dolfin HDF5 checkpoint on the same triangle or tet mesh,
// summed with a cosine of the base frequency. The components are read the way
// SimplexInterpol reads a stamp and held by node, so an evaluation gathers the
// cell's nodes out of each component in turn.

#include <cstdint>
#include <string>
#include <vector>

#include "MeshCore.hpp"
#include "FreqStamps.hpp"
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
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  // What a step reads: velocity, acceleration and their gradients; P, Phi, cell_type stay zero
  void evaluate_motion(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  double get_t_min() { return dolfin_params.template get<double>("t_min"); };
  double get_t_max() { return dolfin_params.template get<double>("t_max"); };
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

  // They depend on the time alone; kept per thread for the times last asked
  const FreqWeights& weights(const double t);
  __attribute__((noinline)) const FreqWeights& fill_weights(const double t);

  std::uint64_t id_ = 0;   // tells this loader from one later at the same address
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
};

#endif
