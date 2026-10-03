// Cases the region and RK4cells tests share; writers in region_cases.cpp
#pragma once
#ifdef USE_DOLFIN

#include <array>
#include <cstring>
#include <functional>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "PointValues.hpp"
#include "case_dir.hpp"
#include "typedefs.hpp"

namespace dolfin { class Mesh; }

// A map of the unit mesh's points; a velocity at stamp k
using Warp = std::function<Vector3d(const Vector3d&)>;
using Field = std::function<Vector3d(const Vector3d&, int)>;

// dolfin's unit square (n x n) or cube (n x n x n), its connectivity built
template<typename Cell>
std::shared_ptr<dolfin::Mesh> unit_mesh(std::size_t n);

// Two stamps at t = 0 and 1 as a checkpoint (h5_params.dat) and, for a P1
// velocity, as XDMF (xdmf_params.dat); the mesh's points warped, the velocity given
template<typename Cell>
void write_stamped(const CaseDir& c, std::size_t n, const std::string& u_el, const std::vector<bool>& per,
                   bool rest, const Warp& warp = {}, const Field& field = {});

// Two frequency components, periodic along every axis
template<typename Cell>
void write_freq(const CaseDir& c, std::size_t n);

// Two stamps at t = 0 and 1 for divfree = true (h5_params.dat), periodic along x or walled
template<typename Cell>
void write_split(const CaseDir& c, std::size_t n, bool periodic);

// A felbm lattice of n^3 nodes, one unit apart and periodic, with a solid
// layer at z index 0, a block, and a lone solid node; two stamps at t = 0 and
// 1 of a smooth velocity, zero on the solid nodes (felbm_params.dat)
void write_lattice(const CaseDir& c, int n);

namespace {

template<typename Cell> constexpr int dim_of = Cell::n_verts - 1;

template<class M>
bool same_bits(const M& a, const M& b){ return std::memcmp(a.data(), b.data(), sizeof(double)*a.size()) == 0; }

template<typename Cell>
std::array<double, 4> random_weights(std::mt19937& rng, const double lo = 0.1){
  std::uniform_real_distribution<double> u(lo, 1.);
  std::array<double, 4> w{};
  double sum = 0.;
  for (int k = 0; k < Cell::n_verts; ++k) sum += (w[k] = u(rng));
  for (int k = 0; k < Cell::n_verts; ++k) w[k] /= sum;
  return w;
}

inline bool lattice_solid(const int n, const int i, const int j, const int k){
  return k == 0 || (i >= n/3 && i <= n/2 && j >= n/4 && j <= 2*n/3 && k >= n/6 && k <= 3*n/4)
      || (i == 3*n/4 && j == 3*n/4 && k == n/2);
}

// Points of the fluid over the lattice, every other one a few periods off
template<typename I>
std::vector<Vector3d> lattice_points(I& intp, const int n, const std::size_t count){
  std::mt19937 rng(17);
  std::uniform_real_distribution<double> uni(0., double(n));
  std::uniform_int_distribution<int> image(-2, 2);
  std::vector<Vector3d> pts;
  while (pts.size() < count){
    Vector3d x(uni(rng), uni(rng), uni(rng));
    if (pts.size() % 2) x += double(n)*Vector3d(image(rng), image(rng), image(rng));
    CellPos pos;
    if (intp.locate(x, 0., pos)) pts.push_back(x);
  }
  return pts;
}

template<typename I>
std::vector<Vector3d> box_points(I& intp, const int D, const std::vector<bool>& per, const std::size_t n){
  std::mt19937 rng(11);
  std::uniform_real_distribution<double> uni(0.001, 0.999);
  std::uniform_int_distribution<int> image(-2, 2);
  const Vector3d lo = intp.get_x_min(), hi = intp.get_x_max();
  std::vector<Vector3d> pts;
  for (std::size_t i = 0; i < n; ++i){
    Vector3d x = Vector3d::Zero();
    for (int d = 0; d < D; ++d){
      x[d] = lo[d] + uni(rng)*(hi[d] - lo[d]);
      if (per[std::size_t(d)] && i % 2) x[d] += image(rng)*(hi[d] - lo[d]);
    }
    pts.push_back(x);
  }
  return pts;
}
}  // namespace

#endif
