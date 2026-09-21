#ifndef __MORTON_HPP
#define __MORTON_HPP

// Morton (Z-order) codes and the deterministic parallel sort that puts points
// or cells in that order. 21 bits an axis in 3D (a 63-bit code), 31 in 2D.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#if defined(__GNUC__) && defined(_OPENMP)
#include <parallel/algorithm>
#define PARTRAC_PARALLEL_SORT 1
#endif

#include "Error.hpp"
#include "typedefs.hpp"

namespace partrac {

constexpr int morton_bits_3d = 21;
constexpr int morton_bits_2d = 31;

// 21 bits of v to every third bit
inline std::uint64_t spread3(std::uint64_t v){
  v &= 0x1fffffULL;
  v = (v | (v << 32)) & 0x1f00000000ffffULL;
  v = (v | (v << 16)) & 0x1f0000ff0000ffULL;
  v = (v | (v <<  8)) & 0x100f00f00f00f00fULL;
  v = (v | (v <<  4)) & 0x10c30c30c30c30c3ULL;
  v = (v | (v <<  2)) & 0x1249249249249249ULL;
  return v;
}

// 31 bits of v to every other bit
inline std::uint64_t spread2(std::uint64_t v){
  v &= 0x7fffffffULL;
  v = (v | (v << 16)) & 0x0000ffff0000ffffULL;
  v = (v | (v <<  8)) & 0x00ff00ff00ff00ffULL;
  v = (v | (v <<  4)) & 0x0f0f0f0f0f0f0f0fULL;
  v = (v | (v <<  2)) & 0x3333333333333333ULL;
  v = (v | (v <<  1)) & 0x5555555555555555ULL;
  return v;
}

inline std::uint64_t morton3(const std::uint32_t x, const std::uint32_t y, const std::uint32_t z){
  return spread3(x) | (spread3(y) << 1) | (spread3(z) << 2);
}

inline std::uint64_t morton2(const std::uint32_t x, const std::uint32_t y){
  return spread2(x) | (spread2(y) << 1);
}

// The box a position is quantised in; a degenerate axis gives index 0
class MortonBox {
public:
  MortonBox() = default;
  MortonBox(const Vector3d& lo, const Vector3d& hi, const int dim) : dim_(dim) {
    const double n = double((std::uint64_t(1) << (dim == 3 ? morton_bits_3d : morton_bits_2d)) - 1);
    top_ = n;
    for (int d = 0; d < 3; ++d){
      lo_[d] = lo[d];
      const double span = hi[d] - lo[d];
      scale_[d] = span > 0. ? n/span : 0.;
    }
  }
  std::uint64_t code(const double* x) const {
    std::array<std::uint32_t, 3> q{{0, 0, 0}};
    for (int d = 0; d < dim_; ++d){
      const double v = (x[d] - lo_[d])*scale_[d];
      q[d] = std::uint32_t(v < 0. ? 0. : (v > top_ ? top_ : std::floor(v)));
    }
    return dim_ == 3 ? morton3(q[0], q[1], q[2]) : morton2(q[0], q[1]);
  }
  std::uint64_t code(const Vector3d& x) const { return code(x.data()); }
  int dim() const { return dim_; }
private:
  std::array<double, 3> lo_{{0., 0., 0.}}, scale_{{0., 0., 0.}};
  int dim_ = 3;
  double top_ = 0.;
};

// A code with where it came from; the order is total, so equal codes keep
// their index order whatever the sort does with them
struct MortonKey {
  std::uint64_t code = 0;
  std::uint32_t index = 0;
  bool operator<(const MortonKey& o) const {
    return code != o.code ? code < o.code : index < o.index;
  }
};

// Sorts in place, on both fields, in parallel where libstdc++ offers it
inline void morton_sort(std::vector<MortonKey>& keys){
#ifdef PARTRAC_PARALLEL_SORT
  __gnu_parallel::sort(keys.begin(), keys.end());
#else
  std::sort(keys.begin(), keys.end());
#endif
}

// Slot i of the returned permutation holds the index that sorts there
inline std::vector<std::uint32_t> morton_permutation(std::vector<MortonKey>& keys){
  morton_sort(keys);
  std::vector<std::uint32_t> perm(keys.size());
  for (std::size_t i = 0; i < keys.size(); ++i) perm[i] = keys[i].index;
  return perm;
}

// n positions, stride doubles apart, in Morton order of the box
inline std::vector<std::uint32_t> morton_order(const double* x, const std::size_t n,
                                               const std::size_t stride, const MortonBox& box){
  if (n > std::numeric_limits<std::uint32_t>::max())
    fail("morton_order: ", n, " points do not fit a 32-bit index");
  std::vector<MortonKey> keys(n);
#pragma omp parallel for schedule(static)
  for (std::ptrdiff_t i = 0; i < std::ptrdiff_t(n); ++i){
    keys[i].code = box.code(x + std::size_t(i)*stride);
    keys[i].index = std::uint32_t(i);
  }
  return morton_permutation(keys);
}

}  // namespace partrac

#endif
