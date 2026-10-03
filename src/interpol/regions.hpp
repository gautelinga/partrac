#ifndef __REGIONS_HPP
#define __REGIONS_HPP

// Regions of a simplex mesh: a cell and the offset that carries a point,
// never wrapped, into that cell's periodic image. Levels are the cell's
// barycentrics at the evaluation point; their rates along a velocity come
// from the constant barycentric gradients. Free functions over the cell and
// facet tables, as the walks in cell_walk.hpp.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>
#include "cell_walk.hpp"
#include "typedefs.hpp"

struct Region {
  int id = -1;
  Vector3d offset = Vector3d::Zero();
};

// The evaluation point of x in R, its levels in lev: the wrapped point where
// the cell holds it, as locate finds it, else x + offset
template<typename Cell, typename Wrap>
inline Vector3d region_point(const std::vector<Cell>& cells, const Region& R, const Vector3d& x,
                             const Wrap& wrap, std::array<double, 4>& lev){
  const Vector3d xx = wrap(x);
  if (cells[R.id].contains(xx, lev))
    return xx;
  const Vector3d xo = x + R.offset;
  cells[R.id].contains(xo, lev);
  return xo;
}

template<int nv, std::size_t N>
inline double min_level(const std::array<double, N>& lev){
  double m = lev[0];
  for (int k = 1; k < nv; ++k) m = std::min(m, lev[k]);
  return m;
}

// The region of x from the cell last known (id, -1 if none): that cell if its
// levels are down to -band, else that cell across a periodic face, else the
// cell relocate finds for the wrapped point; false if none holds x
template<typename Cell, typename Wrap, typename Relocate>
inline bool region_of(const std::vector<Cell>& cells, const Vector3d& period, const int id,
                      const Vector3d& x, const double band, const Wrap& wrap,
                      const Relocate& relocate, Region& R, std::array<double, 4>& lev){
  constexpr int nv = Cell::n_verts;
  const Vector3d xx = wrap(x);
  R.offset = xx - x;
  if (id >= 0){
    R.id = id;
    region_point(cells, R, x, wrap, lev);
    if (min_level<nv>(lev) >= -band)
      return true;
    // On a periodic face: the stored cell on the other side
    for (int d = 0; d < 3; ++d){
      if (period[d] <= 0.) continue;
      for (const double s : {period[d], -period[d]}){
        Region S = R;
        S.offset[d] += s;
        std::array<double, 4> l;
        cells[id].contains(x + S.offset, l);
        if (min_level<nv>(l) >= -band){
          R = S;
          lev = l;
          return true;
        }
      }
    }
  }
  CellPos pos;
  pos.id = id;
  if (!relocate(xx, pos))
    return false;
  R.id = pos.id;
  lev = pos.bary;
  return true;
}

// The region across facet k of R, into next; returns the facet entered, -1 at
// a wall. The entered facet is the one of the neighbour's row naming the cell
// left; of two (a box one cell thick), the level nearest zero at x
template<typename Cell, typename Wrap>
inline int region_across(const std::vector<Cell>& cells, const std::vector<std::int32_t>& across,
                         const Vector3d& period, const Region& R, const int k, const Vector3d& x,
                         const Wrap& wrap, Region& next){
  constexpr int nv = Cell::n_verts;
  const std::int32_t a = across[std::size_t(R.id)*nv + std::size_t(k)];
  if (a == facet_wall)
    return -1;
  std::int32_t back;
  next.offset = R.offset;
  if (a >= 0){
    next.id = a;
    back = R.id;
  }
  else {
    int axis;
    const double shift = periodic_shift(cells[R.id].bary_grad(k), period, axis);
    next.offset[axis] += shift;
    next.id = facet_periodic(a);
    back = facet_periodic(R.id);
  }
  const std::int32_t* row = across.data() + std::size_t(next.id)*nv;
  int entered = -1;
  std::array<double, 4> lev;
  for (int m = 0; m < nv; ++m){
    if (row[m] != back) continue;
    if (entered >= 0){
      // Twice: at the landing point
      region_point(cells, next, x, wrap, lev);
      if (std::abs(lev[m]) < std::abs(lev[entered])) entered = m;
    }
    else entered = m;
  }
  return entered;
}

// The levels' rates along v
template<typename Cell>
inline void level_rates(const Cell& cell, const Vector3d& v, std::array<double, 4>& rate){
  for (int k = 0; k < Cell::n_verts; ++k)
    rate[k] = cell.bary_grad(k).dot(v);
}

#endif
