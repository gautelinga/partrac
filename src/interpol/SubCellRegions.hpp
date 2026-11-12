#ifndef __SUBCELLREGIONS_HPP
#define __SUBCELLREGIONS_HPP

// Regions (regions.hpp) of a simplex mesh loader whose velocity is a
// polynomial per sub-cell of the barycentric split in some cells: every cell
// of the divergence-free split, upgraded P2's wall cells with a table. Such a
// cell's regions are its sub-cells (sub: the macro vertex left out), any other
// cell one region (sub -1), as MeshCore's.
// Levels in the cell's barycentric units: a sub-cell's barycentrics mu
// (split_eval.hpp), the macro facet's divided by n_verts, its least barycentric.
// Split cells: Derived::split_cell(id), marked by mark_split, or every one with Every.

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <string>
#include <vector>
#include "MeshCore.hpp"
#include "regions.hpp"
#include "split_eval.hpp"

template<typename Cell, typename Derived, bool Every = false>
class SubCellRegions : public MeshCore<Cell> {
public:
  SubCellRegions(const std::string& infilename) : MeshCore<Cell>(infilename) {}
  using region_type = SubRegion;
  using typename MeshCore<Cell>::levels_type;
  bool region_of(const int id, const Vector3d& x, const double band, SubRegion& R, levels_type& lev){
    if (!Base::region_of(id, x, band, R, lev))
      return false;
    R.sub = sub_at(R.id, lev);
    return true;
  }
  // The sub-cell of the least barycentric, and the levels at x
  SubRegion region_in(const int id, const Vector3d& x, levels_type& lev) const {
    SubRegion R;
    static_cast<Region&>(R) = Base::region_in(id, x);
    Base::region_point(R, x, lev);
    R.sub = sub_at(R.id, lev);
    return R;
  }
  SubRegion region_in(const int id, const Vector3d& x) const {
    levels_type lev;
    return region_in(id, x, lev);
  }
  Vector3d region_point(const SubRegion& R, const Vector3d& x, levels_type& lev) const {
    if (whole(R))
      return Base::region_point(R, x, lev);
    assert(R.sub >= 0);
    levels_type b;
    const Vector3d xe = Base::region_point(R, x, b);
    region_levels(R.sub, b, lev);
    return xe;
  }
  // A sub-cell's internal plane to the sub-cell of the same cell; a facet,
  // a sub-cell's macro one, to the neighbour, its sub-cell on it if split
  int across(const SubRegion& R, const int k, const Vector3d& x, SubRegion& next) const {
    if (!whole(R) && k < nv - 1){
      assert(R.sub >= 0);
      next = R;
      next.sub = split_eval::sub_beyond(R.sub, k);
      return split_eval::sub_plane(next.sub, R.sub);
    }
    const int e = Base::across(R, whole(R) ? k : R.sub, x, next);
    if (e < 0)
      return -1;
    next.sub = split(next.id) ? e : -1;
    return next.sub < 0 ? e : nv - 1;
  }
  void level_rates(const SubRegion& R, const Vector3d& v, levels_type& rate) const {
    if (whole(R)){
      Base::level_rates(R, v, rate);
      return;
    }
    assert(R.sub >= 0);
    levels_type d;
    Base::level_rates(R, v, d);
    for (int m = 0; m < nv - 1; ++m) rate[m] = d[split_eval::sub_beyond(R.sub, m)] - d[R.sub];
    rate[nv - 1] = d[R.sub];
  }
  Vector3d level_grad(const SubRegion& R, const int k) const {
    if (whole(R))
      return Base::level_grad(R, k);
    assert(R.sub >= 0);
    const Vector3d gi = Base::level_grad(R, R.sub);
    return k < nv - 1 ? Vector3d(Base::level_grad(R, split_eval::sub_beyond(R.sub, k)) - gi) : gi;
  }
  bool is_wall(const SubRegion& R, const int k) const {
    if (whole(R))
      return Base::is_wall(R, k);
    assert(R.sub >= 0);
    return k == nv - 1 && Base::is_wall(R, R.sub);
  }
  // A split cell's an n_verts-th: a path meets each sub-cell once at most;
  // at least the floor (MeshCore)
  double region_size(const int id) const {
    const double s = split(id) ? Base::cell_size(id)/nv : Base::cell_size(id);
    return std::max(s, double(this->size_floor_[std::size_t(id)]));
  }
  // The whole cells as regions, MeshCore's: all of a cell's walls
  const MeshCore<Cell>& whole_cells() const { return *this; }
  // The held nodes the whole cell's, kept for another sub-cell (held_eval.hpp)
  static constexpr bool held_by_cell = !Every;
  template<typename Held>
  void hold_sub(const SubRegion& R, Held& h) const {
    static_assert(!Every, "a sub-cell's own nodes are held");
    h.sub = R.sub;
  }
  // Sub-cell i's levels from the cell's barycentrics: the last bary[i]
  static void region_levels(const int i, const levels_type& bary, levels_type& lev){
    split_eval::sub_levels<nv, 1>(i, bary, lev.data());
  }
  // The cell's barycentrics at the levels of region sub
  static levels_type cell_levels(const int sub, const levels_type& lev){
    levels_type bary = lev;
    if (sub >= 0){
      bary[sub] = lev[nv - 1];
      for (int m = 0; m < nv - 1; ++m) bary[split_eval::sub_beyond(sub, m)] = lev[m] + lev[nv - 1];
    }
    return bary;
  }
  // A sub-cell's barycentrics mu at its levels
  static void sub_bary(const levels_type& lev, double* mu){
    for (int m = 0; m < nv - 1; ++m) mu[m] = lev[m];
    mu[nv - 1] = double(nv)*lev[nv - 1];
  }
protected:
  using Base = MeshCore<Cell>;
  static constexpr int nv = Cell::n_verts;
  // Derived::split_cell of every cell, once its tables are built
  void mark_split(const std::size_t ncells){
    if constexpr (!Every){
      split_cells_.resize(ncells);
      const Derived& d = static_cast<const Derived&>(*this);
      for (std::size_t i = 0; i < ncells; ++i) split_cells_[i] = d.split_cell(int(i));
    }
  }
  bool split(const int id) const {
    if constexpr (Every) return true;
    else return split_cells_[std::size_t(id)];
  }
  bool whole(const SubRegion& R) const { return !Every && R.sub < 0; }
  // Cell id's sub-cell of the least barycentric, lev the cell's barycentrics
  // and then its levels; -1, lev kept, if not split
  int sub_at(const int id, levels_type& lev) const {
    if (!split(id))
      return -1;
    const levels_type b = lev;
    const int i = split_eval::least<nv>(b);
    region_levels(i, b, lev);
    return i;
  }
  std::vector<std::uint8_t> split_cells_;   // by cell, unless Every
};

#endif
