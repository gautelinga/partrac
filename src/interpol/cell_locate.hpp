#ifdef USE_DOLFIN
#ifndef __CELL_LOCATE_HPP
#define __CELL_LOCATE_HPP

#include <dolfin.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>
#include <iostream>
#include "typedefs.hpp"
#include "PointValues.hpp"

// Buffers in evaluate are sized from Cell::n_dofs_max
inline void check_dofs_fit(const Uint ncoeffs_u, const Uint ncoeffs_p,
                           const std::size_t n_dofs_max, const char* what){
  if (ncoeffs_u > n_dofs_max || ncoeffs_p > n_dofs_max){
    std::cout << what << ": element has " << ncoeffs_u << " and " << ncoeffs_p
              << " dofs, against a maximum of " << n_dofs_max << std::endl;
    exit(1);
  }
}

// Dof indices of all cells, flat, fixed stride
class CellDofs {
public:
  void build(const dolfin::GenericDofMap& dofmap,
             const std::vector<dolfin::Cell>& cells, const char* what){
    stride_ = cells.empty() ? 0 : dofmap.cell_dofs(cells[0].index()).size();
    dofs_.resize(cells.size()*stride_);
    for (std::size_t id = 0; id < cells.size(); ++id){
      const auto dofs = dofmap.cell_dofs(cells[id].index());
      if (std::size_t(dofs.size()) != stride_){
        std::cout << what << ": cell " << id << " has " << dofs.size()
                  << " dofs, the first cell " << stride_ << std::endl;
        exit(1);
      }
      for (std::size_t i = 0; i < stride_; ++i){
        if (dofs[i] < 0 || std::uint64_t(dofs[i]) > std::numeric_limits<std::uint32_t>::max()){
          std::cout << what << ": dof index " << dofs[i] << " does not fit 32 bits" << std::endl;
          exit(1);
        }
        dofs_[id*stride_ + i] = std::uint32_t(dofs[i]);
      }
    }
  }
  // At least n dofs per cell
  void check_stride(const std::size_t n, const char* what) const {
    if (stride_ < n){
      std::cout << what << ": " << stride_ << " dofs per cell, evaluate reads " << n << std::endl;
      exit(1);
    }
  }
  const std::uint32_t* operator[](const std::size_t id) const { return dofs_.data() + id*stride_; }
  std::size_t stride() const { return stride_; }
private:
  std::vector<std::uint32_t> dofs_;
  std::size_t stride_ = 0;
};

// Every cell's dofs, sorted, flat, fixed stride
inline std::vector<int> sorted_dof_table(const dolfin::GenericDofMap& dofmap,
                                         const std::size_t ncells, std::size_t& stride){
  stride = ncells ? dofmap.cell_dofs(0).size() : 0;
  std::vector<int> table(ncells * stride);
  for (std::size_t i = 0; i < ncells; ++i){
    const auto d = dofmap.cell_dofs(i);
    if (std::size_t(d.size()) != stride){
      std::cout << "cell " << i << " has " << d.size() << " dofs, the first cell " << stride << std::endl;
      exit(1);
    }
    int* row = table.data() + i*stride;
    std::copy(d.data(), d.data() + stride, row);
    std::sort(row, row + stride);
  }
  return table;
}

// Fraction of consecutive cells sharing a dof; low in a poorly ordered mesh
inline double dof_sharing(const dolfin::GenericDofMap& dofmap, const std::size_t ncells){
  if (ncells < 2) return 1.;
  std::size_t stride = 0;
  const std::vector<int> table = sorted_dof_table(dofmap, ncells, stride);
  std::size_t shared = 0;
  for (std::size_t i = 1; i < ncells; ++i){
    const int* prev = table.data() + (i-1)*stride;
    const int* cur = table.data() + i*stride;
    std::size_t a = 0, b = 0;
    while (a < stride && b < stride){
      if (prev[a] == cur[b]){ ++shared; break; }
      (prev[a] < cur[b]) ? ++a : ++b;
    }
  }
  return double(shared) / (ncells - 1);
}

// Cells in the order of their sorted dofs; returns map[old] -> new
inline std::vector<std::uint32_t> order_cells_by_dofs(const dolfin::GenericDofMap& dofmap, const std::size_t ncells){
  std::size_t stride = 0;
  const std::vector<int> table = sorted_dof_table(dofmap, ncells, stride);
  std::vector<std::uint32_t> by_key(ncells);
  std::iota(by_key.begin(), by_key.end(), 0);
  std::stable_sort(by_key.begin(), by_key.end(),
                   [&](const std::uint32_t a, const std::uint32_t b){
                     const int* ra = table.data() + a*stride;
                     const int* rb = table.data() + b*stride;
                     return std::lexicographical_compare(ra, ra + stride, rb, rb + stride);
                   });
  std::vector<std::uint32_t> map(ncells);
  for (std::size_t l = 0; l < ncells; ++l) map[by_key[l]] = l;
  return map;
}

// Cell order: dolfin's, unless consecutive cells rarely share a dof
inline std::vector<std::uint32_t> cell_order(const dolfin::GenericDofMap& dofmap,
                                            const std::size_t ncells,
                                            const std::string& mode_in,
                                            std::vector<std::uint32_t>& dolfin2local){
  std::vector<std::uint32_t> order(ncells);
  std::iota(order.begin(), order.end(), 0);
  const std::string mode = mode_in.empty() ? "auto" : mode_in;
  if (mode != "auto" && mode != "never" && mode != "always"){
    std::cout << "renumber_cells must be auto, never or always, not " << mode << std::endl;
    exit(1);
  }
  bool renumber = mode == "always";
  if (mode == "auto"){
    const double sharing = dof_sharing(dofmap, ncells);
    renumber = sharing < 0.5;
    std::cout << "Consecutive cells sharing a dof: " << sharing << std::endl;
  }
  if (renumber){
    std::cout << "Cell order: renumbering cells by dofs" << std::endl;
    dolfin2local = order_cells_by_dofs(dofmap, ncells);
    for (std::size_t i = 0; i < ncells; ++i) order[dolfin2local[i]] = i;
  }
  return order;
}

// How locate found a point, for the tests
struct FoundCounts {
  long unsigned int same = 0, walk = 0, tree = 0;
};

// Facet k of a cell faces vertex k; across it: a cell, a wall, or a periodic image
constexpr std::int32_t facet_wall = -1;
inline std::int32_t facet_periodic(const std::int32_t id){ return -2 - id; }

// From the known cell, try its periodic partners, then step across the facet
// of the most negative barycentric; false at a wall or after max_walk steps
template<typename Cell>
inline bool walk_to_cell(const std::vector<Cell>& cells,
                         const std::vector<std::int32_t>& across,
                         const Vector3d& xx,
                         CellPos& pos,
                         FoundCounts* count = nullptr){
  constexpr int max_walk = 8;
  constexpr int nv = Cell::n_verts;
  if (pos.id < 0)
    return false;
  // On failure bary stays the stale cell's
  if (cells[pos.id].contains(xx, pos.bary)){
    if (count) ++count->same;
    return true;
  }
  std::array<double, 4> bary = pos.bary;
  int id = pos.id;
  for (int step = 0; step < max_walk; ++step){
    const std::int32_t* row = across.data() + std::size_t(id)*nv;
    // Periodic partners: the wrapped point lies across the box
    for (int k = 0; k < nv; ++k){
      if (row[k] >= facet_wall) continue;
      const int j = facet_periodic(row[k]);
      std::array<double, 4> b;
      if (cells[j].contains(xx, b)){
        if (count) ++count->walk;
        pos.id = j;
        pos.bary = b;
        return true;
      }
    }
    int k_exit = 0;
    for (int k = 1; k < nv; ++k)
      if (bary[k] < bary[k_exit]) k_exit = k;
    const std::int32_t a = row[k_exit];
    if (a < 0)
      return false;
    id = a;
    if (cells[id].contains(xx, bary)){
      if (count) ++count->walk;
      pos.id = id;
      pos.bary = bary;
      return true;
    }
  }
  return false;
}

// The bounding-box tree: every cell, from nothing known
template<typename Cell>
inline bool tree_to_cell(const std::vector<Cell>& cells,
                         const dolfin::Mesh& mesh,
                         const Uint dim,
                         const Vector3d& xx,
                         CellPos& pos,
                         const std::vector<std::uint32_t>* dolfin2local = nullptr,
                         FoundCounts* count = nullptr){
  const dolfin::Point point(dim, xx.data());
  const unsigned int id = mesh.bounding_box_tree()->compute_first_entity_collision(point);
  if (id == std::numeric_limits<unsigned int>::max())
    return false;
  if (count) ++count->tree;
  pos.id = dolfin2local ? int((*dolfin2local)[id]) : int(id);
  // Tree tolerance: may sit just outside
  cells[pos.id].contains(xx, pos.bary);
  return true;
}

// Walk dx from x, mirroring at walls
template<typename Cell, typename Wrap>
inline bool reflect_in_cells(const std::vector<Cell>& cells,
                             const std::vector<std::int32_t>& across,
                             const int nv,
                             const Vector3d& period,
                             const Vector3d& x,
                             Vector3d& dx,
                             CellPos& pos,
                             const Wrap& wrap){
  constexpr int max_bounces = 8;
  constexpr int max_crossings = 1024;
  if (pos.id < 0)
    return false;
  int id = pos.id;
  Vector3d p = wrap(x);
  Vector3d d = dx;
  Vector3d walked = Vector3d::Zero();
  std::array<double, 4> b0 = pos.bary, b1;
  int bounces = 0;
  for (int crossing = 0; crossing < max_crossings; ++crossing){
    if (cells[id].contains(p + d, b1)){
      pos.id = id;
      pos.bary = b1;
      dx = walked + d;
      return true;
    }
    // First facet the segment leaves
    int k_exit = 0;
    double s = 2.;
    for (int k = 0; k < nv; ++k){
      if (b1[k] < 0.){
        const double sk = (b0[k] > 0. && b1[k] < b0[k]) ? b0[k]/(b0[k] - b1[k]) : 0.;
        if (sk < s){ s = sk; k_exit = k; }
      }
    }
    const Vector3d part = s*d;
    Vector3d rest = d - part;
    p += part;
    walked += part;
    const std::int32_t a = across[std::size_t(id)*nv + k_exit];
    if (a == facet_wall){
      if (++bounces > max_bounces)
        return false;
      const Vector3d n = cells[id].bary_grad(k_exit).normalized();
      const double rn = n.dot(rest);
      if (rn < 0.)
        rest -= 2*rn*n;
    }
    else if (a >= 0){
      id = a;
    }
    else {
      // Periodic image of p
      const Vector3d n = cells[id].bary_grad(k_exit);
      int axis = 0;
      n.cwiseAbs().maxCoeff(&axis);
      p[axis] += (n[axis] > 0. ? 1. : -1.)*period[axis];
      id = facet_periodic(a);
    }
    d = rest;
    cells[id].contains(p, b0);
  }
  return false;
}

#endif
#endif
