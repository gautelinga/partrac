#ifndef __CELL_WALK_HPP
#define __CELL_WALK_HPP

// The facet table's encoding, the per-cell dof table and the two walks over
// them: locating a point from a known cell, and mirroring a move at the walls.
// Nothing here needs the mesh's source, so it compiles without dolfin.

#include <array>
#include <cstdint>
#include <vector>
#include "Error.hpp"
#include "typedefs.hpp"
#include "PointValues.hpp"

#ifdef USE_DOLFIN
namespace dolfin { class Cell; class GenericDofMap; }
#endif

// Buffers in evaluate are sized from Cell::n_dofs_max
inline void check_dofs_fit(const Uint ncoeffs_u, const Uint ncoeffs_p,
                           const std::size_t n_dofs_max, const char* what){
  if (ncoeffs_u > n_dofs_max || ncoeffs_p > n_dofs_max){
    partrac::fail(what, ": element has ", ncoeffs_u, " and ", ncoeffs_p, " dofs, against a maximum of ", n_dofs_max);
  }
}

// Dof indices of all cells, flat, fixed stride
class CellDofs {
public:
#ifdef USE_DOLFIN
  // Defined in dolfin_ref.cpp, the one place that needs the dofmap
  void build(const dolfin::GenericDofMap& dofmap,
             const std::vector<dolfin::Cell>& cells, const char* what);
#endif
  // The node ids of every cell: the nv vertices of the topology, then, for a
  // quadratic field, the cell's edges in dolfin's local order. node_map, when
  // given, renumbers the nodes (vertices first, then edges offset by nverts).
  void fill(const std::vector<std::uint32_t>& topo, const std::vector<std::uint32_t>& edges,
            const std::size_t ncells, const int nv, const int ne, const bool quadratic,
            const std::size_t nverts, const std::uint32_t* node_map){
    stride_ = std::size_t(nv) + (quadratic ? std::size_t(ne) : 0);
    dofs_.resize(ncells*stride_);
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t i = 0; i < std::ptrdiff_t(ncells); ++i){
      std::uint32_t* row = dofs_.data() + std::size_t(i)*stride_;
      for (int k = 0; k < nv; ++k){
        const std::uint32_t v = topo[std::size_t(i)*std::size_t(nv) + std::size_t(k)];
        row[k] = node_map ? node_map[v] : v;
      }
      if (!quadratic) continue;
      for (int e = 0; e < ne; ++e){
        const std::uint32_t n = std::uint32_t(nverts) + edges[std::size_t(i)*std::size_t(ne) + std::size_t(e)];
        row[nv + e] = node_map ? node_map[n] : n;
      }
    }
  }
  // At least n dofs per cell
  void check_stride(const std::size_t n, const char* what) const {
    if (stride_ < n){
      partrac::fail(what, ": ", stride_, " dofs per cell, evaluate reads ", n);
    }
  }
  const std::uint32_t* operator[](const std::size_t id) const { return dofs_.data() + id*stride_; }
  std::size_t stride() const { return stride_; }
  // The flat table, for a loader that caches it
  const std::vector<std::uint32_t>& table() const { return dofs_; }
  void adopt(std::vector<std::uint32_t>&& table, const std::size_t stride){
    dofs_ = std::move(table);
    stride_ = stride;
  }
private:
  std::vector<std::uint32_t> dofs_;
  std::size_t stride_ = 0;
};

// How locate found a point, for the tests
struct FoundCounts {
  long unsigned int same = 0, walk = 0, tree = 0;
};

// Facet k of a cell faces vertex k; across it: a cell, a wall, or a periodic image
constexpr std::int32_t facet_wall = -1;
constexpr std::int32_t facet_periodic(const std::int32_t id){ return -2 - id; }

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
