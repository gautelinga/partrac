#ifdef USE_DOLFIN
#ifndef __CELL_LOCATE_HPP
#define __CELL_LOCATE_HPP

#include <dolfin.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <limits>
#include <vector>
#include <omp.h>
#include <iostream>
#include "typedefs.hpp"

// Buffers in evaluate are sized from Cell::n_dofs_max
inline void check_dofs_fit(const Uint ncoeffs_u, const Uint ncoeffs_p,
                           const std::size_t n_dofs_max, const char* what){
  if (ncoeffs_u > n_dofs_max || ncoeffs_p > n_dofs_max){
    std::cout << what << ": element has " << ncoeffs_u << " and " << ncoeffs_p
              << " dofs, against a maximum of " << n_dofs_max << std::endl;
    exit(1);
  }
}

// Face neighbours, sorted and unique (at most 4)
struct CellNeighbours {
  std::array<std::uint32_t, 4> id{};
  unsigned char n = 0;
  void insert(const Uint c){
    unsigned char k = 0;
    while (k < n && id[k] < c) ++k;
    if (k < n && id[k] == c) return;
    assert(n < 4);
    for (unsigned char m = n; m > k; --m) id[m] = id[m-1];
    id[k] = std::uint32_t(c); ++n;
  }
  std::size_t size() const { return n; }
  const std::uint32_t* begin() const { return id.data(); }
  const std::uint32_t* end() const { return id.data() + n; }
};

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

// Per-thread locate counters, cache-line aligned
struct alignas(64) FoundCounts {
  long unsigned int same = 0, nneigh = 0, other = 0;
};

inline void print_found_counts(std::vector<FoundCounts>& found){
  long unsigned int same = 0, nneigh = 0, other = 0;
  for (const auto & f : found){
    same += f.same; nneigh += f.nneigh; other += f.other;
  }
  long int found_sum = same + nneigh + other;
  double frac_same = double(same) / found_sum;
  double frac_nneigh = double(nneigh) / found_sum;
  double frac_other = 1. - frac_same - frac_nneigh;
  std::cout << "Found in same cell: " << frac_same << ", nearest neighbour cell: " << frac_nneigh << ", other cell: " << frac_other << std::endl;
  std::fill(found.begin(), found.end(), FoundCounts{});
}

template<typename Cell>
inline bool locate_in_cells(const std::vector<Cell>& cells,
                            const std::vector<CellNeighbours>& cell2cells,
                            const dolfin::Mesh& mesh,
                            const Uint dim,
                            const Vector3d& xx,
                            CellPos& pos,
                            std::vector<FoundCounts>& found){
  FoundCounts& count = found[omp_get_thread_num()];
  if (pos.id >= 0){
    // On failure bary stays the stale cell's
    if (cells[pos.id].contains(xx, pos.bary)){
      ++count.same;
      return true;
    }
    std::array<double, 4> bary;
    for ( auto neigh_id : cell2cells[pos.id] ){
      if (cells[neigh_id].contains(xx, bary)){
        ++count.nneigh;
        pos.id = neigh_id;
        pos.bary = bary;
        return true;
      }
    }
  }
  const dolfin::Point point(dim, xx.data());
  const unsigned int id = mesh.bounding_box_tree()->compute_first_entity_collision(point);
  if (id == std::numeric_limits<unsigned int>::max())
    return false;
  ++count.other;
  pos.id = id;
  // Tree tolerance: may sit just outside
  cells[id].contains(xx, pos.bary);
  return true;
}

#endif
#endif
