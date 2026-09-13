#ifdef USE_DOLFIN
#ifndef __CELL_LOCATE_HPP
#define __CELL_LOCATE_HPP

#include <dolfin.h>
#include <array>
#include <cassert>
#include <limits>
#include <vector>
#include <omp.h>
#include <iostream>
#include "typedefs.hpp"

// The cell a point is in, tried in the order that is nearly always right: the
// cell it was in last, then that cell's neighbours, then the mesh's tree.
// Five interpolators carried this as five copies; the cell type is the only
// thing that differed. Cell needs contains(const Vector3d&). The counters are
// one per thread and are what print_found reports.
// The buffers in evaluate are sized from Cell::n_dofs_max, so an element
// richer than the basis routines know about would run off the end of them.
inline void check_dofs_fit(const Uint ncoeffs_u, const Uint ncoeffs_p,
                           const std::size_t n_dofs_max, const char* what){
  if (ncoeffs_u > n_dofs_max || ncoeffs_p > n_dofs_max){
    std::cout << what << ": element has " << ncoeffs_u << " and " << ncoeffs_p
              << " dofs, against a maximum of " << n_dofs_max << std::endl;
    exit(1);
  }
}

// A cell's face-neighbours: at most three for a triangle, four for a tet, a
// periodic image standing in for a missing one. Kept sorted and without
// repeats, as the std::set it replaces was, so a point on a shared edge is
// still found in the same cell.
struct CellNeighbours {
  std::array<Uint, 4> id{};
  unsigned char n = 0;
  void insert(const Uint c){
    unsigned char k = 0;
    while (k < n && id[k] < c) ++k;
    if (k < n && id[k] == c) return;
    assert(n < 4);
    for (unsigned char m = n; m > k; --m) id[m] = id[m-1];
    id[k] = c; ++n;
  }
  std::size_t size() const { return n; }
  const Uint* begin() const { return id.data(); }
  const Uint* end() const { return id.data() + n; }
};

template<typename Cell>
inline bool locate_in_cells(const std::vector<Cell>& cells,
                            const std::vector<CellNeighbours>& cell2cells,
                            const dolfin::Mesh& mesh,
                            const Uint dim,
                            const Vector3d& xx,
                            int& id_prev,
                            std::vector<long unsigned int>& found_same,
                            std::vector<long unsigned int>& found_nneigh,
                            std::vector<long unsigned int>& found_other){
  const int tid = omp_get_thread_num();
  if (id_prev >= 0){
    if (cells[id_prev].contains(xx)){
      ++found_same[tid];
      return true;
    }
    for ( auto neigh_id : cell2cells[id_prev] ){
      if (cells[neigh_id].contains(xx)){
        ++found_nneigh[tid];
        id_prev = neigh_id;
        return true;
      }
    }
  }
  const dolfin::Point point(dim, xx.data());
  const unsigned int id = mesh.bounding_box_tree()->compute_first_entity_collision(point);
  if (id == std::numeric_limits<unsigned int>::max())
    return false;
  ++found_other[tid];
  id_prev = id;
  return true;
}

#endif
#endif
