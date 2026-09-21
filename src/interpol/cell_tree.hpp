#ifndef __CELL_TREE_HPP
#define __CELL_TREE_HPP

// A cell tree (Garth & Joy 2010) over tets or triangles addressed by their
// index in the caller's arrays, with Shewchuk's exact predicates at the leaves.
// The tree holds the topology and the coordinates by pointer; the caller keeps
// them alive. Queries are const and run from any number of threads.

#include <cstdint>
#include <utility>
#include <vector>

#include "Error.hpp"
#include "morton.hpp"
#include "predicates/predicates.hpp"
#include "typedefs.hpp"

namespace partrac {

// exactinit, once; every caller of cell_contains_exact needs it run first,
// and the tree's constructor does it before any thread starts
void init_exact_predicates();

// Exact containment. Putting x in vertex slot k gives bary_k times the cell's
// own determinant, so x is in the cell when no two of the nv determinants
// disagree in sign; a zero is a point on the boundary and still contained.
// All nv zero is a degenerate cell, which contains nothing.
inline bool cell_contains_exact(const std::uint32_t* row, const int nv, const double* coords,
                                const std::size_t stride, const Vector3d& x){
  double* v[4];
  for (int k = 0; k < nv; ++k) v[k] = const_cast<double*>(coords + std::size_t(row[k])*stride);
  double* const p = const_cast<double*>(x.data());
  bool pos = false, neg = false;
  for (int k = 0; k < nv; ++k){
    double* save = v[k];
    v[k] = p;
    const double det = (nv == 4) ? orient3d(v[0], v[1], v[2], v[3]) : orient2d(v[0], v[1], v[2]);
    v[k] = save;
    if (det > 0.){ if (neg) return false; pos = true; }
    else if (det < 0.){ if (pos) return false; neg = true; }
  }
  return pos || neg;
}

// A node covers a contiguous range of the tree's cell list. An interior node
// splits it along axis into a left and a right part, overlapping between the
// two planes; a leaf holds the range itself.
struct CellTreeNode {
  static constexpr std::uint8_t leaf_axis = 3;
  float lmax = 0.f;        // the left part's largest coordinate along axis
  float rmin = 0.f;        // the right part's smallest, both padded outward
  std::uint32_t a = 0;     // interior: the left child, the right is the next
  std::uint32_t n = 0;     // leaf: the number of cells from a
  std::uint8_t axis = leaf_axis;
  bool is_leaf() const { return axis == leaf_axis; }
};

class CellTree {
public:
  static constexpr std::uint32_t leaf_cells = 8;      // cells a leaf holds at most
  static constexpr int n_buckets = 16;                // split candidates a node tries
  static constexpr std::uint32_t task_cells = 4096;   // subtrees below this are tasks
  static constexpr int max_depth = 48;                // a deeper range becomes a leaf

  // topology: ncells rows of nv vertex indices; coords: npoints rows of
  // coord_stride doubles, which is the cell's own dimension, nv - 1; neither is
  // copied. verbose prints the tree's quality on one line.
  CellTree(const std::uint32_t* topology, std::size_t ncells, int nv,
           const double* coords, std::size_t npoints,
           std::size_t coord_stride = 3, bool verbose = false);

  // The containing cell of lowest id, or -1; independent of the traversal
  int locate(const Vector3d& x, std::size_t* tested = nullptr) const;

  // The exact test on one cell, for a caller that already knows it
  bool contains(const std::uint32_t cell, const Vector3d& x) const {
    return cell_contains_exact(topology_ + std::size_t(cell)*std::size_t(nv_), nv_,
                               coords_, stride_, x);
  }

  // Cells tested for a uniformly random point in the domain, on average
  double expected_cells_tested() const { return quality_; }
  int depth() const { return depth_; }
  std::size_t size() const { return ncells_; }
  const Vector3d& x_min() const { return x_min_; }
  const Vector3d& x_max() const { return x_max_; }
  const MortonBox& box() const { return box_; }
  const std::vector<std::uint32_t>& cell_list() const { return cells_; }
  const std::vector<CellTreeNode>& nodes() const { return nodes_; }
  // Every leaf as (first cell in cell_list, count)
  std::vector<std::pair<std::uint32_t, std::uint32_t>> leaves() const;

private:
  const std::uint32_t* topology_ = nullptr;
  const double* coords_ = nullptr;
  std::size_t ncells_ = 0, stride_ = 3;
  int nv_ = 4, dim_ = 3, depth_ = 0;
  double quality_ = 0.;
  Vector3d x_min_ = Vector3d::Zero(), x_max_ = Vector3d::Zero();
  MortonBox box_;
  std::vector<std::uint32_t> cells_;   // the caller's cell ids, in tree order
  std::vector<CellTreeNode> nodes_;
};

}  // namespace partrac

#endif
