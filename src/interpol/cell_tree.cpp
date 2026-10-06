#include "cell_tree.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>

namespace partrac {

// exactinit before any predicate call, once, before any thread runs
void init_exact_predicates(){
  static const bool once = [](){ exactinit(); return true; }();
  (void) once;
}

namespace {

constexpr float f_inf = std::numeric_limits<float>::infinity();

// float that brackets the double it came from
inline float f_up(const double v){ const float f = float(v); return f < v ? std::nextafterf(f, f_inf) : f; }
inline float f_down(const double v){ const float f = float(v); return f > v ? std::nextafterf(f, -f_inf) : f; }

// One of the ~16 candidates a node's split search classifies its cells into
struct Bucket {
  float lo[3], hi[3];
  std::uint32_t n;
  void clear(){ for (int d = 0; d < 3; ++d){ lo[d] = f_inf; hi[d] = -f_inf; } n = 0; }
  void add(const float* b){
    for (int d = 0; d < 3; ++d){
      if (b[d] < lo[d]) lo[d] = b[d];
      if (b[3+d] > hi[d]) hi[d] = b[3+d];
    }
    ++n;
  }
  void add(const Bucket& o){
    for (int d = 0; d < 3; ++d){
      if (o.lo[d] < lo[d]) lo[d] = o.lo[d];
      if (o.hi[d] > hi[d]) hi[d] = o.hi[d];
    }
    n += o.n;
  }
  double volume(const int dim) const {
    double v = 1.;
    for (int d = 0; d < dim; ++d) v *= double(hi[d]) - double(lo[d]);
    return v > 0. ? v : 0.;
  }
};

// The top-down build; every node covers a contiguous range of cells, so a
// split is two passes over the range and the subtrees are independent
struct Builder {
  CellTreeNode* nodes = nullptr;
  std::uint32_t* cells = nullptr;
  const float* bbox = nullptr;       // 6 floats a cell, by the caller's cell id
  std::uint32_t capacity = 0;
  int dim = 3;
  std::atomic<std::uint32_t> used{1};
  std::atomic<int> deepest{0};
  std::atomic<bool> full{false};

  const float* box_of(const std::uint32_t c) const { return bbox + std::size_t(c)*6; }
  float centre(const std::uint32_t c, const int d) const {
    const float* b = box_of(c);
    return 0.5f*(b[d] + b[3+d]);
  }
  void note_depth(const int d){
    int cur = deepest.load(std::memory_order_relaxed);
    while (d > cur && !deepest.compare_exchange_weak(cur, d, std::memory_order_relaxed)) {}
  }
  void make_leaf(const std::uint32_t node, const std::uint32_t first, const std::uint32_t n){
    nodes[node].axis = CellTreeNode::leaf_axis;
    nodes[node].a = first;
    nodes[node].n = n;
  }
  void split(std::uint32_t node, std::uint32_t first, std::uint32_t last, int depth);
};

void Builder::split(const std::uint32_t node, const std::uint32_t first,
                    const std::uint32_t last, const int depth){
  note_depth(depth);
  const std::uint32_t n = last - first;
  if (n <= CellTree::leaf_cells || depth >= CellTree::max_depth
      || full.load(std::memory_order_relaxed)){
    make_leaf(node, first, n);
    return;
  }
  // The range of the cells' box centres, which the buckets divide
  float clo[3] = {f_inf, f_inf, f_inf}, chi[3] = {-f_inf, -f_inf, -f_inf};
  for (std::uint32_t i = first; i < last; ++i){
    for (int d = 0; d < dim; ++d){
      const float c = centre(cells[i], d);
      if (c < clo[d]) clo[d] = c;
      if (c > chi[d]) chi[d] = c;
    }
  }
  float scale[3] = {0.f, 0.f, 0.f};
  for (int d = 0; d < dim; ++d)
    scale[d] = chi[d] > clo[d] ? CellTree::n_buckets/(chi[d] - clo[d]) : 0.f;
  const auto bucket_of = [&](const std::uint32_t c, const int d){
    int i = int((centre(c, d) - clo[d])*scale[d]);
    return i < 0 ? 0 : (i >= CellTree::n_buckets ? CellTree::n_buckets - 1 : i);
  };

  Bucket bk[3][CellTree::n_buckets];
  for (int d = 0; d < dim; ++d)
    for (int b = 0; b < CellTree::n_buckets; ++b) bk[d][b].clear();
  for (std::uint32_t i = first; i < last; ++i){
    const float* b = box_of(cells[i]);
    for (int d = 0; d < dim; ++d) bk[d][bucket_of(cells[i], d)].add(b);
  }

  // vol(L)*N_L + vol(R)*N_R, the cells a random point in the node pays for
  int axis = -1, cut = 0;
  double best = std::numeric_limits<double>::infinity();
  for (int d = 0; d < dim; ++d){
    if (scale[d] == 0.f) continue;
    Bucket suf[CellTree::n_buckets + 1];
    suf[CellTree::n_buckets].clear();
    for (int b = CellTree::n_buckets - 1; b >= 0; --b){
      suf[b] = suf[b+1];
      suf[b].add(bk[d][b]);
    }
    Bucket pre;
    pre.clear();
    for (int b = 0; b < CellTree::n_buckets - 1; ++b){
      pre.add(bk[d][b]);
      if (pre.n == 0 || suf[b+1].n == 0) continue;
      const double cost = pre.volume(dim)*pre.n + suf[b+1].volume(dim)*suf[b+1].n;
      if (cost < best){ best = cost; axis = d; cut = b + 1; }
    }
  }

  std::uint32_t nl = 0;
  if (axis >= 0){
    const int c = cut, d = axis;
    std::uint32_t* mid = std::partition(cells + first, cells + last,
                                        [&](const std::uint32_t id){ return bucket_of(id, d) < c; });
    nl = std::uint32_t(mid - (cells + first));
  }
  if (axis < 0 || nl == 0 || nl == n){
    // Every cell in one bucket: the median along the widest axis, keeping the
    // order total so the two halves do not depend on the input order
    axis = 0;
    for (int d = 1; d < dim; ++d)
      if (chi[d] - clo[d] > chi[axis] - clo[axis]) axis = d;
    const int d = axis;
    nl = n/2;
    std::nth_element(cells + first, cells + first + nl, cells + last,
                     [&](const std::uint32_t a, const std::uint32_t b){
                       const float ca = centre(a, d), cb = centre(b, d);
                       return ca != cb ? ca < cb : a < b;
                     });
  }

  // The planes: the left part's largest and the right part's smallest along
  // the axis, each padded outward by one ulp so no cell is lost to rounding
  float lmax = -f_inf, rmin = f_inf;
  for (std::uint32_t i = first; i < first + nl; ++i)
    lmax = std::max(lmax, box_of(cells[i])[3+axis]);
  for (std::uint32_t i = first + nl; i < last; ++i)
    rmin = std::min(rmin, box_of(cells[i])[axis]);

  const std::uint32_t kid = used.fetch_add(2);
  if (std::size_t(kid) + 1 >= std::size_t(capacity)){
    full.store(true, std::memory_order_relaxed);
    make_leaf(node, first, n);
    return;
  }
  nodes[node].axis = std::uint8_t(axis);
  nodes[node].a = kid;
  nodes[node].n = 0;
  nodes[node].lmax = std::nextafterf(lmax, f_inf);
  nodes[node].rmin = std::nextafterf(rmin, -f_inf);

  const std::uint32_t mid = first + nl;
#pragma omp task default(shared) firstprivate(kid, first, mid, depth) if (mid - first >= CellTree::task_cells)
  split(kid, first, mid, depth + 1);
  split(kid + 1, mid, last, depth + 1);
}

}  // namespace

CellTree::CellTree(const std::uint32_t* topology, const std::size_t ncells, const int nv,
                   const double* coords, const std::size_t npoints,
                   const std::size_t coord_stride, const bool verbose)
  : topology_(topology), coords_(coords), ncells_(ncells), stride_(coord_stride),
    nv_(nv), dim_(nv - 1)
{
  init_exact_predicates();
  if (nv != 3 && nv != 4)
    fail("cell tree: ", nv, " vertices a cell, which is neither a triangle nor a tet");
  // The dimension is the cell's, so the coordinates have to come that wide
  if (coord_stride != std::size_t(dim_))
    fail("cell tree: a cell of ", nv, " vertices is ", dim_, "D, but the coordinates come ",
         coord_stride, " to a point");
  if (ncells > std::size_t(std::numeric_limits<std::int32_t>::max()))
    fail("cell tree: ", ncells, " cells do not fit the 32-bit cell id of a position");
  if (ncells == 0 || topology == nullptr || coords == nullptr) return;

  std::uint32_t vmax = 0;
#pragma omp parallel for reduction(max:vmax) schedule(static)
  for (std::ptrdiff_t i = 0; i < std::ptrdiff_t(ncells*std::size_t(nv)); ++i)
    vmax = std::max(vmax, topology[i]);
  if (std::size_t(vmax) >= npoints)
    fail("cell tree: vertex ", vmax, " of the topology is outside ", npoints, " coordinates");

  // Every cell's box, and the domain's
  std::vector<float> bbox(6*ncells, 0.f);
  double lo0 = HUGE_VAL, lo1 = HUGE_VAL, lo2 = HUGE_VAL;
  double hi0 = -HUGE_VAL, hi1 = -HUGE_VAL, hi2 = -HUGE_VAL;
#pragma omp parallel for reduction(min:lo0,lo1,lo2) reduction(max:hi0,hi1,hi2) schedule(static)
  for (std::ptrdiff_t c = 0; c < std::ptrdiff_t(ncells); ++c){
    const std::uint32_t* row = topology + std::size_t(c)*std::size_t(nv);
    double lo[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL}, hi[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
    for (int k = 0; k < nv; ++k){
      const double* v = coords + std::size_t(row[k])*coord_stride;
      for (int d = 0; d < dim_; ++d){
        lo[d] = std::min(lo[d], v[d]);
        hi[d] = std::max(hi[d], v[d]);
      }
    }
    float* b = bbox.data() + std::size_t(c)*6;
    for (int d = 0; d < dim_; ++d){
      b[d] = f_down(lo[d]);
      b[3+d] = f_up(hi[d]);
    }
    lo0 = std::min(lo0, lo[0]);  hi0 = std::max(hi0, hi[0]);
    lo1 = std::min(lo1, lo[1]);  hi1 = std::max(hi1, hi[1]);
    if (dim_ == 3){ lo2 = std::min(lo2, lo[2]);  hi2 = std::max(hi2, hi[2]); }
  }
  x_min_ = Vector3d(lo0, lo1, dim_ == 3 ? lo2 : 0.);
  x_max_ = Vector3d(hi0, hi1, dim_ == 3 ? hi2 : 0.);
  box_ = MortonBox(x_min_, x_max_, dim_);

  // Cells in Morton order of their centroid, so every node's range is contiguous
  std::vector<MortonKey> keys(ncells);
#pragma omp parallel for schedule(static)
  for (std::ptrdiff_t c = 0; c < std::ptrdiff_t(ncells); ++c){
    const std::uint32_t* row = topology + std::size_t(c)*std::size_t(nv);
    double m[3] = {0., 0., 0.};
    for (int k = 0; k < nv; ++k){
      const double* v = coords + std::size_t(row[k])*coord_stride;
      for (int d = 0; d < dim_; ++d) m[d] += v[d];
    }
    for (int d = 0; d < dim_; ++d) m[d] /= nv;
    keys[c].code = box_.code(m);
    keys[c].index = std::uint32_t(c);
  }
  cells_ = morton_permutation(keys);
  keys.clear();
  keys.shrink_to_fit();

  // Room for a leaf per two cells; a build that runs out stops splitting
  const std::size_t capacity = std::max<std::size_t>(64, ncells/2 + 8);
  nodes_.assign(capacity, CellTreeNode());
  Builder b;
  b.nodes = nodes_.data();
  b.cells = cells_.data();
  b.bbox = bbox.data();
  b.capacity = std::uint32_t(capacity);
  b.dim = dim_;
#pragma omp parallel
#pragma omp single
  b.split(0, 0, std::uint32_t(ncells), 0);
  const std::size_t used = std::min<std::size_t>(b.used.load(), capacity);
  nodes_.resize(used);
  nodes_.shrink_to_fit();
  depth_ = b.deepest.load();

  // Cells tested for a uniformly random point: the leaves the point can fall in
  double domain = 1.;
  for (int d = 0; d < dim_; ++d) domain *= x_max_[d] - x_min_[d];
  std::size_t nleaves = 0;
  double q = 0.;
  for (const CellTreeNode& nd : nodes_){
    if (!nd.is_leaf()) continue;
    ++nleaves;
    Bucket leaf;
    leaf.clear();
    for (std::uint32_t i = 0; i < nd.n; ++i) leaf.add(bbox.data() + std::size_t(cells_[nd.a + i])*6);
    q += leaf.volume(dim_)*nd.n;
  }
  quality_ = domain > 0. ? q/domain : 0.;

  // A build that ran out of nodes left its deepest ranges unsplit, which a
  // locate pays for: said whether or not the statistics are asked for
  if (b.full.load())
    std::cerr << "Cell tree: the node pool of " << capacity << " ran out on " << ncells_
              << " cells; the deepest ranges are scanned" << std::endl;
  if (verbose){
    std::cout << "Cell tree: " << ncells_ << " cells, " << used << " nodes, " << nleaves
              << " leaves, depth " << depth_ << ", cells tested per random point "
              << quality_ << std::endl;
  }
}

std::vector<std::pair<std::uint32_t, std::uint32_t>> CellTree::leaves() const {
  std::vector<std::pair<std::uint32_t, std::uint32_t>> out;
  for (const CellTreeNode& nd : nodes_)
    if (nd.is_leaf() && nd.n > 0) out.emplace_back(nd.a, nd.n);
  return out;
}

int CellTree::locate(const Vector3d& x, std::size_t* tested) const {
  int best = -1;
  std::size_t ntest = 0;
  if (!nodes_.empty()){
    std::uint32_t stack[max_depth + 2];
    int top = 0;
    stack[top++] = 0;
    while (top > 0){
      const CellTreeNode& nd = nodes_[stack[--top]];
      if (nd.is_leaf()){
        for (std::uint32_t i = 0; i < nd.n; ++i){
          const std::uint32_t c = cells_[nd.a + i];
          ++ntest;
          if (cell_contains_exact(topology_ + std::size_t(c)*std::size_t(nv_), nv_,
                                  coords_, stride_, x)
              && (best < 0 || int(c) < best))
            best = int(c);
        }
      }
      else {
        const double v = x[nd.axis];
        if (v <= double(nd.lmax)) stack[top++] = nd.a;
        if (v >= double(nd.rmin)) stack[top++] = nd.a + 1;
      }
    }
  }
  if (tested) *tested = ntest;
  return best;
}

}  // namespace partrac
