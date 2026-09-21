// The dolfin-free cell tree, on hand-built structured meshes. What is checked
// is that the answer is the mesh's and nothing else's: the same containing
// cell as a scan over every cell with the same exact predicate, the lowest id
// where several cells contain the point (a facet, an edge, a vertex), and the
// same answer whatever the thread count the tree was built with. The quality
// number the build reports is recomputed here from the leaves, since it is the
// one number a later slow locate would be read against. The barycentrics are
// taken through Tet and Triangle, which is where every caller of the tree gets
// them: the tree only names the cell.
#include <catch2/catch.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>
#include <omp.h>

#include "typedefs.hpp"
#include "cell_tree.hpp"
#include "morton.hpp"
#include "Tet.hpp"
#include "Triangle.hpp"

using partrac::CellTree;
using partrac::MortonBox;
using partrac::MortonKey;

namespace {

// A structured box of simplices on a dyadic lattice, so that vertices, edge
// midpoints and the facet points below are exact doubles
struct BoxMesh {
  int nv = 4;
  std::size_t n = 8;
  std::vector<double> coords;        // gdim doubles a point, as the loaders store them
  std::vector<std::uint32_t> topo;
  std::size_t gdim() const { return std::size_t(nv) - 1; }
  std::size_t ncells() const { return topo.size()/std::size_t(nv); }
  std::size_t npoints() const { return coords.size()/gdim(); }
  const std::uint32_t* row(const std::size_t c) const { return topo.data() + c*std::size_t(nv); }
  const double* point(const std::uint32_t v) const { return coords.data() + std::size_t(v)*gdim(); }
  Vector3d vertex(const std::uint32_t v) const {
    const double* p = point(v);
    return Vector3d(p[0], p[1], nv == 4 ? p[2] : 0.);
  }
  Vector3d centroid(const std::size_t c) const {
    Vector3d m = Vector3d::Zero();
    for (int k = 0; k < nv; ++k) m += vertex(row(c)[k]);
    return m/nv;
  }
};

// The barycentrics of x in cell c, the way every caller of the tree gets them:
// from the cell itself, once the tree has named it
std::array<double, 4> barycentrics(const BoxMesh& m, const std::size_t c, const Vector3d& x){
  std::array<double, 4> bary{};
  const std::uint32_t* r = m.row(c);
  if (m.nv == 4) Tet(m.point(r[0]), m.point(r[1]), m.point(r[2]), m.point(r[3])).contains(x, bary);
  else           Triangle(m.point(r[0]), m.point(r[1]), m.point(r[2])).contains(x, bary);
  return bary;
}

// n x n x n cubes, each cut into the six tets of the Kuhn subdivision
BoxMesh tet_box(const std::size_t n){
  BoxMesh m;
  m.nv = 4;
  m.n = n;
  const double h = 1.0/double(n);
  const std::size_t s = n + 1;
  const auto id = [&](const std::size_t i, const std::size_t j, const std::size_t k){
    return std::uint32_t((i*s + j)*s + k);
  };
  m.coords.resize(3*s*s*s);
  for (std::size_t i = 0; i < s; ++i)
    for (std::size_t j = 0; j < s; ++j)
      for (std::size_t k = 0; k < s; ++k){
        const std::size_t v = id(i, j, k);
        m.coords[3*v] = double(i)*h;
        m.coords[3*v+1] = double(j)*h;
        m.coords[3*v+2] = double(k)*h;
      }
  const int perm[6][3] = {{0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0}};
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      for (std::size_t k = 0; k < n; ++k)
        for (int p = 0; p < 6; ++p){
          std::size_t c[3] = {i, j, k};
          m.topo.push_back(id(c[0], c[1], c[2]));
          for (int step = 0; step < 3; ++step){
            ++c[perm[p][step]];
            m.topo.push_back(id(c[0], c[1], c[2]));
          }
        }
  return m;
}

// n x n squares, each cut into two triangles
BoxMesh triangle_box(const std::size_t n){
  BoxMesh m;
  m.nv = 3;
  m.n = n;
  const double h = 1.0/double(n);
  const std::size_t s = n + 1;
  const auto id = [&](const std::size_t i, const std::size_t j){ return std::uint32_t(i*s + j); };
  m.coords.assign(2*s*s, 0.);
  for (std::size_t i = 0; i < s; ++i)
    for (std::size_t j = 0; j < s; ++j){
      const std::size_t v = id(i, j);
      m.coords[2*v] = double(i)*h;
      m.coords[2*v+1] = double(j)*h;
    }
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j){
      const std::uint32_t a = id(i, j), b = id(i+1, j), c = id(i+1, j+1), d = id(i, j+1);
      m.topo.insert(m.topo.end(), {a, b, c});
      m.topo.insert(m.topo.end(), {a, c, d});
    }
  return m;
}

CellTree build(const BoxMesh& m, const bool verbose = false){
  return CellTree(m.topo.data(), m.ncells(), m.nv, m.coords.data(), m.npoints(), m.gdim(), verbose);
}

// Every cell containing x, by the same exact test the tree uses
std::vector<int> containing(const BoxMesh& m, const Vector3d& x){
  std::vector<int> out;
  for (std::size_t c = 0; c < m.ncells(); ++c)
    if (partrac::cell_contains_exact(m.row(c), m.nv, m.coords.data(), m.gdim(), x))
      out.push_back(int(c));
  return out;
}

int lowest(const std::vector<int>& ids){ return ids.empty() ? -1 : ids.front(); }

}  // namespace

TEST_CASE("every cell centroid is located in that cell", "[celltree]") {
  const BoxMesh tets = tet_box(8);
  const CellTree tree = build(tets, true);   // the one place the quality line is printed
  for (std::size_t c = 0; c < tets.ncells(); ++c)
    REQUIRE( tree.locate(tets.centroid(c)) == int(c) );

  const BoxMesh tris = triangle_box(32);
  const CellTree tree2 = build(tris);
  for (std::size_t c = 0; c < tris.ncells(); ++c)
    REQUIRE( tree2.locate(tris.centroid(c)) == int(c) );
}

TEST_CASE("a random point lands in the same cell as a scan over every cell", "[celltree]") {
  std::mt19937_64 rng(20240920);
  std::uniform_real_distribution<double> u(0., 1.);

  const BoxMesh tets = tet_box(8);
  const CellTree tree = build(tets);
  for (int i = 0; i < 3000; ++i){
    const Vector3d x(u(rng), u(rng), u(rng));
    const int id = tree.locate(x);
    REQUIRE( id == lowest(containing(tets, x)) );
    REQUIRE( id >= 0 );                       // the box is filled by its tets
    REQUIRE( tree.contains(std::uint32_t(id), x) );
    // the barycentrics the evaluation reads, in floating point as before
    const std::array<double, 4> bary = barycentrics(tets, std::size_t(id), x);
    const double sum = bary[0] + bary[1] + bary[2] + bary[3];
    REQUIRE( sum == Approx(1.) );
    for (int k = 0; k < 4; ++k) REQUIRE( bary[k] > -1e-12 );
  }

  const BoxMesh tris = triangle_box(32);
  const CellTree tree2 = build(tris);
  for (int i = 0; i < 3000; ++i){
    const Vector3d x(u(rng), u(rng), 0.);
    const int id = tree2.locate(x);
    REQUIRE( id == lowest(containing(tris, x)) );
    REQUIRE( id >= 0 );
    const std::array<double, 4> bary = barycentrics(tris, std::size_t(id), x);
    REQUIRE( bary[0] + bary[1] + bary[2] == Approx(1.) );
    for (int k = 0; k < 3; ++k) REQUIRE( bary[k] > -1e-12 );
  }
}

TEST_CASE("a point on a shared facet, edge or vertex goes to the lowest cell id",
          "[celltree]") {
  const BoxMesh m = tet_box(8);
  const double h = 1.0/8.0;

  // A lattice vertex, an edge midpoint and a point in the interior of a facet,
  // all of them exact: the coordinates are eighths and the two halvings below
  // are exact for those.
  std::vector<Vector3d> points;
  points.push_back(Vector3d(4*h, 4*h, 4*h));                    // a vertex
  points.push_back(Vector3d(4.5*h, 4*h, 4*h));                  // an edge midpoint
  for (std::size_t c = 0; c < m.ncells() && points.size() < 3; ++c)
    for (int k = 0; k < 4 && points.size() < 3; ++k){
      // the facet opposite vertex k, halfway along a median of it
      std::array<Vector3d, 3> v;
      int nf = 0;
      for (int j = 0; j < 4; ++j) if (j != k) v[nf++] = m.vertex(m.row(c)[j]);
      const Vector3d p = 0.5*(0.5*(v[0] + v[1]) + v[2]);
      if (containing(m, p).size() == 2) points.push_back(p);    // an interior facet
    }
  REQUIRE( points.size() == 3 );
  REQUIRE( containing(m, points[0]).size() > 6 );               // shared by many
  REQUIRE( containing(m, points[1]).size() > 2 );
  REQUIRE( containing(m, points[2]).size() == 2 );

  // The answer is the mesh's, so it does not change with the thread count the
  // tree was built with, nor with how the traversal reached the cells
  const int threads_before = omp_get_max_threads();
  std::vector<std::array<int, 3>> answers;
  for (const int threads : {1, 4}){
    omp_set_num_threads(threads);
    const CellTree tree = build(m);
    std::array<int, 3> got{};
    for (std::size_t i = 0; i < 3; ++i) got[i] = tree.locate(points[i]);
    answers.push_back(got);
  }
  omp_set_num_threads(threads_before);
  for (std::size_t i = 0; i < 3; ++i){
    REQUIRE( answers[0][i] == lowest(containing(m, points[i])) );
    REQUIRE( answers[1][i] == answers[0][i] );
  }
}

TEST_CASE("the tree's quality is the expected number of cells tested", "[celltree]") {
  const BoxMesh m = tet_box(8);
  const CellTree tree = build(m);

  // The same sum, from the leaves, over the test's own boxes: the volume a
  // uniformly random point falls in, times the cells it then pays for
  double domain = 1.;
  for (int d = 0; d < 3; ++d) domain *= tree.x_max()[d] - tree.x_min()[d];
  double expected = 0.;
  std::size_t in_leaves = 0;
  for (const auto& leaf : tree.leaves()){
    Vector3d lo = Vector3d::Constant(HUGE_VAL), hi = Vector3d::Constant(-HUGE_VAL);
    for (std::uint32_t i = 0; i < leaf.second; ++i){
      const std::uint32_t c = tree.cell_list()[leaf.first + i];
      for (int k = 0; k < 4; ++k){
        const Vector3d v = m.vertex(m.row(c)[k]);
        lo = lo.cwiseMin(v);
        hi = hi.cwiseMax(v);
      }
    }
    expected += (hi - lo).prod()*leaf.second;
    in_leaves += leaf.second;
  }
  expected /= domain;
  REQUIRE( in_leaves == m.ncells() );          // every cell sits in exactly one leaf
  REQUIRE( tree.expected_cells_tested() == Approx(expected).epsilon(0.10) );
  // a tree worth having tests a leaf's worth of cells, not a mesh's
  REQUIRE( tree.expected_cells_tested() < 2*CellTree::leaf_cells );
  REQUIRE( tree.depth() > 0 );
}

TEST_CASE("a point outside every cell is not located", "[celltree]") {
  const BoxMesh tets = tet_box(8);
  const CellTree tree = build(tets);
  REQUIRE( tree.locate(Vector3d(5., 5., 5.)) == -1 );        // outside the tree's box
  REQUIRE( tree.locate(Vector3d(0.5, 0.5, -1e-9)) == -1 );   // just outside a facet
  REQUIRE( tree.locate(Vector3d(0.5, 0.5, 2.)) == -1 );      // above it, inside in x and y

  const BoxMesh tris = triangle_box(32);
  const CellTree tree2 = build(tris);
  REQUIRE( tree2.locate(Vector3d(-0.5, 0.5, 0.)) == -1 );
}

TEST_CASE("coordinates that are not the cell's own dimension are refused", "[celltree]") {
  // The tree's dimension is the cell type's, so the caller's coordinate stride
  // has to be that too: a triangle read out of three-wide points would take
  // every other pair, and nothing downstream would say so
  const BoxMesh tris = triangle_box(4);
  REQUIRE_THROWS_AS( CellTree(tris.topo.data(), tris.ncells(), tris.nv, tris.coords.data(),
                              tris.npoints(), 3), partrac::Error );
  const BoxMesh tets = tet_box(2);
  REQUIRE_THROWS_AS( CellTree(tets.topo.data(), tets.ncells(), tets.nv, tets.coords.data(),
                              tets.npoints(), 2), partrac::Error );
}

TEST_CASE("the Morton sort is deterministic where codes are equal", "[celltree]") {
  // A coarse lattice inside a fine box, so that thousands of points share a
  // code and the sort has to fall back on the index to order them
  const std::size_t n = 20000;
  std::mt19937_64 rng(7);
  std::uniform_int_distribution<int> lattice(0, 7);
  std::vector<double> x(3*n);
  for (std::size_t i = 0; i < n; ++i)
    for (int d = 0; d < 3; ++d) x[3*i+d] = 0.125*lattice(rng);
  const MortonBox box(Vector3d::Zero(), Vector3d::Ones(), 3);

  const int threads_before = omp_get_max_threads();
  omp_set_num_threads(1);
  const std::vector<std::uint32_t> one = partrac::morton_order(x.data(), n, 3, box);
  omp_set_num_threads(4);
  const std::vector<std::uint32_t> four = partrac::morton_order(x.data(), n, 3, box);
  omp_set_num_threads(threads_before);
  REQUIRE( one == four );

  std::size_t equal_codes = 0;
  for (std::size_t i = 1; i < n; ++i){
    const std::uint64_t a = box.code(x.data() + 3*one[i-1]), b = box.code(x.data() + 3*one[i]);
    REQUIRE( a <= b );
    if (a == b){
      ++equal_codes;
      REQUIRE( one[i-1] < one[i] );      // equal codes keep their index order
    }
  }
  REQUIRE( equal_codes > n/2 );          // the case this test is about
}
