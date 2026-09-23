// The dolfin-free mesh tables: the parallel key sort, the edge numbering in
// dolfin's local order, the grouped dof-to-node scatter and the facet table
// with its periodic partners. The fixtures are structured boxes built here, so
// the reference is independent of the code under test: an O(N^2) facet match
// written in this file, and -- where the build has dolfin -- the facet table
// the existing loaders get from build_facet_neighbours on the same mesh. A
// mistake in any of these tables is silent in a run: a particle lands in the
// wrong cell, or reads a neighbour's dof.
#include <catch2/catch.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <random>
#include <set>
#include <vector>
#include <omp.h>

#include "typedefs.hpp"
#include "Error.hpp"
#include "mesh_tables.hpp"

#ifdef USE_DOLFIN
#include <dolfin.h>
#include "Tet.hpp"
#include "Triangle.hpp"
#include "dolfin_ref.hpp"
#endif

namespace {

using mesh_tables::local_edges;
using mesh_tables::n_local_edges;

// A mesh as the tables take it: the topology, the coordinates and the box
template<int NV>
struct Box {
  std::size_t n = 0;
  std::vector<std::uint32_t> topo;
  std::vector<double> coords;
  Vector3d lo = Vector3d::Zero(), hi = Vector3d::Zero();
  std::size_t ncells() const { return topo.size()/NV; }
  std::size_t nverts() const { return coords.size()/(NV-1); }
};

// n x n x n cubes on the unit cube, each cut into the six Kuhn tets (the
// tetrahedra 0 -> e_a -> e_a+e_b -> 1 over the six axis orders), which agree on
// every shared face, so the mesh is conforming and its outer faces pair up
// across the box
Box<4> tet_box(const std::size_t n){
  Box<4> b;
  b.n = n;
  b.hi = Vector3d(1., 1., 1.);
  const std::size_t m = n + 1;
  b.coords.resize(3*m*m*m);
  auto vid = [m](const std::size_t i, const std::size_t j, const std::size_t k){
    return (i*m + j)*m + k;
  };
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t j = 0; j < m; ++j)
      for (std::size_t k = 0; k < m; ++k){
        double* c = b.coords.data() + 3*vid(i, j, k);
        c[0] = double(i)/double(n);
        c[1] = double(j)/double(n);
        c[2] = double(k)/double(n);
      }
  std::array<int, 3> ax = {0, 1, 2};
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      for (std::size_t k = 0; k < n; ++k){
        std::array<int, 3> p = ax;
        std::sort(p.begin(), p.end());
        do {
          std::array<int, 3> o = {0, 0, 0};
          for (int s = 0; s < 4; ++s){
            b.topo.push_back(std::uint32_t(vid(i + std::size_t(o[0]), j + std::size_t(o[1]),
                                               k + std::size_t(o[2]))));
            if (s < 3) o[p[s]] = 1;
          }
        } while (std::next_permutation(p.begin(), p.end()));
      }
  return b;
}

// n x n squares on the unit square, each cut along its (0,0)-(1,1) diagonal
Box<3> triangle_box(const std::size_t n){
  Box<3> b;
  b.n = n;
  b.hi = Vector3d(1., 1., 0.);
  const std::size_t m = n + 1;
  b.coords.resize(2*m*m);
  auto vid = [m](const std::size_t i, const std::size_t j){ return i*m + j; };
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t j = 0; j < m; ++j){
      double* c = b.coords.data() + 2*vid(i, j);
      c[0] = double(i)/double(n);
      c[1] = double(j)/double(n);
    }
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j){
      const std::uint32_t v00 = std::uint32_t(vid(i, j)), v10 = std::uint32_t(vid(i+1, j));
      const std::uint32_t v11 = std::uint32_t(vid(i+1, j+1)), v01 = std::uint32_t(vid(i, j+1));
      b.topo.insert(b.topo.end(), {v00, v10, v11});
      b.topo.insert(b.topo.end(), {v00, v11, v01});
    }
  return b;
}

// The facet table the slow way: every pair of facets compared, then every pair
// of exterior midpoints along each periodic axis
template<int NV>
std::vector<std::int32_t> brute_facets(const Box<NV>& b, const std::vector<bool>& periodic,
                                       const double tol){
  const Uint gdim = NV - 1;
  const std::size_t nc = b.ncells();
  std::vector<std::int32_t> across(nc*NV, mesh_tables::facet_wall);
  // Sorted vertices of each facet, and its midpoint
  std::vector<std::array<std::uint32_t, NV-1>> fv(nc*NV);
  std::vector<Vector3d> mid(nc*NV, Vector3d::Zero());
  for (std::size_t i = 0; i < nc; ++i)
    for (std::size_t k = 0; k < NV; ++k){
      auto& w = fv[i*NV + k];
      int m = 0;
      for (std::size_t j = 0; j < NV; ++j)
        if (j != k) w[m++] = b.topo[i*NV + j];
      std::sort(w.begin(), w.end());
      for (std::size_t j = 0; j < NV-1; ++j)
        for (Uint d = 0; d < gdim; ++d) mid[i*NV + k][d] += b.coords[std::size_t(w[j])*gdim + d];
      mid[i*NV + k] /= double(NV - 1);
    }
  for (std::size_t a = 0; a < fv.size(); ++a)
    for (std::size_t c = a + 1; c < fv.size(); ++c)
      if (fv[a] == fv[c]){
        across[a] = std::int32_t(c/NV);
        across[c] = std::int32_t(a/NV);
      }
  // Exterior facets, each on the first periodic axis whose face it lies on
  std::vector<std::vector<std::size_t>> lo(gdim), hi(gdim);
  for (std::size_t s = 0; s < fv.size(); ++s){
    if (across[s] != mesh_tables::facet_wall) continue;
    for (Uint a = 0; a < gdim; ++a){
      if (!periodic[a]) continue;
      if (mid[s][a] < b.lo[a] + tol){ lo[a].push_back(s); break; }
      if (mid[s][a] > b.hi[a] - tol){ hi[a].push_back(s); break; }
    }
  }
  for (Uint a = 0; a < gdim; ++a)
    for (const std::size_t l : lo[a]){
      Vector3d q = mid[l];
      q[a] += b.hi[a] - b.lo[a];
      for (const std::size_t h : hi[a])
        if ((q - mid[h]).norm() < tol){
          across[l] = mesh_tables::facet_periodic(std::int32_t(h/NV));
          across[h] = mesh_tables::facet_periodic(std::int32_t(l/NV));
          break;
        }
    }
  return across;
}

// The unique edges of the topology, counted with a set
template<int NV>
std::size_t brute_edge_count(const Box<NV>& b){
  std::set<std::pair<std::uint32_t, std::uint32_t>> s;
  for (std::size_t i = 0; i < b.ncells(); ++i)
    for (std::size_t a = 0; a < NV; ++a)
      for (std::size_t c = a + 1; c < NV; ++c){
        const std::uint32_t p = b.topo[i*NV + a], q = b.topo[i*NV + c];
        s.insert({std::min(p, q), std::max(p, q)});
      }
  return s.size();
}

}  // namespace

TEST_CASE("The local edge order is dolfin's, the one quadbasis expects", "[mesh_tables]") {
  // dolfin's order: the vertex pairs a < b, reverse lexicographic
  constexpr auto et = local_edges<4>();
  REQUIRE(et == std::array<std::array<int, 2>, 6>{{{2,3},{1,3},{1,2},{0,3},{0,2},{0,1}}});
  constexpr auto tr = local_edges<3>();
  REQUIRE(tr == std::array<std::array<int, 2>, 3>{{{1,2},{0,2},{0,1}}});
#ifdef USE_DOLFIN
  // Tet::mid_ holds the quadbasis slot of the midpoint of edges 01, 02, 03, 12,
  // 13, 23: edge e of the table above must be the slot 4 + e the gather fills
  const std::array<std::array<int, 2>, 6> asc = {{{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}}};
  for (std::size_t i = 0; i < 6; ++i){
    const std::size_t e = std::size_t(std::find(et.begin(), et.end(), asc[i]) - et.begin());
    REQUIRE(Tet::mid_[i] == int(4 + e));
  }
  const std::array<std::array<int, 2>, 3> asc2 = {{{0,1},{0,2},{1,2}}};
  for (std::size_t i = 0; i < 3; ++i){
    const std::size_t e = std::size_t(std::find(tr.begin(), tr.end(), asc2[i]) - tr.begin());
    REQUIRE(Triangle::mid_[i] == int(3 + e));
  }
#endif
}

TEST_CASE("The tet facet table matches an O(N^2) match on a periodic box", "[mesh_tables]") {
  const Box<4> b = tet_box(3);
  const std::vector<bool> periodic = {true, true, true};
  std::vector<std::int32_t> across;
  mesh_tables::build_facet_neighbours<4>(b.topo, b.ncells(), b.coords, 3, periodic,
                                         b.lo, b.hi, 1e-12, across);
  REQUIRE(across == brute_facets<4>(b, periodic, 1e-12));
  // Nothing is left a wall when every axis is periodic, and the table is symmetric
  for (std::size_t s = 0; s < across.size(); ++s){
    REQUIRE(across[s] != mesh_tables::facet_wall);
    const std::int32_t a = across[s];
    const std::size_t other = std::size_t(a >= 0 ? a : mesh_tables::facet_periodic(a));
    bool back = false;
    for (std::size_t k = 0; k < 4; ++k){
      const std::int32_t c = across[other*4 + k];
      const std::size_t j = std::size_t(c >= 0 ? c : mesh_tables::facet_periodic(c));
      if (j == s/4 && (c >= 0) == (a >= 0)) back = true;
    }
    REQUIRE(back);
  }
  // With walls instead, the six faces of the box are exterior
  std::vector<std::int32_t> walls;
  mesh_tables::build_facet_neighbours<4>(b.topo, b.ncells(), b.coords, 3, {false, false, false},
                                         b.lo, b.hi, 1e-12, walls);
  REQUIRE(walls == brute_facets<4>(b, {false, false, false}, 1e-12));
  REQUIRE(std::count(walls.begin(), walls.end(), mesh_tables::facet_wall) == 12*3*3);
}

TEST_CASE("The triangle facet table matches an O(N^2) match on a periodic box", "[mesh_tables]") {
  const Box<3> b = triangle_box(4);
  for (const auto& periodic : std::vector<std::vector<bool>>{{true, true, false},
                                                             {true, false, false},
                                                             {false, false, false}}){
    std::vector<std::int32_t> across;
    mesh_tables::build_facet_neighbours<3>(b.topo, b.ncells(), b.coords, 2, periodic,
                                           b.lo, b.hi, 1e-12, across);
    REQUIRE(across == brute_facets<3>(b, periodic, 1e-12));
  }
}

#ifdef USE_DOLFIN
TEST_CASE("The facet table equals the one the dolfin loaders build", "[mesh_tables]") {
  // The same mesh both ways: the topology handed to the tables is read back out
  // of dolfin, so a cell's vertices are in the order that decides the slots
  const std::vector<bool> periodic = {true, true, true};
  const Vector3d lo(0., 0., 0.), hi(1., 1., 1.);
  auto mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitCubeMesh(3, 3, 3));
  mesh->init();
  std::vector<dolfin::Cell> dolfin_cells;
  for (dolfin::CellIterator c(*mesh); !c.end(); ++c) dolfin_cells.push_back(*c);
  std::vector<std::int32_t> ref;
  build_facet_neighbours(ref, mesh, dolfin_cells, nullptr, periodic, lo, hi, 3, 1e-12);

  std::vector<std::uint32_t> topo;
  for (const auto& c : dolfin_cells){
    const auto v = c.entities(0);
    for (std::size_t k = 0; k < 4; ++k) topo.push_back(std::uint32_t(v[k]));
  }
  std::vector<std::int32_t> mine;
  mesh_tables::build_facet_neighbours<4>(topo, dolfin_cells.size(), mesh->coordinates(), 3,
                                         periodic, lo, hi, 1e-12, mine);
  REQUIRE(mine == ref);
}

TEST_CASE("The triangle facet table equals the one the dolfin loaders build", "[mesh_tables]") {
  const std::vector<bool> periodic = {true, true, false};
  const Vector3d lo(0., 0., 0.), hi(1., 1., 0.);
  auto mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitSquareMesh(4, 4));
  mesh->init();
  std::vector<dolfin::Cell> dolfin_cells;
  for (dolfin::CellIterator c(*mesh); !c.end(); ++c) dolfin_cells.push_back(*c);
  std::vector<std::int32_t> ref;
  build_facet_neighbours(ref, mesh, dolfin_cells, nullptr, periodic, lo, hi, 2, 1e-12);

  std::vector<std::uint32_t> topo;
  for (const auto& c : dolfin_cells){
    const auto v = c.entities(0);
    for (std::size_t k = 0; k < 3; ++k) topo.push_back(std::uint32_t(v[k]));
  }
  std::vector<std::int32_t> mine;
  mesh_tables::build_facet_neighbours<3>(topo, dolfin_cells.size(), mesh->coordinates(), 2,
                                         periodic, lo, hi, 1e-12, mine);
  REQUIRE(mine == ref);
}
#endif

TEST_CASE("Every cell gets ten distinct nodes and the edges are counted once", "[mesh_tables]") {
  const Box<4> b = tet_box(3);
  std::vector<std::uint32_t> edges;
  const std::size_t ne = mesh_tables::build_edge_table<4>(b.topo, b.ncells(), edges);
  REQUIRE(ne == brute_edge_count<4>(b));
  REQUIRE(edges.size() == 6*b.ncells());
  constexpr auto loc = local_edges<4>();
  std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> seen;
  for (std::size_t i = 0; i < b.ncells(); ++i){
    // The ten P2 nodes of the cell: four vertices and six edges
    std::set<std::size_t> nodes;
    for (std::size_t j = 0; j < 4; ++j) nodes.insert(b.topo[i*4 + j]);
    for (std::size_t e = 0; e < 6; ++e) nodes.insert(b.nverts() + edges[i*6 + e]);
    REQUIRE(nodes.size() == 10);
    // Column e is the edge of the vertex pair local_edges[e], the same id everywhere
    for (std::size_t e = 0; e < 6; ++e){
      const std::uint32_t p = b.topo[i*4 + std::size_t(loc[e][0])];
      const std::uint32_t q = b.topo[i*4 + std::size_t(loc[e][1])];
      const auto key = std::make_pair(std::min(p, q), std::max(p, q));
      const auto it = seen.find(key);
      if (it == seen.end()) seen[key] = edges[i*6 + e];
      else REQUIRE(it->second == edges[i*6 + e]);
    }
  }
  REQUIRE(seen.size() == ne);
  // Distinct pairs get distinct ids
  std::set<std::uint32_t> ids;
  for (const auto& kv : seen) ids.insert(kv.second);
  REQUIRE(ids.size() == ne);

  // The triangle case is the same code with three columns
  const Box<3> t = triangle_box(4);
  std::vector<std::uint32_t> tedges;
  const std::size_t tne = mesh_tables::build_edge_table<3>(t.topo, t.ncells(), tedges);
  REQUIRE(tne == brute_edge_count<3>(t));
  REQUIRE(tedges.size() == 3*t.ncells());
}

namespace {

// A stored P2 dof table for the box: the nodes are numbered in some order of
// their own (reversed here, as a solver's would not be ours), a node's
// components are adjacent, and the columns are blocked by component
struct Stored {
  std::size_t n_nodes = 0, ncomp = 0, n_total = 0;
  std::vector<std::uint32_t> dofs;
  std::vector<double> vec;
  std::vector<std::size_t> stored_of;   // our node -> the file's node
};

double node_value(const std::size_t node, const std::size_t c){
  return 1. + 0.25*double(node) - 3.5*double(c);
}

Stored make_stored(const Box<4>& b, const std::vector<std::uint32_t>& edges,
                   const std::size_t nedges, const std::size_t ncomp, const bool quadratic){
  Stored s;
  s.n_nodes = quadratic ? 10 : 4;
  s.ncomp = ncomp;
  s.n_total = b.nverts() + (quadratic ? nedges : 0);
  s.stored_of.resize(s.n_total);
  for (std::size_t i = 0; i < s.n_total; ++i) s.stored_of[i] = s.n_total - 1 - i;
  s.vec.assign(s.n_total*ncomp, 0.);
  for (std::size_t i = 0; i < s.n_total; ++i)
    for (std::size_t c = 0; c < ncomp; ++c)
      s.vec[s.stored_of[i]*ncomp + c] = node_value(i, c);
  s.dofs.assign(b.ncells()*s.n_nodes*ncomp, 0);
  for (std::size_t i = 0; i < b.ncells(); ++i)
    for (std::size_t j = 0; j < s.n_nodes; ++j){
      const std::size_t node = j < 4 ? b.topo[i*4 + j] : b.nverts() + edges[i*6 + (j - 4)];
      for (std::size_t c = 0; c < ncomp; ++c)
        s.dofs[(i*ncomp + c)*s.n_nodes + j] = std::uint32_t(s.stored_of[node]*ncomp + c);
    }
  return s;
}

}  // namespace

TEST_CASE("The grouped scatter gives every node one value, or fails", "[mesh_tables]") {
  const Box<4> b = tet_box(2);
  std::vector<std::uint32_t> edges;
  const std::size_t nedges = mesh_tables::build_edge_table<4>(b.topo, b.ncells(), edges);
  const std::vector<std::uint32_t> no_edges;

  SECTION("P2 vector: the values come back in vertex-then-edge order") {
    const Stored s = make_stored(b, edges, nedges, 3, true);
    std::vector<double> values;
    mesh_tables::DofNodes map;
    mesh_tables::scatter_dofs_to_nodes<4>(b.topo, edges, b.ncells(), b.nverts(), nedges, 3,
                                          s.dofs, s.vec, values, map);
    REQUIRE(values.size() == s.n_total*3);
    for (std::size_t i = 0; i < s.n_total; ++i)
      for (std::size_t c = 0; c < 3; ++c){
        REQUIRE(values[i*3 + c] == node_value(i, c));
        // The mapping moves a later stamp's vector without the dof table; here
        // the space is unconstrained, so each dof feeds exactly one node slot
        const std::size_t d = s.stored_of[i]*3 + c;
        REQUIRE(map.start[d + 1] - map.start[d] == 1);
        REQUIRE(map.slot[map.start[d]] == i*3 + c);
      }
  }

  SECTION("P1 scalar: the same code with four nodes and one component") {
    const Stored s = make_stored(b, edges, nedges, 1, false);
    std::vector<double> values;
    mesh_tables::DofNodes map;
    mesh_tables::scatter_dofs_to_nodes<4>(b.topo, no_edges, b.ncells(), b.nverts(), 0, 1,
                                          s.dofs, s.vec, values, map);
    REQUIRE(values.size() == b.nverts());
    for (std::size_t i = 0; i < b.nverts(); ++i) REQUIRE(values[i] == node_value(i, 0));
  }

  SECTION("A cell pointed at another node's dof is a disagreement, not a last writer") {
    Stored s = make_stored(b, edges, nedges, 3, true);
    // A vertex of the last cell, redirected to the dof of another node: the
    // cells sharing that vertex then disagree
    const std::size_t i = b.ncells() - 1;
    const std::size_t node = b.topo[i*4];
    const std::size_t other = (node + 1) % s.n_total;
    s.dofs[(i*3 + 0)*10 + 0] = std::uint32_t(s.stored_of[other]*3 + 0);
    std::vector<double> values;
    mesh_tables::DofNodes map;
    REQUIRE_THROWS_AS(mesh_tables::scatter_dofs_to_nodes<4>(b.topo, edges, b.ncells(), b.nverts(),
                                                            nedges, 3, s.dofs, s.vec, values, map),
                      partrac::Error);
  }

  SECTION("A dof outside the vector fails rather than reading past it") {
    Stored s = make_stored(b, edges, nedges, 3, true);
    s.dofs[0] = std::uint32_t(s.vec.size());
    std::vector<double> values;
    mesh_tables::DofNodes map;
    REQUIRE_THROWS_AS(mesh_tables::scatter_dofs_to_nodes<4>(b.topo, edges, b.ncells(), b.nverts(),
                                                            nedges, 3, s.dofs, s.vec, values, map),
                      partrac::Error);
  }
}

TEST_CASE("A dof shared by a node and its image feeds both", "[mesh_tables]") {
  // The smallest periodic-reduced case there is: one P1 tet whose vertices 0
  // and 3 are a periodic pair, so the space stores three dofs for four nodes
  // and dof 0 serves two of them. A dof -> node permutation can hold only one
  // of the two, and every later stamp then leaves the other where it was.
  const std::vector<std::uint32_t> topo = {0, 1, 2, 3};
  const std::vector<std::uint32_t> no_edges;
  const std::vector<std::uint32_t> cell_dofs = {0, 1, 2, 0};
  const std::vector<double> vec = {5., 6., 7.};
  std::vector<double> values;
  mesh_tables::DofNodes map;
  mesh_tables::scatter_dofs_to_nodes<4>(topo, no_edges, 1, 4, 0, 1, cell_dofs, vec, values, map);
  REQUIRE(values == std::vector<double>{5., 6., 7., 5.});

  REQUIRE(map.start.size() == vec.size() + 1);
  REQUIRE(map.slot.size() == 4);
  REQUIRE(map.start[1] - map.start[0] == 2);
  REQUIRE(map.slot[map.start[0]] == 0);
  REQUIRE(map.slot[map.start[0] + 1] == 3);
  REQUIRE(map.start[2] - map.start[1] == 1);
  REQUIRE(map.start[3] - map.start[2] == 1);

  // A second stamp read through the mapping, as simplex_load::read_vector does:
  // the image node must move with its master, not keep the first stamp
  const std::vector<double> later = {50., 60., 70.};
  for (std::size_t i = 0; i < later.size(); ++i)
    for (std::size_t q = map.start[i]; q < map.start[i + 1]; ++q)
      values[map.slot[q]] = later[i];
  REQUIRE(values == std::vector<double>{50., 60., 70., 50.});
}

TEST_CASE("The key sort is stable and gives the same order at 1 and 4 threads", "[mesh_tables]") {
  // Many equal keys: the thread count must not decide which of them comes first
  const std::size_t n = 200000;
  std::vector<std::uint64_t> base(n);
  std::mt19937_64 rng(7);
  for (std::size_t i = 0; i < n; ++i) base[i] = rng() % 97;
  const int saved = omp_get_max_threads();
  std::vector<std::vector<std::uint32_t>> out;
  for (const int nt : {1, 4}){
    omp_set_num_threads(nt);
    std::vector<std::uint64_t> key = base;
    std::vector<std::uint32_t> pay;
    mesh_tables::sort_by_key_indexed(key, pay);
    bool ordered = true, carried = true;
    for (std::size_t i = 0; i < n; ++i){
      if (key[i] != base[pay[i]]) carried = false;
      // Ties in the order of the input index, so the order is total
      if (i && (key[i-1] > key[i] || (key[i-1] == key[i] && pay[i-1] >= pay[i]))) ordered = false;
    }
    REQUIRE(carried);
    REQUIRE(ordered);
    out.push_back(pay);
  }
  omp_set_num_threads(saved);
  REQUIRE(out[0] == out[1]);
}

TEST_CASE("A vertex on a max face is matched to its image on the min face", "[mesh_tables]") {
  // A P1 field written from a constrained space gives a vertex and its images
  // one value, so the XDMF loader reads them through this matching; a corner
  // of a doubly periodic box has to reach the master through the chain
  const std::size_t n = 4;
  const double h = 1./double(n);
  std::vector<double> coords;
  for (std::size_t i = 0; i <= n; ++i)
    for (std::size_t j = 0; j <= n; ++j){
      coords.push_back(double(j)*h);
      coords.push_back(double(i)*h);
    }
  const std::size_t nverts = coords.size()/2;
  const Vector3d x_min(0., 0., 0.), x_max(1., 1., 0.);
  const auto at = [&](const std::size_t i, const std::size_t j){ return i*(n+1) + j; };

  SECTION("periodic in x only"){
    const std::vector<bool> per = {true, false, false};
    const std::vector<std::uint32_t> master =
      mesh_tables::match_periodic_vertices(coords, nverts, 2, per, x_min, x_max, 1e-12);
    for (std::size_t i = 0; i <= n; ++i){
      REQUIRE(master[at(i, n)] == at(i, 0));      // the max face reads the min face
      REQUIRE(master[at(i, 0)] == at(i, 0));      // a master is its own
      REQUIRE(master[at(i, 1)] == at(i, 1));      // the inside is untouched
    }
  }
  SECTION("periodic in x and y: a corner reaches the origin"){
    const std::vector<bool> per = {true, true, false};
    const std::vector<std::uint32_t> master =
      mesh_tables::match_periodic_vertices(coords, nverts, 2, per, x_min, x_max, 1e-12);
    REQUIRE(master[at(n, n)] == at(0, 0));
    REQUIRE(master[at(n, 2)] == at(0, 2));
    REQUIRE(master[at(2, n)] == at(2, 0));
  }
  SECTION("no periodic direction is the identity"){
    const std::vector<bool> per = {false, false, false};
    const std::vector<std::uint32_t> master =
      mesh_tables::match_periodic_vertices(coords, nverts, 2, per, x_min, x_max, 1e-12);
    for (std::size_t v = 0; v < nverts; ++v) REQUIRE(master[v] == v);
  }
}
