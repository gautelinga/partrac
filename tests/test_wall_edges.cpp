// The near-wall edge table (near_wall::build_wall_edges) against a reference
// kept here: the straightforward construction, facets sorted by their vertex
// ids through a comparator and the side normals summed in a std::map. The
// table must match it bit for bit, since the near-wall P2 rule reads it in
// every wall cell and a changed summation order changes the velocities there.
// The fixtures are jittered lattices built here, every normal off the axes so
// that a sum depends on its order: a box with walls on every side (edges and
// corners), a slab one cell thick between two walls, periodic along it, where
// every cell touches a wall, a fully periodic box with no wall, and the same
// in 2D; each with its cells renumbered and their vertices rotated, which must
// leave every cell's entry unchanged.
#include <catch2/catch.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <map>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

#include "typedefs.hpp"
#include "Tet.hpp"
#include "Triangle.hpp"
#include "cell_walk.hpp"
#include "mesh_tables.hpp"
#include "near_wall.hpp"

namespace {

using near_wall::WallEdges;
using near_wall::WallRule;

// The reference construction
namespace ref {

Vector3d point_of(const std::vector<double>& coords, const std::size_t v, const Uint gdim){
  Vector3d p = Vector3d::Zero();
  for (Uint d = 0; d < gdim; ++d) p[d] = coords[v*gdim + d];
  return p;
}

template<int NV>
std::array<std::uint32_t, NV-1> facet_vertices(const std::uint32_t* row, const int k){
  std::array<std::uint32_t, NV-1> w{};
  int m = 0;
  for (int j = 0; j < NV; ++j)
    if (j != k) w[m++] = row[j];
  std::sort(w.begin(), w.end());
  return w;
}

template<int NV>
Vector3d wall_normal(const std::uint32_t* row, const int k, const std::vector<double>& coords,
                     const Uint gdim){
  const std::array<std::uint32_t, NV-1> w = facet_vertices<NV>(row, k);
  const Vector3d P0 = point_of(coords, row[k], gdim);
  const Vector3d P1 = point_of(coords, w[0], gdim);
  const Vector3d P2 = point_of(coords, w[1], gdim);
  double n[3];
  if constexpr (NV == 4){
    const Vector3d P3 = point_of(coords, w[2], gdim);
    const double V0[3] = {P0[0]-P1[0], P0[1]-P1[1], P0[2]-P1[2]};
    const double V1[3] = {P2[0]-P1[0], P2[1]-P1[1], P2[2]-P1[2]};
    const double V2[3] = {P3[0]-P1[0], P3[1]-P1[1], P3[2]-P1[2]};
    n[0] = V1[1]*V2[2] - V1[2]*V2[1];
    n[1] = V1[2]*V2[0] - V1[0]*V2[2];
    n[2] = V1[0]*V2[1] - V1[1]*V2[0];
    const double len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    for (int i = 0; i < 3; ++i) n[i] /= len;
    if (n[0]*V0[0] + n[1]*V0[1] + n[2]*V0[2] > 0.)
      for (int i = 0; i < 3; ++i) n[i] *= -1.;
  }
  else {
    double t[3] = {P2[0]-P1[0], P2[1]-P1[1], P2[2]-P1[2]};
    const double tlen = std::sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
    for (int i = 0; i < 3; ++i) t[i] /= tlen;
    for (int i = 0; i < 3; ++i) n[i] = P2[i] - P0[i];
    const double a = n[0]*t[0] + n[1]*t[1] + n[2]*t[2];
    for (int i = 0; i < 3; ++i) n[i] -= a*t[i];
    const double len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    for (int i = 0; i < 3; ++i) n[i] /= len;
  }
  return Vector3d(n[0], n[1], n[2]);
}

template<typename Cell>
void build_wall_edges(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                      const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                      const std::vector<std::int32_t>& facet_neigh,
                      const std::vector<std::uint32_t>& vclass,
                      std::vector<std::int32_t>& wall_index,
                      std::vector<WallEdges<Cell>>& wall_cells){
  constexpr int nv = Cell::n_verts;
  constexpr auto edge_ends = WallRule<Cell>::edge_ends;
  const auto klass = [&](const std::uint32_t v){ return std::size_t(vclass.empty() ? v : vclass[v]); };

  std::vector<std::size_t> slots;
  for (std::size_t s = 0; s < ncells*nv; ++s)
    if (facet_neigh[s] == facet_wall) slots.push_back(s);
  std::sort(slots.begin(), slots.end(), [&](const std::size_t a, const std::size_t b){
    return facet_vertices<nv>(topo.data() + (a/nv)*nv, int(a % nv))
         < facet_vertices<nv>(topo.data() + (b/nv)*nv, int(b % nv));
  });

  std::map<std::pair<std::size_t, std::size_t>, Vector3d> side_normal;
  std::vector<Vector3d> n_sum(nverts, Vector3d::Zero());
  std::vector<Vector3d> n_first(nverts, Vector3d::Zero());
  std::vector<std::uint32_t> n_count(nverts, 0);
  std::vector<bool> corner(nverts, false);
  for (const std::size_t s : slots){
    const std::uint32_t* row = topo.data() + (s/nv)*nv;
    const int k = int(s % nv);
    const Vector3d n = wall_normal<nv>(row, k, coords, gdim);
    for (int j = 0; j < nv; ++j){
      if (j == k) continue;
      const std::size_t iv = klass(row[j]);
      if (n_count[iv] == 0) n_first[iv] = n;
      else if ((nv == 3 && n_count[iv] > 1) || n_first[iv].dot(n) < 0.5) corner[iv] = true;
      n_sum[iv] += n;
      ++n_count[iv];
      side_normal.emplace(std::make_pair(klass(row[k]), klass(row[j])),
                          Vector3d::Zero().eval()).first->second += n;
    }
  }

  wall_index.assign(ncells, -1);
  wall_cells.clear();
  for (std::size_t i = 0; i < ncells; ++i){
    const std::uint32_t* row = topo.data() + i*nv;
    WallEdges<Cell> w;
    w.wall = 0;
    for (int k = 0; k < nv; ++k)
      if (n_count[klass(row[k])] > 0) w.wall |= 1 << k;
    if (w.wall == 0) continue;
    for (int e = 0; e < nv*(nv-1)/2; ++e){
      for (int o = 0; o < 2; ++o){
        typename WallRule<Cell>::End& we = w.ends[2*e + o];
        const std::uint32_t iw = row[edge_ends[std::size_t(e)][std::size_t(o)]];
        const std::uint32_t iv = row[edge_ends[std::size_t(e)][std::size_t(1-o)]];
        const std::size_t cw = klass(iw);
        if constexpr (nv == 3) we = {0.5, 0., 0., 0.5};
        else                   we = {0., 0., 0., 0., 0., 0.};
        if (n_count[cw] == 0 || corner[cw]) continue;
        const Vector3d edge = point_of(coords, iv, gdim) - point_of(coords, iw, gdim);
        const auto side = side_normal.find({klass(iv), cw});
        const Vector3d n = (side != side_normal.end() ? side->second : n_sum[cw]).normalized();
        const double delta = edge.dot(n);
        Vector3d q = -0.25*n;
        if constexpr (nv == 3){
          if (std::abs(delta) > 0.1*edge.norm())
            q += (edge - delta*n)/(2*delta);
          we = {0.5 + q[0]*n[0], q[0]*n[1], q[1]*n[0], 0.5 + q[1]*n[1]};
        }
        else {
          if (std::abs(delta) > 0.1*edge.norm())
            q += (edge - delta*n)/(4*delta);
          we = {q[0], q[1], q[2], n[0], n[1], n[2]};
        }
      }
    }
    wall_index[i] = std::int32_t(wall_cells.size());
    wall_cells.push_back(w);
  }
}

}  // namespace ref

// A lattice as the loaders hand it over: topology, coordinates, facet table, vertex classes
template<int NV>
struct Mesh {
  std::vector<std::uint32_t> topo;
  std::vector<double> coords;
  std::vector<std::int32_t> across;
  std::vector<std::uint32_t> vclass;
  std::size_t ncells() const { return topo.size()/NV; }
  std::size_t nverts() const { return coords.size()/(NV-1); }
};

// Smooth displacement of the unit box, periodic along every axis: zero across
// a periodic axis's faces, so the images still match, and off the plane of a
// wall, so no wall normal is axis-aligned
double jitter(const double* x, const int c, const int dim, const std::vector<bool>& periodic,
              const double a){
  const double pi = std::acos(-1.);
  double f = periodic[std::size_t(c)] ? std::sin(2*pi*x[c]) : std::cos(3*x[c] + 1.);
  for (int d = 1; d < dim; ++d)
    f *= 1. + 0.5*std::sin(2*pi*x[(c + d) % dim] + d);
  return a*f;
}

// The facet table and the vertex classes, as a loader builds them
template<int NV>
void finish(Mesh<NV>& m, const std::vector<bool>& periodic){
  const Uint dim = NV - 1;
  const Vector3d lo = Vector3d::Zero(), hi(1., 1., dim == 3 ? 1. : 0.);
  mesh_tables::build_facet_neighbours<NV>(m.topo, m.ncells(), m.coords, dim, periodic, lo, hi,
                                          1e-10, m.across);
  m.vclass = mesh_tables::match_periodic_vertices(m.coords, m.nverts(), dim, periodic, lo, hi, 1e-10);
}

// n[0] x n[1] x n[2] cubes, six Kuhn tets each, jittered
Mesh<4> tet_lattice(const std::array<std::size_t, 3> n, const std::vector<bool>& periodic){
  Mesh<4> m;
  const std::size_t mx = n[0] + 1, my = n[1] + 1, mz = n[2] + 1;
  auto vid = [&](const std::size_t i, const std::size_t j, const std::size_t k){
    return std::uint32_t((i*my + j)*mz + k);
  };
  m.coords.resize(3*mx*my*mz);
  for (std::size_t i = 0; i < mx; ++i)
    for (std::size_t j = 0; j < my; ++j)
      for (std::size_t k = 0; k < mz; ++k){
        const double x[3] = {double(i)/double(n[0]), double(j)/double(n[1]), double(k)/double(n[2])};
        double* c = m.coords.data() + 3*vid(i, j, k);
        for (int d = 0; d < 3; ++d)
          c[d] = x[d] + jitter(x, d, 3, periodic, 0.08/double(n[std::size_t(d)]));
      }
  for (std::size_t i = 0; i < n[0]; ++i)
    for (std::size_t j = 0; j < n[1]; ++j)
      for (std::size_t k = 0; k < n[2]; ++k){
        std::array<int, 3> p = {0, 1, 2};
        do {
          std::array<std::size_t, 3> o = {0, 0, 0};
          for (int s = 0; s < 4; ++s){
            m.topo.push_back(vid(i + o[0], j + o[1], k + o[2]));
            if (s < 3) o[std::size_t(p[std::size_t(s)])] = 1;
          }
        } while (std::next_permutation(p.begin(), p.end()));
      }
  finish(m, periodic);
  return m;
}

// n[0] x n[1] squares, two triangles each, jittered
Mesh<3> triangle_lattice(const std::array<std::size_t, 2> n, const std::vector<bool>& periodic){
  Mesh<3> m;
  const std::size_t mx = n[0] + 1, my = n[1] + 1;
  auto vid = [&](const std::size_t i, const std::size_t j){ return std::uint32_t(i*my + j); };
  m.coords.resize(2*mx*my);
  for (std::size_t i = 0; i < mx; ++i)
    for (std::size_t j = 0; j < my; ++j){
      const double x[2] = {double(i)/double(n[0]), double(j)/double(n[1])};
      double* c = m.coords.data() + 2*vid(i, j);
      for (int d = 0; d < 2; ++d)
        c[d] = x[d] + jitter(x, d, 2, periodic, 0.08/double(n[std::size_t(d)]));
    }
  for (std::size_t i = 0; i < n[0]; ++i)
    for (std::size_t j = 0; j < n[1]; ++j){
      const std::uint32_t v00 = vid(i, j), v10 = vid(i+1, j), v11 = vid(i+1, j+1), v01 = vid(i, j+1);
      m.topo.insert(m.topo.end(), {v00, v10, v11});
      m.topo.insert(m.topo.end(), {v00, v11, v01});
    }
  finish(m, periodic);
  return m;
}

// The cells in a random order, each one's vertices rotated; perm[new] = old
template<int NV>
Mesh<NV> renumbered(const Mesh<NV>& m, const std::vector<bool>& periodic,
                    std::vector<std::size_t>& perm, std::vector<int>& rot){
  std::mt19937 rng(7);
  perm.resize(m.ncells());
  std::iota(perm.begin(), perm.end(), std::size_t(0));
  std::shuffle(perm.begin(), perm.end(), rng);
  rot.resize(m.ncells());
  Mesh<NV> r;
  r.coords = m.coords;
  for (std::size_t c = 0; c < m.ncells(); ++c){
    rot[c] = int(rng() % NV);
    for (int k = 0; k < NV; ++k)
      r.topo.push_back(m.topo[perm[c]*NV + std::size_t((k + rot[c]) % NV)]);
  }
  finish(r, periodic);
  return r;
}

// Bitwise equality of two entries: every end, and the wall bits
template<typename Cell>
bool same_bits(const WallEdges<Cell>& a, const WallEdges<Cell>& b){
  return std::memcmp(a.ends.data(), b.ends.data(), sizeof(a.ends)) == 0 && a.wall == b.wall;
}

// Both constructions on m, compared bit for bit; returns the number of wall cells
template<typename Cell>
std::size_t check_against_reference(const Mesh<Cell::n_verts>& m){
  const Uint dim = Cell::n_verts - 1;
  std::vector<std::int32_t> index, ref_index;
  std::vector<WallEdges<Cell>> cells, ref_cells;
  near_wall::build_wall_edges<Cell>(m.topo, m.coords, m.ncells(), m.nverts(), dim, m.across,
                                    m.vclass, index, cells, false);
  ref::build_wall_edges<Cell>(m.topo, m.coords, m.ncells(), m.nverts(), dim, m.across,
                              m.vclass, ref_index, ref_cells);
  REQUIRE(index == ref_index);
  REQUIRE(cells.size() == ref_cells.size());
  std::size_t differ = 0;
  for (std::size_t i = 0; i < cells.size(); ++i)
    differ += !same_bits(cells[i], ref_cells[i]);
  REQUIRE(differ == 0);
  return cells.size();
}

// A renumbered mesh gives each cell the entry it had, its ends rotated with its vertices
template<typename Cell>
void check_renumbering(const Mesh<Cell::n_verts>& m, const std::vector<bool>& periodic){
  constexpr int nv = Cell::n_verts;
  const Uint dim = nv - 1;
  std::vector<std::size_t> perm;
  std::vector<int> rot;
  const Mesh<nv> r = renumbered(m, periodic, perm, rot);
  check_against_reference<Cell>(r);
  std::vector<std::int32_t> index, r_index;
  std::vector<WallEdges<Cell>> cells, r_cells;
  near_wall::build_wall_edges<Cell>(m.topo, m.coords, m.ncells(), m.nverts(), dim, m.across,
                                    m.vclass, index, cells, false);
  near_wall::build_wall_edges<Cell>(r.topo, r.coords, r.ncells(), r.nverts(), dim, r.across,
                                    r.vclass, r_index, r_cells, false);
  constexpr auto edge_ends = WallRule<Cell>::edge_ends;
  std::size_t differ = 0;
  for (std::size_t c = 0; c < r.ncells(); ++c){
    const std::int32_t a = r_index[c], b = index[perm[c]];
    REQUIRE((a < 0) == (b < 0));
    if (a < 0) continue;
    const WallEdges<Cell>& x = r_cells[std::size_t(a)];
    const WallEdges<Cell>& y = cells[std::size_t(b)];
    // local vertex k of the renumbered cell is (k + rot) of the original
    for (int k = 0; k < nv; ++k)
      differ += (x.wall >> k & 1) != (y.wall >> ((k + rot[c]) % nv) & 1);
    for (std::size_t e = 0; e < edge_ends.size(); ++e)
      for (std::size_t o = 0; o < 2; ++o){
        const int u = (edge_ends[e][o] + rot[c]) % nv, v = (edge_ends[e][1-o] + rot[c]) % nv;
        std::size_t f = 0, fo = 0;
        for (std::size_t g = 0; g < edge_ends.size(); ++g){
          if (edge_ends[g][0] == u && edge_ends[g][1] == v){ f = g; fo = 0; }
          if (edge_ends[g][1] == u && edge_ends[g][0] == v){ f = g; fo = 1; }
        }
        differ += std::memcmp(&x.ends[2*e + o], &y.ends[2*f + fo], sizeof(x.ends[0])) != 0;
      }
  }
  REQUIRE(differ == 0);
}

}  // namespace

TEST_CASE("The wall edge table equals the reference construction bit for bit, in 3D", "[interpol][wall]") {
  SECTION("a box walled on every side: edges and corners") {
    const std::vector<bool> periodic = {false, false, false};
    const Mesh<4> m = tet_lattice({5, 4, 3}, periodic);
    REQUIRE(check_against_reference<Tet>(m) > 0);
    check_renumbering<Tet>(m, periodic);
  }
  SECTION("a slab one cell thick between two walls, periodic along it: every cell a wall cell") {
    const std::vector<bool> periodic = {true, true, false};
    const Mesh<4> m = tet_lattice({6, 5, 1}, periodic);
    REQUIRE(check_against_reference<Tet>(m) == m.ncells());
    check_renumbering<Tet>(m, periodic);
  }
  SECTION("a channel periodic in x") {
    const std::vector<bool> periodic = {true, false, false};
    const Mesh<4> m = tet_lattice({4, 3, 3}, periodic);
    REQUIRE(check_against_reference<Tet>(m) > 0);
    check_renumbering<Tet>(m, periodic);
  }
  SECTION("a fully periodic box: no wall cell") {
    const std::vector<bool> periodic = {true, true, true};
    const Mesh<4> m = tet_lattice({3, 3, 3}, periodic);
    REQUIRE(check_against_reference<Tet>(m) == 0);
  }
}

TEST_CASE("The wall edge table equals the reference construction bit for bit, in 2D", "[interpol][wall]") {
  SECTION("a box walled on every side") {
    const std::vector<bool> periodic = {false, false, false};
    const Mesh<3> m = triangle_lattice({7, 5}, periodic);
    REQUIRE(check_against_reference<Triangle>(m) > 0);
    check_renumbering<Triangle>(m, periodic);
  }
  SECTION("a strip one cell thick between two walls, periodic along it: every cell a wall cell") {
    const std::vector<bool> periodic = {true, false, false};
    const Mesh<3> m = triangle_lattice({8, 1}, periodic);
    REQUIRE(check_against_reference<Triangle>(m) == m.ncells());
    check_renumbering<Triangle>(m, periodic);
  }
  SECTION("a fully periodic box: no wall cell") {
    const std::vector<bool> periodic = {true, true, false};
    const Mesh<3> m = triangle_lattice({4, 4}, periodic);
    REQUIRE(check_against_reference<Triangle>(m) == 0);
  }
}
