#include "near_wall.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <parallel/algorithm>

#include "Error.hpp"
#include "cell_walk.hpp"
#include "mesh_tables.hpp"
#include "wall_split.hpp"

namespace near_wall {

namespace {

// A vertex as a point, the third coordinate zero in 2D
Vector3d point_of(const std::vector<double>& coords, const std::size_t v, const Uint gdim){
  Vector3d p = Vector3d::Zero();
  for (Uint d = 0; d < gdim; ++d) p[d] = coords[v*gdim + d];
  return p;
}

// Whether slot s's facet is a wall: walls' bit, by cell, else a facet_wall entry
inline bool is_wall(const std::vector<std::int32_t>& facet_neigh, const std::vector<std::uint8_t>& walls,
                    const std::size_t s, const int nv){
  return walls.empty() ? facet_neigh[s] == facet_wall : (walls[s/std::size_t(nv)] >> (s % std::size_t(nv)) & 1) != 0;
}

// A vertex's class: its periodic master, itself without classes
inline std::uint32_t class_of(const std::vector<std::uint32_t>& vclass, const std::uint32_t v){
  return vclass.empty() ? v : vclass[v];
}

// The facet's vertices, ascending, as the entity numbering has them
template<int NV>
std::array<std::uint32_t, NV-1> facet_vertices(const std::uint32_t* row, const int k){
  std::array<std::uint32_t, NV-1> w{};
  int m = 0;
  for (int j = 0; j < NV; ++j)
    if (j != k) w[m++] = row[j];
  std::sort(w.begin(), w.end());
  return w;
}

// Outward unit normal of the facet w facing vertex k, in the arithmetic the
// element libraries use: the cross product of two facet edges in 3D, the part
// of an edge across the facet in 2D, signed by the vertex it faces
template<int NV>
Vector3d wall_normal(const std::uint32_t* row, const int k, const std::array<std::uint32_t, NV-1>& w,
                     const std::vector<double>& coords, const Uint gdim){
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

// Wall normals at the vertices and per (fluid end, wall end); if asked
// also every facet's classes, apex class and normal
template<int nv>
struct Geometry {
  struct Side { std::uint32_t wall; Vector3d n; };
  struct FacetSum { std::array<std::uint32_t, nv> key; Vector3d n; };   // classes, apex class last
  std::vector<FacetSum> facet;   // in the facet order
  std::vector<Vector3d> n_sum;
  std::vector<std::uint8_t> on_wall;
  std::vector<std::size_t> first;
  std::vector<Side> side;
};

template<int nv>
Geometry<nv> wall_geometry(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                           const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                           const std::vector<std::int32_t>& facet_neigh, const std::vector<std::uint8_t>& walls,
                           const std::vector<std::uint32_t>& vclass, const bool facet_sums){
  Geometry<nv> g;

  // The wall facets, ordered by their vertex ids: a facet's place in that order
  // is its place in the entity numbering, so a vertex sums its normals in one
  // sequence whatever the cells are numbered by
  struct Facet { std::array<std::uint32_t, nv-1> key; std::size_t slot; };
  std::size_t nw = 0;
  for (std::size_t s = 0; s < ncells*nv; ++s) nw += is_wall(facet_neigh, walls, s, nv);
  if (nw > std::numeric_limits<std::uint32_t>::max())
    partrac::fail(nw, " wall facets: more than the wall edge table can number");
  std::vector<Facet> facets(nw);
  for (std::size_t s = 0, f = 0; s < ncells*nv; ++s)
    if (is_wall(facet_neigh, walls, s, nv)) facets[f++].slot = s;
#pragma omp parallel for schedule(static)
  for (std::size_t f = 0; f < nw; ++f)
    facets[f].key = facet_vertices<nv>(topo.data() + (facets[f].slot/nv)*nv, int(facets[f].slot % nv));
  __gnu_parallel::sort(facets.begin(), facets.end(), [](const Facet& a, const Facet& b){
    return a.key < b.key || (a.key == b.key && a.slot < b.slot);
  });
  std::vector<Vector3d> normal(nw);
#pragma omp parallel for schedule(static)
  for (std::size_t f = 0; f < nw; ++f){
    const std::size_t s = facets[f].slot;
    normal[f] = wall_normal<nv>(topo.data() + (s/nv)*nv, int(s % nv), facets[f].key, coords, gdim);
  }
  if (facet_sums){
    g.facet.resize(nw);
    for (std::size_t f = 0; f < nw; ++f){
      const std::uint32_t* row = topo.data() + (facets[f].slot/nv)*nv;
      for (int j = 0; j < nv - 1; ++j) g.facet[f].key[std::size_t(j)] = class_of(vclass, facets[f].key[std::size_t(j)]);
      g.facet[f].key[nv-1] = class_of(vclass, row[facets[f].slot % nv]);
      g.facet[f].n = normal[f];
    }
  }

  // Wall normals at the vertices, summed over the facets that meet there in that
  // one order; the side entries (wall end, facet) bucketed by the fluid end
  std::vector<Vector3d>& n_sum = g.n_sum;
  std::vector<std::uint8_t>& on_wall = g.on_wall;
  std::vector<std::size_t>& first = g.first;
  n_sum.assign(nverts, Vector3d::Zero());
  on_wall.assign(nverts, 0);
  first.assign(nverts + 1, 0);
  for (std::size_t f = 0; f < nw; ++f){
    const std::uint32_t* row = topo.data() + (facets[f].slot/nv)*nv;
    const int k = int(facets[f].slot % nv);
    const Vector3d& n = normal[f];
    for (int j = 0; j < nv; ++j){
      if (j == k) continue;
      const std::size_t iv = class_of(vclass, row[j]);
      n_sum[iv] += n;
      on_wall[iv] = 1;
    }
    first[class_of(vclass, row[k]) + 1] += nv - 1;
  }
  for (std::size_t v = 0; v < nverts; ++v) first[v + 1] += first[v];
  struct Entry { std::uint32_t wall, facet; };
  std::vector<Entry> entry(first[nverts]);
  {
    std::vector<std::size_t> pos(first.begin(), first.end() - 1);
    for (std::size_t f = 0; f < nw; ++f){
      const std::uint32_t* row = topo.data() + (facets[f].slot/nv)*nv;
      const int k = int(facets[f].slot % nv);
      // the facet's own cell holds the apex it faces
      std::size_t& p = pos[class_of(vclass, row[k])];
      for (int j = 0; j < nv; ++j)
        if (j != k) entry[p++] = {class_of(vclass, row[j]), std::uint32_t(f)};
    }
  }
  std::vector<Facet>().swap(facets);

  // The side normals per fluid end and wall end: each pair's facets summed in
  // the facet order, normalized; a fluid end's pairs ascending by wall end
  std::vector<typename Geometry<nv>::Side>& side = g.side;
  side.reserve(entry.size());
  for (std::size_t v = 0, b = 0; v < nverts; ++v){
    const std::size_t e = first[v + 1];
    for (std::size_t i = b + 1; i < e; ++i)   // stable, the buckets are small
      for (std::size_t j = i; j > b && entry[j-1].wall > entry[j].wall; --j)
        std::swap(entry[j-1], entry[j]);
    first[v] = side.size();
    for (std::size_t i = b; i < e; ){
      Vector3d sum = Vector3d::Zero();
      const std::uint32_t w = entry[i].wall;
      for (; i < e && entry[i].wall == w; ++i) sum += normal[entry[i].facet];
      side.push_back({w, sum.normalized()});
    }
    b = e;
  }
  first[nverts] = side.size();
  return g;
}

// The side normal of the fluid end cv and wall end cw, as summed and normalized,
// else the sum of the normals at cw
struct EndNormal { const Vector3d& n; bool side; };
template<int nv>
EndNormal edge_end_normal(const Geometry<nv>& g, const std::uint32_t cw, const std::uint32_t cv){
  std::size_t p = g.first[cv];
  while (p < g.first[cv + 1] && g.side[p].wall < cw) ++p;
  if (p < g.first[cv + 1] && g.side[p].wall == cw) return {g.side[p].n, true};
  return {g.n_sum[cw], false};
}

// The cells with a wall vertex, numbered in cell order; verbose prints their count
template<int nv>
std::size_t number_wall_cells(const std::vector<std::uint32_t>& topo, const std::size_t ncells,
                              const std::vector<std::uint8_t>& on_wall, const std::vector<std::uint32_t>& vclass,
                              std::vector<std::int32_t>& wall_index, const bool verbose){
  wall_index.assign(ncells, -1);
  std::size_t ncw = 0;
  for (std::size_t i = 0; i < ncells; ++i){
    const std::uint32_t* row = topo.data() + i*std::size_t(nv);
    for (int k = 0; k < nv; ++k)
      if (on_wall[class_of(vclass, row[k])]){ wall_index[i] = std::int32_t(ncw++); break; }
  }
  if (verbose)
    std::cout << "Wall cells: " << ncw << std::endl;
  return ncw;
}

// fill(i, w, row, wall) for wall cell w, cell i, in parallel; wall: its wall vertices by bit
template<int nv, typename Fill>
void each_wall_cell(const std::vector<std::uint32_t>& topo, const std::vector<std::uint8_t>& on_wall,
                    const std::vector<std::uint32_t>& vclass, const std::vector<std::int32_t>& wall_index,
                    Fill&& fill){
#pragma omp parallel for schedule(dynamic, 4096)
  for (std::size_t i = 0; i < wall_index.size(); ++i){
    if (wall_index[i] < 0) continue;
    const std::uint32_t* row = topo.data() + i*nv;
    std::uint8_t wall = 0;
    for (int k = 0; k < nv; ++k)
      if (on_wall[class_of(vclass, row[k])]) wall |= std::uint8_t(1 << k);
    fill(i, std::size_t(wall_index[i]), row, wall);
  }
}

// Normals summed over the facets at a wall edge, by classes (lo, hi, apex);
// apex any_apex: over every facet there
using SumKey = std::array<std::uint32_t, 3>;
constexpr std::uint32_t any_apex = std::numeric_limits<std::uint32_t>::max();
struct NormalSum { SumKey key; Vector3d n; };

// Sorted by key, the sums of a key merged in their order
void merge_sums(std::vector<NormalSum>& v){
  std::stable_sort(v.begin(), v.end(), [](const NormalSum& a, const NormalSum& b){ return a.key < b.key; });
  std::size_t m = 0;
  for (std::size_t i = 0; i < v.size(); ){
    NormalSum s{v[i].key, Vector3d::Zero()};
    for (; i < v.size() && v[i].key == s.key; ++i) s.n += v[i].n;
    v[m++] = s;
  }
  v.resize(m);
}

const NormalSum* find_sum(const std::vector<NormalSum>& v, const SumKey& k){
  const auto it = std::lower_bound(v.begin(), v.end(), k, [](const NormalSum& a, const SumKey& b){ return a.key < b; });
  return it != v.end() && it->key == k ? &*it : nullptr;
}

// Unit, D values to out, or false where too short
template<int D>
bool unit_normal(const Vector3d& v, double* out){
  const double len = v.norm();
  if (!(len > 1e-12)) return false;
  for (int d = 0; d < D; ++d) out[d] = v[d]/len;
  return true;
}

// What a wall cell's wall terms are made from
template<int nv>
struct TermSources {
  const Geometry<nv>& g;
  std::vector<NormalSum> edge_sum, face_sum;
  const std::vector<std::int32_t>& facet_neigh;
  const std::vector<std::uint32_t>& vclass;
  const std::vector<double>& coords;
  const Uint gdim;
};

// Cell i's wall terms from its entities at rest, by bit, their normals to out
// in that order: edge ends from a wall vertex, but along a wall edge; side
// faces of a wall edge, but wall facets; in 3D the bubble's facets
template<typename Cell>
std::uint32_t cell_terms(const TermSources<Cell::n_verts>& s, const std::uint32_t* row, const unsigned rest,
                         double* out){
  constexpr int nv = Cell::n_verts, D = nv - 1;
  using S = wall_split::Shape<nv>;
  static constexpr auto edges = split_eval::lex_edges<nv>();
  static constexpr auto apex = wall_split::edge_apices<nv>();
  const auto cls = [&](const int k){ return class_of(s.vclass, row[k]); };
  const auto vertex = [&](const int v){ return (rest >> (nv == 4 ? S::n_bits + v : nv + v) & 1) != 0; };
  // 2D: the wall edge is the facet off it
  const auto edge = [&](const int e){
    return (rest >> (nv == 4 ? nv + e : 3 - edges[std::size_t(e)][0] - edges[std::size_t(e)][1]) & 1) != 0;
  };
  std::uint32_t terms = 0;
  // The side normal, else the mean at the wall vertex
  for (int e = 0; e < S::n_edges; ++e)
    for (int o = 0; o < 2; ++o){
      const int w = edges[std::size_t(e)][std::size_t(o)], v = edges[std::size_t(e)][std::size_t(1 - o)];
      if (!vertex(w) || edge(e)) continue;
      if (unit_normal<D>(edge_end_normal(s.g, cls(w), cls(v)).n, out)){
        terms |= 1u << (2*e + o);
        out += D;
      }
    }
  // The facets there with this apex, else all of them
  for (int e = 0; e < S::n_edges; ++e)
    for (int j = 0; j < nv - 2; ++j){
      const int a = edges[std::size_t(e)][0], b = edges[std::size_t(e)][1], v = apex[std::size_t(e)][std::size_t(j)];
      if (!edge(e) || (nv == 4 && (rest >> (6 - a - b - v) & 1))) continue;
      const std::uint32_t ca = cls(a), cb = cls(b);
      const SumKey k{std::min(ca, cb), std::max(ca, cb), any_apex};
      const NormalSum* es = find_sum(s.edge_sum, k);
      if (!es) continue;
      const NormalSum* fs = find_sum(s.face_sum, {k[0], k[1], cls(v)});
      if (unit_normal<D>(fs ? fs->n : es->n, out)){
        terms |= 1u << (S::triple0 + e*(nv - 2) + j);
        out += D;
      }
    }
  // Each facet's normal, then its column of the facets' pseudo-inverse
  if constexpr (nv == 4){
    int ks[nv], m = 0;
    for (int k = 0; k < nv; ++k)
      if (rest >> k & 1) ks[m++] = k;
    if (m == 0) return terms;
    Eigen::Matrix<double, Eigen::Dynamic, D, Eigen::RowMajor, nv, D> N(m, D);
    for (int r = 0; r < m; ++r){
      const Vector3d n = wall_normal<nv>(row, ks[r], facet_vertices<nv>(row, ks[r]), s.coords, s.gdim);
      for (int d = 0; d < D; ++d) N(r, d) = n[d];
    }
    const Eigen::Matrix<double, D, Eigen::Dynamic, 0, D, nv> P = N.completeOrthogonalDecomposition().pseudoInverse();
    for (int r = 0; r < m; ++r){
      for (int d = 0; d < D; ++d){
        out[d] = N(r, d);
        out[D + d] = P(d, r);
      }
      out += 2*D;
      terms |= 1u << (S::bubble0 + ks[r]);
    }
  }
  return terms;
}

}  // namespace

template<typename Cell>
void build_wall_edges(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                      const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                      const std::vector<std::int32_t>& facet_neigh, const std::vector<std::uint8_t>& wall_facets,
                      const std::vector<std::uint32_t>& vclass,
                      std::vector<std::int32_t>& wall_index,
                      std::vector<WallEdges<Cell>>& wall_cells, const bool verbose){
  constexpr int nv = Cell::n_verts;
  constexpr auto edge_ends = WallRule<Cell>::edge_ends;
  const Geometry<nv> g = wall_geometry<nv>(topo, coords, ncells, nverts, gdim, facet_neigh, wall_facets, vclass, false);

  wall_cells.clear();
  wall_cells.resize(number_wall_cells<nv>(topo, ncells, g.on_wall, vclass, wall_index, verbose));
  each_wall_cell<nv>(topo, g.on_wall, vclass, wall_index,
                     [&](std::size_t, const std::size_t wc, const std::uint32_t* row, const std::uint8_t wall){
    WallEdges<Cell>& w = wall_cells[wc];
    w.wall = wall;
    for (int e = 0; e < nv*(nv-1)/2; ++e){
      for (int o = 0; o < 2; ++o){
        typename WallRule<Cell>::End& we = w.ends[2*e + o];
        const std::uint32_t iw = row[edge_ends[std::size_t(e)][std::size_t(o)]];
        const std::uint32_t iv = row[edge_ends[std::size_t(e)][std::size_t(1-o)]];
        const std::uint32_t cw = class_of(vclass, iw);
        if constexpr (nv == 3) we = {0.5, 0., 0., 0.5};   // linear
        else                   we = {0., 0., 0., 0., 0., 0.};
        if (!g.on_wall[cw]) continue;
        const Vector3d edge = point_of(coords, iv, gdim) - point_of(coords, iw, gdim);
        // over several wall facets: their mean normal
        const EndNormal en = edge_end_normal(g, cw, class_of(vclass, iv));
        const Vector3d n = en.side ? en.n : en.n.normalized();
        const double delta = edge.dot(n);
        // v_n ~ delta^2, divergence-free wall cell
        Vector3d q = -0.25*n;
        if constexpr (nv == 3){
          if (std::abs(delta) > 0.1*edge.norm())
            q += (edge - delta*n)/(2*delta);
          // M = I/2 + q n^T
          we = {0.5 + q[0]*n[0], q[0]*n[1], q[1]*n[0], 0.5 + q[1]*n[1]};
        }
        else {
          if (std::abs(delta) > 0.1*edge.norm())
            q += (edge - delta*n)/(4*delta);
          we = {q[0], q[1], q[2], n[0], n[1], n[2]};
        }
      }
    }
  });
}

template<typename Cell>
void rebuild_wall_edges_at_rest(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                                const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                                const std::vector<std::int32_t>& facet_neigh, const std::vector<std::uint32_t>& vclass,
                                std::vector<std::int32_t>& wall_index, std::vector<WallEdges<Cell>>& wall_cells)
{
  constexpr int nv = Cell::n_verts;
  constexpr unsigned all = (1u << nv) - 1;
  std::vector<std::uint8_t> walls(ncells, 0);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < ncells; ++i){
    if (wall_index[i] < 0) continue;
    const unsigned rest = wall_cells[std::size_t(wall_index[i])].wall;
    for (int k = 0; k < nv; ++k)
      if (facet_neigh[i*nv + std::size_t(k)] == facet_wall && (rest | 1u << k) == all)
        walls[i] = std::uint8_t(walls[i] | 1u << k);
  }
  build_wall_edges<Cell>(topo, coords, ncells, nverts, gdim, facet_neigh, walls, vclass, wall_index, wall_cells,
                         false);
}

template<typename Cell>
void build_tables(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                  const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                  const std::vector<bool>& periodic, const Vector3d& x_min, const Vector3d& x_max,
                  const double periodic_tol, const std::vector<std::int32_t>& facet_neigh,
                  const std::vector<std::uint8_t>& wall_facets, const bool p2, std::vector<std::uint32_t>& vclass,
                  std::vector<std::int32_t>& wall_index, std::vector<WallEdges<Cell>>& wall_cells,
                  std::vector<wall_split::WallCell>& walls, const bool verbose)
{
  if (p2){
    if (vclass.empty() && std::find(periodic.begin(), periodic.end(), true) != periodic.end())
      vclass = mesh_tables::match_periodic_vertices(coords, nverts, gdim, periodic, x_min, x_max,
                                                    periodic_tol);
    build_wall_patterns<Cell>(topo, ncells, nverts, facet_neigh, wall_facets, vclass, wall_index, walls, verbose);
    return;
  }
  if (vclass.empty())
    vclass = mesh_tables::match_periodic_vertices(coords, nverts, gdim, periodic, x_min, x_max,
                                                  periodic_tol);
  build_wall_edges<Cell>(topo, coords, ncells, nverts, gdim, facet_neigh, wall_facets, vclass, wall_index,
                         wall_cells, verbose);
}

template<typename Cell>
void build_wall_patterns(const std::vector<std::uint32_t>& topo, const std::size_t ncells,
                         const std::size_t nverts, const std::vector<std::int32_t>& facet_neigh,
                         const std::vector<std::uint8_t>& wall_facets, const std::vector<std::uint32_t>& vclass,
                         std::vector<std::int32_t>& wall_index, std::vector<wall_split::WallCell>& walls,
                         const bool verbose)
{
  constexpr int nv = Cell::n_verts;
  constexpr auto ev = wall_split::entity_verts<nv>();
  const auto pair_key = [](const std::uint32_t a, const std::uint32_t b){
    return std::uint64_t(std::min(a, b)) << 32 | std::uint64_t(std::max(a, b));
  };
  // The wall vertices by class; in 3D the wall edges as class pairs, sorted
  std::vector<std::uint8_t> on_wall(nverts, 0);
  std::vector<std::uint64_t> edges;
  for (std::size_t s = 0; s < ncells*nv; ++s){
    if (!is_wall(facet_neigh, wall_facets, s, nv)) continue;
    const std::uint32_t* row = topo.data() + (s/nv)*nv;
    const int k = int(s % nv);
    std::uint32_t c[nv - 1];
    for (int j = 0, m = 0; j < nv; ++j)
      if (j != k) c[m++] = class_of(vclass, row[j]);
    for (int i = 0; i < nv - 1; ++i) on_wall[c[i]] = 1;
    if constexpr (nv == 4)
      for (int a = 0; a < 3; ++a)
        for (int b = a + 1; b < 3; ++b) edges.push_back(pair_key(c[a], c[b]));
  }
  __gnu_parallel::sort(edges.begin(), edges.end());
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

  // Each cell's pattern, held in its wall index until the cells are numbered
  wall_index.assign(ncells, 0);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < ncells; ++i){
    const std::uint32_t* row = topo.data() + i*nv;
    unsigned p = 0;
    for (int k = 0; k < nv; ++k)
      if (is_wall(facet_neigh, wall_facets, i*nv + std::size_t(k), nv)) p |= 1u << k;
    std::uint32_t c[nv];
    unsigned w = 0;
    for (int k = 0; k < nv; ++k){
      c[k] = class_of(vclass, row[k]);
      if (on_wall[c[k]]) w |= 1u << k;
    }
    for (std::size_t q = 0; q < ev.size(); ++q){
      const int a = ev[q][0], b = ev[q][1];
      if (!(w >> a & 1) || !(w >> b & 1)) continue;
      if (a == b || std::binary_search(edges.begin(), edges.end(), pair_key(c[a], c[b])))
        p |= 1u << (nv + int(q));
    }
    if constexpr (nv == 4) p |= w << wall_split::Shape<nv>::n_bits;
    wall_index[i] = std::int32_t(p);
  }
  walls.clear();
  for (std::size_t i = 0; i < ncells; ++i){
    const std::int32_t p = wall_index[i];
    if (!p){
      wall_index[i] = -1;
      continue;
    }
    if (walls.size() >= std::size_t(std::numeric_limits<std::int32_t>::max()))
      partrac::fail(ncells, " cells: more wall cells than the wall index can number");
    wall_index[i] = std::int32_t(walls.size());
    wall_split::WallCell c;
    c.rest = std::uint16_t(p);
    walls.push_back(c);
  }
  if (verbose)
    std::cout << "Wall cells: " << walls.size() << std::endl;
}

template<typename Cell>
void build_wall_terms(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                      const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                      const std::vector<std::int32_t>& facet_neigh, const std::vector<std::uint32_t>& vclass,
                      std::vector<std::int32_t>& wall_index, wall_split::Tables<Cell::n_verts>& t,
                      const bool verbose)
{
  constexpr int nv = Cell::n_verts, D = nv - 1;
  using S = wall_split::Shape<nv>;
  // Nothing at rest: read as without the rule; the others renumbered in cell order
  std::size_t live = 0;
  for (std::size_t i = 0; i < ncells; ++i){
    if (wall_index[i] < 0) continue;
    const wall_split::WallCell c = t.cell[std::size_t(wall_index[i])];
    if (!c.rest){
      wall_index[i] = -1;
      continue;
    }
    t.cell[live] = c;
    wall_index[i] = std::int32_t(live++);
  }
  t.cell.resize(live);
  // The wall facets at rest
  std::vector<std::uint8_t> walls(ncells, 0);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < ncells; ++i)
    if (wall_index[i] >= 0) walls[i] = std::uint8_t(t.cell[std::size_t(wall_index[i])].rest & ((1u << nv) - 1));
  const Geometry<nv> g = wall_geometry<nv>(topo, coords, ncells, nverts, gdim, facet_neigh, walls, vclass, true);

  // The facets' normals summed by wall edge, and by wall edge and apex, in the facet order
  TermSources<nv> src{g, {}, {}, facet_neigh, vclass, coords, gdim};
  src.edge_sum.reserve(g.facet.size()*std::size_t((nv - 1)*(nv - 2)/2));
  src.face_sum.reserve(src.edge_sum.capacity());
  for (const auto& f : g.facet)
    for (int i = 0; i < nv - 1; ++i)
      for (int j = i + 1; j < nv - 1; ++j){
        const std::uint32_t a = f.key[std::size_t(i)], b = f.key[std::size_t(j)];
        const std::uint32_t lo = std::min(a, b), hi = std::max(a, b);
        src.edge_sum.push_back({{lo, hi, any_apex}, f.n});
        src.face_sum.push_back({{lo, hi, f.key[nv-1]}, f.n});
      }
  merge_sums(src.edge_sum);
  merge_sums(src.face_sum);

  // The terms, then the places in the pool, then the normals
  std::vector<std::uint32_t> count(t.cell.size(), 0);
#pragma omp parallel for schedule(dynamic, 4096)
  for (std::size_t i = 0; i < ncells; ++i){
    if (wall_index[i] < 0) continue;
    wall_split::WallCell& c = t.cell[std::size_t(wall_index[i])];
    std::array<double, 2*S::n_terms*D> buf;
    c.terms = cell_terms<Cell>(src, topo.data() + i*nv, c.rest, buf.data());
    count[std::size_t(wall_index[i])] = std::uint32_t(__builtin_popcount(c.terms)
                                                      + __builtin_popcount(c.terms >> S::bubble0));
  }
  std::size_t nn = 0;
  for (std::size_t w = 0; w < t.cell.size(); ++w){
    t.cell[w].normals = std::uint32_t(nn);
    nn += count[w];
    if (nn > std::numeric_limits<std::uint32_t>::max())
      partrac::fail(t.cell.size(), " wall cells: more wall terms than their table can number");
  }
  std::vector<std::uint32_t>().swap(count);
  t.normals.assign(nn*std::size_t(D), 0.);
#pragma omp parallel for schedule(dynamic, 4096)
  for (std::size_t i = 0; i < ncells; ++i){
    if (wall_index[i] < 0) continue;
    const wall_split::WallCell& c = t.cell[std::size_t(wall_index[i])];
    cell_terms<Cell>(src, topo.data() + i*nv, c.rest, t.normals.data() + std::size_t(c.normals)*D);
  }
  t.build(verbose);
}

void warn_velocity_walls(const std::vector<std::int32_t>& wall_index)
{
  if (std::none_of(wall_index.begin(), wall_index.end(), [](const std::int32_t w){ return w >= 0; })) return;
  std::cerr << "Warning: wall_upgrade takes the walls from the velocity at rest, which misses no-slip nodes "
               "solved only to a tolerance; wall_facets names them" << std::endl;
}

double rest_tolerance(const std::vector<double>& u)
{
  double u_max = 0.;
  for (const double x : u) u_max = std::max(u_max, std::abs(x));
  return 1e-12*u_max;
}

double rest_tolerance_by_vertex(const std::vector<double>& u, const std::vector<std::uint32_t>& vclass,
                                const std::size_t nverts, const Uint dim)
{
  double u_max = 0.;
  for (std::size_t v = 0; v < nverts; ++v){
    if (vclass[v] != v) continue;
    for (Uint c = 0; c < dim; ++c) u_max = std::max(u_max, std::abs(u[v*dim + c]));
  }
  return 1e-12*u_max;
}

template<typename Cell>
bool check_moving_walls(const CellDofs& dofs, const std::vector<std::int32_t>& wall_index, const bool p2,
                        const std::vector<wall_split::WallCell>& walls, const std::vector<WallEdges<Cell>>& edges,
                        const std::vector<const std::vector<double>*>& us, const std::vector<std::uint32_t>& vclass,
                        const std::string& where, const std::vector<double>& amps)
{
  constexpr int D = Cell::n_verts - 1;
  static constexpr auto slots = wall_split::entity_slots<Cell>();
  const auto magnitude = [](const double* x){
    double s = 0.;
    for (int d = 0; d < D; ++d) s += x[d]*x[d];
    return std::sqrt(s);
  };
  // Each field's largest |u|
  std::vector<double> u_max(us.size(), 0.);
  for (std::size_t k = 0; k < us.size(); ++k){
    const std::vector<double>& u = *us[k];
    const std::size_t n = u.size()/D;
    double m = 0.;
#pragma omp parallel for schedule(static) reduction(max:m)
    for (std::size_t v = 0; v < n; ++v)
      if (vclass.empty() || (v < vclass.size() && vclass[v] == v)) m = std::max(m, magnitude(&u[v*D]));
    u_max[k] = m;
  }
  // The sums judged: each field alone, else one of all weighted by |amps|
  std::vector<std::vector<std::size_t>> sums;
  if (amps.empty()) for (std::size_t k = 0; k < us.size(); ++k) sums.push_back({k});
  else {
    sums.emplace_back();
    for (std::size_t k = 0; k < us.size(); ++k) sums.back().push_back(k);
  }
  const auto weight = [&](const std::size_t k){ return amps.empty() ? 1. : std::abs(amps[k]); };
  bool flow = false;
  double ratio = 0.;
  std::vector<std::uint8_t> moving;
  for (const std::vector<std::size_t>& sum : sums){
    // Its scale, and at a node the bound of its |u|
    double scale = 0.;
    std::size_t n = 0;
    for (const std::size_t k : sum){
      scale += weight(k)*u_max[k];
      n = std::max(n, us[k]->size()/D);
    }
    if (!(scale > 0.)) continue;
    flow = true;
    moving.resize(std::max(moving.size(), n), 0);
    const double tol = 1e-6*scale;
    double r = 0.;
#pragma omp parallel for schedule(dynamic, 4096) reduction(max:r)
    for (std::size_t i = 0; i < wall_index.size(); ++i){
      const std::int32_t w = wall_index[i];
      if (w < 0) continue;
      unsigned nodes = 0;
      if (p2)
        for (unsigned m = walls[std::size_t(w)].rest; m; m &= m - 1) nodes |= slots[std::size_t(__builtin_ctz(m))];
      else nodes = edges[std::size_t(w)].wall;
      const std::uint32_t* row = dofs[i];
      for (; nodes; nodes &= nodes - 1){
        const std::size_t q = row[__builtin_ctz(nodes)];
        double a = 0.;
        for (const std::size_t k : sum) a += weight(k)*magnitude(&(*us[k])[q*D]);
        if (!(a > tol)) continue;
#pragma omp atomic write
        moving[q] = 1;
        r = std::max(r, a/scale);
      }
    }
    ratio = std::max(ratio, r);
  }
  const std::size_t count = std::size_t(std::count(moving.begin(), moving.end(), std::uint8_t(1)));
  if (count)
    std::cerr << "Warning: the walls wall_facets names move at " << count << " velocity nodes in " << where
              << ", |u| up to " << ratio << " of the largest |u|; they are kept as no-slip walls" << std::endl;
  return flow;
}

void label_cell_type(std::vector<int>& cell_type, const std::vector<std::int32_t>& across,
                     const Uint nv, const std::vector<std::uint8_t>& walls){
  const std::size_t ncells = across.size() / nv;
  cell_type.assign(ncells, 0);
  for (std::size_t i = 0; i < ncells; ++i){
    if (!walls.empty()){
      if (walls[i]) cell_type[i] = 1;
      continue;
    }
    for (Uint k = 0; k < nv; ++k)
      if (across[i*nv + k] == facet_wall) cell_type[i] = 1;
  }
  for (std::size_t i = 0; i < ncells; ++i){
    if (cell_type[i] != 1) continue;
    for (Uint k = 0; k < nv; ++k){
      const std::int32_t a = across[i*nv + k];
      if (a == facet_wall) continue;
      const std::int32_t j = a >= 0 ? a : facet_periodic(a);
      if (cell_type[j] == 0) cell_type[j] = 2;
    }
  }
}

template void build_wall_edges<Triangle>(const std::vector<std::uint32_t>&, const std::vector<double>&,
                                         const std::size_t, const std::size_t, const Uint,
                                         const std::vector<std::int32_t>&, const std::vector<std::uint8_t>&,
                                         const std::vector<std::uint32_t>&,
                                         std::vector<std::int32_t>&, std::vector<WallEdges<Triangle>>&,
                                         const bool);
template void build_wall_edges<Tet>(const std::vector<std::uint32_t>&, const std::vector<double>&,
                                    const std::size_t, const std::size_t, const Uint,
                                    const std::vector<std::int32_t>&, const std::vector<std::uint8_t>&,
                                    const std::vector<std::uint32_t>&,
                                    std::vector<std::int32_t>&, std::vector<WallEdges<Tet>>&,
                                    const bool);

template void build_tables<Triangle>(const std::vector<std::uint32_t>&, const std::vector<double>&,
                                     const std::size_t, const std::size_t, const Uint, const std::vector<bool>&,
                                     const Vector3d&, const Vector3d&, const double,
                                     const std::vector<std::int32_t>&, const std::vector<std::uint8_t>&, const bool,
                                     std::vector<std::uint32_t>&,
                                     std::vector<std::int32_t>&, std::vector<WallEdges<Triangle>>&,
                                     std::vector<wall_split::WallCell>&, const bool);
template void build_tables<Tet>(const std::vector<std::uint32_t>&, const std::vector<double>&,
                                const std::size_t, const std::size_t, const Uint, const std::vector<bool>&,
                                const Vector3d&, const Vector3d&, const double,
                                const std::vector<std::int32_t>&, const std::vector<std::uint8_t>&, const bool,
                                std::vector<std::uint32_t>&,
                                std::vector<std::int32_t>&, std::vector<WallEdges<Tet>>&,
                                std::vector<wall_split::WallCell>&, const bool);
template void build_wall_patterns<Triangle>(const std::vector<std::uint32_t>&, const std::size_t, const std::size_t,
                                            const std::vector<std::int32_t>&, const std::vector<std::uint8_t>&,
                                            const std::vector<std::uint32_t>&,
                                            std::vector<std::int32_t>&, std::vector<wall_split::WallCell>&, const bool);
template void build_wall_patterns<Tet>(const std::vector<std::uint32_t>&, const std::size_t, const std::size_t,
                                       const std::vector<std::int32_t>&, const std::vector<std::uint8_t>&,
                                       const std::vector<std::uint32_t>&,
                                       std::vector<std::int32_t>&, std::vector<wall_split::WallCell>&, const bool);

template void build_wall_terms<Triangle>(const std::vector<std::uint32_t>&, const std::vector<double>&,
                                         const std::size_t, const std::size_t, const Uint,
                                         const std::vector<std::int32_t>&, const std::vector<std::uint32_t>&,
                                         std::vector<std::int32_t>&, wall_split::Tables<3>&, const bool);
template void rebuild_wall_edges_at_rest<Triangle>(const std::vector<std::uint32_t>&, const std::vector<double>&,
                                                   const std::size_t, const std::size_t, const Uint,
                                                   const std::vector<std::int32_t>&, const std::vector<std::uint32_t>&,
                                                   std::vector<std::int32_t>&, std::vector<WallEdges<Triangle>>&);
template void rebuild_wall_edges_at_rest<Tet>(const std::vector<std::uint32_t>&, const std::vector<double>&,
                                              const std::size_t, const std::size_t, const Uint,
                                              const std::vector<std::int32_t>&, const std::vector<std::uint32_t>&,
                                              std::vector<std::int32_t>&, std::vector<WallEdges<Tet>>&);
template bool check_moving_walls<Triangle>(const CellDofs&, const std::vector<std::int32_t>&, const bool,
                                           const std::vector<wall_split::WallCell>&,
                                           const std::vector<WallEdges<Triangle>>&,
                                           const std::vector<const std::vector<double>*>&,
                                           const std::vector<std::uint32_t>&, const std::string&,
                                           const std::vector<double>&);
template bool check_moving_walls<Tet>(const CellDofs&, const std::vector<std::int32_t>&, const bool,
                                      const std::vector<wall_split::WallCell>&, const std::vector<WallEdges<Tet>>&,
                                      const std::vector<const std::vector<double>*>&,
                                      const std::vector<std::uint32_t>&, const std::string&,
                                      const std::vector<double>&);
template void build_wall_terms<Tet>(const std::vector<std::uint32_t>&, const std::vector<double>&,
                                    const std::size_t, const std::size_t, const Uint,
                                    const std::vector<std::int32_t>&, const std::vector<std::uint32_t>&,
                                    std::vector<std::int32_t>&, wall_split::Tables<4>&, const bool);

}  // namespace near_wall
