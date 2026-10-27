#include "near_wall.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <parallel/algorithm>

#include "Error.hpp"
#include "cell_walk.hpp"
#include "mesh_tables.hpp"
#include "p12_eval.hpp"

namespace near_wall {

namespace {

// A vertex as a point, the third coordinate zero in 2D
Vector3d point_of(const std::vector<double>& coords, const std::size_t v, const Uint gdim){
  Vector3d p = Vector3d::Zero();
  for (Uint d = 0; d < gdim; ++d) p[d] = coords[v*gdim + d];
  return p;
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

// The wall geometry both rules read: the normals at the wall vertices, the
// corners, and the side normals per fluid end and wall end; with facet_sums
// also every facet's classes, apex class and normal
template<int nv>
struct Geometry {
  struct Side { std::uint32_t wall; Vector3d n; };
  struct FacetSum { std::array<std::uint32_t, nv> key; Vector3d n; };   // classes, apex class last
  std::vector<Vector3d> n_sum;
  std::vector<std::uint32_t> n_count;
  std::vector<bool> corner;
  std::vector<std::size_t> first;
  std::vector<Side> side;
  std::vector<FacetSum> facet;   // in the facet order
};

template<int nv>
Geometry<nv> wall_geometry(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                           const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                           const std::vector<std::int32_t>& facet_neigh,
                           const std::vector<std::uint32_t>& vclass, const bool facet_sums){
  Geometry<nv> g;

  // The wall facets, ordered by their vertex ids: a facet's place in that order
  // is its place in the entity numbering, so a vertex sums its normals in one
  // sequence whatever the cells are numbered by
  struct Facet { std::array<std::uint32_t, nv-1> key; std::size_t slot; };
  const std::size_t nw = std::size_t(std::count(facet_neigh.begin(), facet_neigh.end(), facet_wall));
  if (nw > std::numeric_limits<std::uint32_t>::max())
    partrac::fail(nw, " wall facets: more than the wall edge table can number");
  std::vector<Facet> facets(nw);
  for (std::size_t s = 0, f = 0; s < ncells*nv; ++s)
    if (facet_neigh[s] == facet_wall) facets[f++].slot = s;
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
  std::vector<std::uint32_t>& n_count = g.n_count;
  std::vector<bool>& corner = g.corner;
  std::vector<std::size_t>& first = g.first;
  n_sum.assign(nverts, Vector3d::Zero());
  std::vector<Vector3d> n_first(nverts, Vector3d::Zero());
  n_count.assign(nverts, 0);
  corner.assign(nverts, false);
  first.assign(nverts + 1, 0);
  for (std::size_t f = 0; f < nw; ++f){
    const std::uint32_t* row = topo.data() + (facets[f].slot/nv)*nv;
    const int k = int(facets[f].slot % nv);
    const Vector3d& n = normal[f];
    for (int j = 0; j < nv; ++j){
      if (j == k) continue;
      const std::size_t iv = class_of(vclass, row[j]);
      if (n_count[iv] == 0) n_first[iv] = n;
      else if ((nv == 3 && n_count[iv] > 1) || n_first[iv].dot(n) < 0.5) corner[iv] = true;   // over 60 degrees
      n_sum[iv] += n;
      ++n_count[iv];
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
                              const std::vector<std::uint32_t>& n_count, const std::vector<std::uint32_t>& vclass,
                              std::vector<std::int32_t>& wall_index, const bool verbose){
  wall_index.assign(ncells, -1);
  std::size_t ncw = 0;
  for (std::size_t i = 0; i < ncells; ++i){
    const std::uint32_t* row = topo.data() + i*std::size_t(nv);
    for (int k = 0; k < nv; ++k)
      if (n_count[class_of(vclass, row[k])] > 0){ wall_index[i] = std::int32_t(ncw++); break; }
  }
  if (verbose)
    std::cout << "Wall cells: " << ncw << std::endl;
  return ncw;
}

// fill(i, w, row, wall) for wall cell w, cell i, in parallel; wall: its wall vertices by bit
template<int nv, typename Fill>
void each_wall_cell(const std::vector<std::uint32_t>& topo, const std::vector<std::uint32_t>& n_count,
                    const std::vector<std::uint32_t>& vclass, const std::vector<std::int32_t>& wall_index,
                    Fill&& fill){
#pragma omp parallel for schedule(dynamic, 4096)
  for (std::size_t i = 0; i < wall_index.size(); ++i){
    if (wall_index[i] < 0) continue;
    const std::uint32_t* row = topo.data() + i*nv;
    std::uint8_t wall = 0;
    for (int k = 0; k < nv; ++k)
      if (n_count[class_of(vclass, row[k])] > 0) wall |= std::uint8_t(1 << k);
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

// Unit, or false where too short
template<std::size_t D>
bool unit_normal(const Vector3d& v, std::array<double, D>& out){
  const double len = v.norm();
  if (!(len > 1e-12)) return false;
  for (std::size_t d = 0; d < D; ++d) out[d] = v[int(d)]/len;
  return true;
}

// What a wall cell's terms are made from
template<int nv>
struct TermSources {
  const Geometry<nv>& g;
  std::vector<NormalSum> edge_sum, face_sum;
  const std::vector<std::int32_t>& facet_neigh;
  const std::vector<std::uint32_t>& vclass;
  const std::vector<double>& coords;
  const Uint gdim;
};

// A wall cell's terms and normals, unpacked
template<typename Cell>
struct CellTerms {
  using R = UpgradeRule<Cell>;
  using Normal = typename UpgradeTable<Cell>::Normal;
  static constexpr int f0 = 2*R::n_edges;
  std::array<Normal, 2*R::n_edges> end;
  std::array<Normal, R::n_edges*R::n_apex> face;
  std::array<Normal, R::n_bubble ? Cell::n_verts : 0> bn, bp;
  std::uint32_t terms = 0;
  std::uint16_t triples = 0;
  std::uint8_t facets = 0;
  // Its normals in the pool
  std::size_t n_normals() const {
    return std::size_t(__builtin_popcount(terms & ((1u << f0) - 1)) + __builtin_popcount(triples)
                       + 2*__builtin_popcount(facets));
  }
};

// Cell i's terms and normals
template<typename Cell>
void cell_terms(const TermSources<Cell::n_verts>& s, const std::size_t i, const std::uint32_t* row,
                CellTerms<Cell>& w){
  constexpr int nv = Cell::n_verts, D = nv - 1;
  using R = UpgradeRule<Cell>;
  constexpr auto edge_ends = WallRule<Cell>::edge_ends;
  constexpr int f0 = 2*R::n_edges, b0 = f0 + R::n_faces;
  const Geometry<nv>& g = s.g;
  const auto cls = [&](const int k){ return class_of(s.vclass, row[k]); };
  // Edge ends: the side normal, else the mean at the wall vertex; none at a corner
  for (int e = 0; e < R::n_edges; ++e)
    for (int o = 0; o < 2; ++o){
      const std::uint32_t cw = cls(at(at(edge_ends, e), o)), cv = cls(at(at(edge_ends, e), 1 - o));
      if (g.n_count[cw] == 0 || g.corner[cw]) continue;
      if (unit_normal(edge_end_normal(g, cw, cv).n, at(w.end, 2*e + o)))
        w.terms |= 1u << (2*e + o);
    }
  // Faces over wall edges: the facets there with this apex, else all of them
  for (int e = 0; e < R::n_edges; ++e)
    for (int j = 0; j < R::n_apex; ++j){
      const int a = at(edge_ends, e)[0], b = at(edge_ends, e)[1];
      const int v = at(at(R::apex, e), j);
      const std::uint32_t ca = cls(a), cb = cls(b);
      const SumKey k{std::min(ca, cb), std::max(ca, cb), any_apex};
      const NormalSum* es = find_sum(s.edge_sum, k);
      if (!es || g.corner[ca] || g.corner[cb]) continue;
      const NormalSum* fs = find_sum(s.face_sum, {k[0], k[1], cls(v)});
      if (unit_normal(fs ? fs->n : es->n, at(w.face, e*R::n_apex + j))){
        const int off = off_vertex<Cell>(a, b, v);
        w.triples |= std::uint16_t(1u << (e*R::n_apex + j));
        w.terms |= 1u << (f0 + R::face_slot(off, a));
        w.terms |= 1u << (f0 + R::face_slot(off, b));
      }
    }
  // 3D: the bubble's rows, the cell's wall facets without a corner vertex
  if constexpr (R::n_bubble > 0){
    int ks[nv], m = 0;
    for (int k = 0; k < nv; ++k){
      if (s.facet_neigh[i*nv + std::size_t(k)] != facet_wall) continue;
      bool at_corner = false;
      for (int j = 0; j < nv; ++j)
        if (j != k && g.corner[cls(j)]) at_corner = true;
      if (!at_corner) ks[m++] = k;
    }
    if (m > 0){
      // At most nv rows: no allocation
      Eigen::Matrix<double, Eigen::Dynamic, D, Eigen::RowMajor, nv, D> N(m, D);
      for (int r = 0; r < m; ++r){
        const Vector3d n = wall_normal<nv>(row, ks[r], facet_vertices<nv>(row, ks[r]), s.coords, s.gdim);
        for (int d = 0; d < D; ++d){
          N(r, d) = n[d];
          at(w.bn, ks[r])[std::size_t(d)] = n[d];
        }
      }
      const Eigen::Matrix<double, D, Eigen::Dynamic, 0, D, nv> P = N.completeOrthogonalDecomposition().pseudoInverse();
      for (int r = 0; r < m; ++r){
        for (int d = 0; d < D; ++d) at(w.bp, ks[r])[std::size_t(d)] = P(d, r);
        w.facets |= std::uint8_t(1u << ks[r]);
      }
      w.terms |= 1u << b0;
    }
  }
}

}  // namespace

template<typename Cell>
void build_wall_edges(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                      const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                      const std::vector<std::int32_t>& facet_neigh,
                      const std::vector<std::uint32_t>& vclass,
                      std::vector<std::int32_t>& wall_index,
                      std::vector<WallEdges<Cell>>& wall_cells, const bool verbose){
  constexpr int nv = Cell::n_verts;
  constexpr auto edge_ends = WallRule<Cell>::edge_ends;
  const Geometry<nv> g = wall_geometry<nv>(topo, coords, ncells, nverts, gdim, facet_neigh, vclass, false);

  wall_cells.clear();
  wall_cells.resize(number_wall_cells<nv>(topo, ncells, g.n_count, vclass, wall_index, verbose));
  each_wall_cell<nv>(topo, g.n_count, vclass, wall_index,
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
        if (g.n_count[cw] == 0 || g.corner[cw]) continue;
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
void upgrade_coefficients(const double* u, double* c, const UpgradeCell& w, const std::uint32_t terms,
                          const typename UpgradeTable<Cell>::Normal* n, const unsigned rest)
{
  constexpr int nv = Cell::n_verts, D = nv - 1;
  constexpr int n2 = int(Cell::n_dofs_max);
  using R = UpgradeRule<Cell>;
  constexpr auto edge_ends = WallRule<Cell>::edge_ends;
  constexpr int f0 = 2*R::n_edges, b0 = f0 + R::n_faces;

  // The P2 slope along each edge end from a wall vertex at rest, as 4 u_m - u_v
  std::array<std::array<double, D>, 2*R::n_edges> s;
  for (int e = 0; e < R::n_edges; ++e)
    for (int o = 0; o < 2; ++o){
      const int a = edge_ends[std::size_t(e)][std::size_t(o)], b = edge_ends[std::size_t(e)][std::size_t(1-o)];
      if (!(rest >> a & 1)) continue;
      for (int d = 0; d < D; ++d)
        s[std::size_t(2*e + o)][std::size_t(d)] = 4.*u[d*n2 + Cell::mid_[std::size_t(e)]] - u[d*n2 + b];
    }
  for (std::uint32_t m = terms & ~((1u << f0) - 1); m; m &= m - 1)
    for (int d = 0; d < D; ++d) c[D*__builtin_ctz(m) + d] = 0.;
  // Edge ends
  for (std::uint32_t m = terms & ((1u << f0) - 1); m; m &= m - 1, ++n){
    const int t = __builtin_ctz(m);
    double p = 0.;
    if (rest >> edge_ends[std::size_t(t/2)][std::size_t(t%2)] & 1)
      for (int d = 0; d < D; ++d) p += (*n)[std::size_t(d)]*s[std::size_t(t)][std::size_t(d)];
    for (int d = 0; d < D; ++d) c[D*t + d] = -p*(*n)[std::size_t(d)];
  }
  // Faces over wall edges
  for (unsigned m = w.triples; m; m &= m - 1, ++n){
    const int t = __builtin_ctz(m), e = t/R::n_apex;
    const int a = edge_ends[std::size_t(e)][0], b = edge_ends[std::size_t(e)][1];
    if (!(rest >> a & 1) || !(rest >> b & 1)) continue;
    const int v = R::apex[std::size_t(e)][std::size_t(t % R::n_apex)];
    const auto& sa = s[std::size_t(end_of<Cell>(a, v))];
    const auto& sb = s[std::size_t(end_of<Cell>(b, v))];
    double pa = 0., pb = 0.;
    for (int d = 0; d < D; ++d){
      pa += (*n)[std::size_t(d)]*sa[std::size_t(d)];
      pb += (*n)[std::size_t(d)]*sb[std::size_t(d)];
    }
    const int k = off_vertex<Cell>(a, b, v);
    double* ca = c + D*(f0 + R::face_slot(k, a));
    double* cb = c + D*(f0 + R::face_slot(k, b));
    for (int d = 0; d < D; ++d){
      ca[d] -= (2.*pa + pb)*(*n)[std::size_t(d)];
      cb[d] -= (pa + 2.*pb)*(*n)[std::size_t(d)];
    }
  }
  // 3D: the bubble, from the wall facets at rest
  if constexpr (R::n_bubble > 0){
    double* cb = c + D*b0;
    for (unsigned m = w.facets; m; m &= m - 1, n += 2){
      const int k = __builtin_ctz(m);
      const unsigned others = ((1u << nv) - 1) & ~(1u << k);
      if ((rest & others) != others) continue;
      const auto& bn = n[0];
      const auto& bp = n[1];
      double t = 0.;
      for (int i = 0; i < nv; ++i){
        if (i == k) continue;
        const auto& si = s[std::size_t(end_of<Cell>(i, k))];
        for (int d = 0; d < D; ++d) t += bn[std::size_t(d)]*si[std::size_t(d)];
      }
      for (int d = 0; d < D; ++d) cb[d] += -2.*t*bp[std::size_t(d)];
    }
  }
}

template<typename Cell>
void field_terms(const UpgradeTable<Cell>& upgrade, const CellDofs& u_dofs, const double* u,
                 const double tol, std::vector<double>& coef)
{
  constexpr int D = Cell::n_verts - 1;
  constexpr int nt = n_terms<Cell>;
  const std::size_t ncw = upgrade.cells.size();
  coef.resize(upgrade.n_coef);
#pragma omp parallel for schedule(dynamic, 1024)
  for (std::size_t w = 0; w < ncw; ++w){
    const UpgradeCell& wc = upgrade.cells[w];
    const TermSpan sp = upgrade.spans[w];
    std::array<double, Cell::n_dofs_max*3> uc;
    gather_cell_nodes<Cell::n_verts, Cell::n_dofs_max, D>(u_dofs[wc.cell], u_dofs.stride(), u, uc.data());
    std::array<double, D*nt> c;
    const unsigned rest = rest_vertices<Cell>(uc.data(), Cell::n_dofs_max, wc.wall, tol);
    if (rest)
      upgrade_coefficients<Cell>(uc.data(), c.data(), wc, sp.mask, upgrade.normals.data() + wc.normals, rest);
    // The cell's terms only, in their order; none at rest: zero
    double* out = coef.data() + sp.offset;
    for (std::uint32_t m = sp.mask; m; m &= m - 1, out += D)
      for (int d = 0; d < D; ++d) out[d] = rest ? c[std::size_t(D*__builtin_ctz(m) + d)] : 0.;
  }
}

template<typename Cell>
void build_wall_terms(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                      const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                      const std::vector<std::int32_t>& facet_neigh,
                      const std::vector<std::uint32_t>& vclass,
                      std::vector<std::int32_t>& wall_index,
                      UpgradeTable<Cell>& table, const bool verbose){
  constexpr int nv = Cell::n_verts, D = nv - 1;
  using Normal = typename UpgradeTable<Cell>::Normal;
  constexpr int f0 = CellTerms<Cell>::f0;
  const Geometry<nv> g = wall_geometry<nv>(topo, coords, ncells, nverts, gdim, facet_neigh, vclass, true);

  // The facets' normals summed by wall edge, and by wall edge and apex, in the facet order
  TermSources<nv> src{g, {}, {}, facet_neigh, vclass, coords, gdim};
  src.edge_sum.reserve(g.facet.size()*std::size_t((nv - 1)*(nv - 2)/2));
  src.face_sum.reserve(src.edge_sum.capacity());
  for (const auto& f : g.facet)
    for (int i = 0; i < nv - 1; ++i)
      for (int j = i + 1; j < nv - 1; ++j){
        const std::uint32_t a = at(f.key, i), b = at(f.key, j);
        const std::uint32_t lo = std::min(a, b), hi = std::max(a, b);
        src.edge_sum.push_back({{lo, hi, any_apex}, f.n});
        src.face_sum.push_back({{lo, hi, f.key[nv-1]}, f.n});
      }
  merge_sums(src.edge_sum);
  merge_sums(src.face_sum);

  // The masks, then the places in the pool and in a stamp's table, then the normals
  const std::size_t ncw = number_wall_cells<nv>(topo, ncells, g.n_count, vclass, wall_index, verbose);
  table.cells.assign(ncw, UpgradeCell{});
  table.spans.assign(ncw, TermSpan{});
  std::vector<std::uint32_t> n_normals(ncw);
  each_wall_cell<nv>(topo, g.n_count, vclass, wall_index,
                     [&](const std::size_t i, const std::size_t wc, const std::uint32_t* row, const std::uint8_t wall){
    CellTerms<Cell> w;
    cell_terms<Cell>(src, i, row, w);
    table.cells[wc] = {std::uint32_t(i), 0, w.triples, w.facets, wall};
    table.spans[wc].mask = w.terms;
    n_normals[wc] = std::uint32_t(w.n_normals());
  });
  std::size_t nn = 0, nc = 0;
  for (std::size_t w = 0; w < ncw; ++w){
    table.cells[w].normals = std::uint32_t(nn);
    table.spans[w].offset = std::uint32_t(nc);
    nn += n_normals[w];
    nc += std::size_t(D*__builtin_popcount(table.spans[w].mask));
    if (nn > std::numeric_limits<std::uint32_t>::max() || nc > std::numeric_limits<std::uint32_t>::max())
      partrac::fail(ncw, " wall cells: more wall terms than their tables can number");
  }
  table.n_coef = std::uint32_t(nc);
  std::vector<std::uint32_t>().swap(n_normals);
  table.normals.resize(nn);
  each_wall_cell<nv>(topo, g.n_count, vclass, wall_index,
                     [&](const std::size_t i, const std::size_t wc, const std::uint32_t* row, std::uint8_t){
    CellTerms<Cell> w;
    cell_terms<Cell>(src, i, row, w);
    Normal* out = table.normals.data() + table.cells[wc].normals;
    for (std::uint32_t m = w.terms & ((1u << f0) - 1); m; m &= m - 1) *out++ = at(w.end, __builtin_ctz(m));
    for (unsigned m = w.triples; m; m &= m - 1) *out++ = at(w.face, __builtin_ctz(m));
    if constexpr (UpgradeRule<Cell>::n_bubble > 0)
      for (unsigned m = w.facets; m; m &= m - 1){
        *out++ = at(w.bn, __builtin_ctz(m));
        *out++ = at(w.bp, __builtin_ctz(m));
      }
  });
}

template<typename Cell>
void build_tables(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                  const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                  const std::vector<bool>& periodic, const Vector3d& x_min, const Vector3d& x_max,
                  const double periodic_tol, const std::vector<std::int32_t>& facet_neigh,
                  const bool p2, std::vector<std::uint32_t>& vclass,
                  std::vector<std::int32_t>& wall_index, std::vector<WallEdges<Cell>>& wall_cells,
                  UpgradeTable<Cell>& upgrade, const bool verbose)
{
  if (vclass.empty())
    vclass = mesh_tables::match_periodic_vertices(coords, nverts, gdim, periodic, x_min, x_max,
                                                  periodic_tol);
  if (p2)
    build_wall_terms<Cell>(topo, coords, ncells, nverts, gdim, facet_neigh, vclass, wall_index,
                           upgrade, verbose);
  else
    build_wall_edges<Cell>(topo, coords, ncells, nverts, gdim, facet_neigh, vclass, wall_index,
                           wall_cells, verbose);
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

void label_cell_type(std::vector<int>& cell_type, const std::vector<std::int32_t>& across,
                     const Uint nv){
  const std::size_t ncells = across.size() / nv;
  cell_type.assign(ncells, 0);
  for (std::size_t i = 0; i < ncells; ++i)
    for (Uint k = 0; k < nv; ++k)
      if (across[i*nv + k] == facet_wall) cell_type[i] = 1;
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
                                         const std::vector<std::int32_t>&,
                                         const std::vector<std::uint32_t>&,
                                         std::vector<std::int32_t>&, std::vector<WallEdges<Triangle>>&,
                                         const bool);
template void build_wall_edges<Tet>(const std::vector<std::uint32_t>&, const std::vector<double>&,
                                    const std::size_t, const std::size_t, const Uint,
                                    const std::vector<std::int32_t>&,
                                    const std::vector<std::uint32_t>&,
                                    std::vector<std::int32_t>&, std::vector<WallEdges<Tet>>&,
                                    const bool);

template void upgrade_coefficients<Triangle>(const double*, double*, const UpgradeCell&, const std::uint32_t,
                                             const UpgradeTable<Triangle>::Normal*, const unsigned);
template void upgrade_coefficients<Tet>(const double*, double*, const UpgradeCell&, const std::uint32_t,
                                        const UpgradeTable<Tet>::Normal*, const unsigned);
template void field_terms<Triangle>(const UpgradeTable<Triangle>&, const CellDofs&, const double*, const double,
                                    std::vector<double>&);
template void field_terms<Tet>(const UpgradeTable<Tet>&, const CellDofs&, const double*, const double,
                               std::vector<double>&);
template void build_wall_terms<Triangle>(const std::vector<std::uint32_t>&, const std::vector<double>&,
                                         const std::size_t, const std::size_t, const Uint,
                                         const std::vector<std::int32_t>&,
                                         const std::vector<std::uint32_t>&,
                                         std::vector<std::int32_t>&, UpgradeTable<Triangle>&, const bool);
template void build_wall_terms<Tet>(const std::vector<std::uint32_t>&, const std::vector<double>&,
                                    const std::size_t, const std::size_t, const Uint,
                                    const std::vector<std::int32_t>&,
                                    const std::vector<std::uint32_t>&,
                                    std::vector<std::int32_t>&, UpgradeTable<Tet>&, const bool);

template void build_tables<Triangle>(const std::vector<std::uint32_t>&, const std::vector<double>&,
                                     const std::size_t, const std::size_t, const Uint, const std::vector<bool>&,
                                     const Vector3d&, const Vector3d&, const double,
                                     const std::vector<std::int32_t>&, const bool, std::vector<std::uint32_t>&,
                                     std::vector<std::int32_t>&, std::vector<WallEdges<Triangle>>&,
                                     UpgradeTable<Triangle>&, const bool);
template void build_tables<Tet>(const std::vector<std::uint32_t>&, const std::vector<double>&,
                                const std::size_t, const std::size_t, const Uint, const std::vector<bool>&,
                                const Vector3d&, const Vector3d&, const double,
                                const std::vector<std::int32_t>&, const bool, std::vector<std::uint32_t>&,
                                std::vector<std::int32_t>&, std::vector<WallEdges<Tet>>&,
                                UpgradeTable<Tet>&, const bool);

}  // namespace near_wall
