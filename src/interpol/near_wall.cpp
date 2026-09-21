#include "near_wall.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <utility>

#include "cell_walk.hpp"

namespace near_wall {

namespace {

// A vertex as a point, the third coordinate zero in 2D
Vector3d point_of(const std::vector<double>& coords, const std::size_t v, const Uint gdim){
  Vector3d p = Vector3d::Zero();
  for (Uint d = 0; d < gdim; ++d) p[d] = coords[v*gdim + d];
  return p;
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

// Outward unit normal of the facet facing vertex k, in the arithmetic the
// element libraries use: the cross product of two facet edges in 3D, the part
// of an edge across the facet in 2D, signed by the vertex it faces
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
  const auto klass = [&](const std::uint32_t v){ return std::size_t(vclass.empty() ? v : vclass[v]); };

  // The wall facets, ordered by their vertex ids: a facet's place in that order
  // is its place in the entity numbering, so a vertex sums its normals in one
  // sequence whatever the cells are numbered by
  std::vector<std::size_t> slots;
  for (std::size_t s = 0; s < ncells*nv; ++s)
    if (facet_neigh[s] == facet_wall) slots.push_back(s);
  std::sort(slots.begin(), slots.end(), [&](const std::size_t a, const std::size_t b){
    return facet_vertices<nv>(topo.data() + (a/nv)*nv, int(a % nv))
         < facet_vertices<nv>(topo.data() + (b/nv)*nv, int(b % nv));
  });

  // Wall normals at the vertices, and the side normals per fluid end and wall
  // end, summed over the facets that meet there in that one order
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
      else if ((nv == 3 && n_count[iv] > 1) || n_first[iv].dot(n) < 0.5) corner[iv] = true;   // over 60 degrees
      n_sum[iv] += n;
      ++n_count[iv];
      // the facet's own cell holds the apex it faces
      side_normal.emplace(std::make_pair(klass(row[k]), klass(row[j])),
                          Vector3d::Zero().eval()).first->second += n;
    }
  }

  wall_index.assign(ncells, -1);
  wall_cells.clear();
  for (std::size_t i = 0; i < ncells; ++i){
    const std::uint32_t* row = topo.data() + i*nv;
    WallEdges<Cell> w{};
    for (int k = 0; k < nv; ++k)
      if (n_count[klass(row[k])] > 0) w.wall |= 1 << k;
    if (w.wall == 0) continue;

    for (int e = 0; e < nv*(nv-1)/2; ++e){
      for (int o = 0; o < 2; ++o){
        typename WallRule<Cell>::End& we = w.ends[2*e + o];
        const std::uint32_t iw = row[edge_ends[std::size_t(e)][std::size_t(o)]];
        const std::uint32_t iv = row[edge_ends[std::size_t(e)][std::size_t(1-o)]];
        const std::size_t cw = klass(iw);
        if constexpr (nv == 3) we = {0.5, 0., 0., 0.5};   // linear
        else                   we = {0., 0., 0., 0., 0., 0.};
        if (n_count[cw] == 0 || corner[cw]) continue;
        const Vector3d edge = point_of(coords, iv, gdim) - point_of(coords, iw, gdim);
        // over several wall facets: their mean normal
        const auto side = side_normal.find({klass(iv), cw});
        const Vector3d n = (side != side_normal.end() ? side->second : n_sum[cw]).normalized();
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
    wall_index[i] = std::int32_t(wall_cells.size());
    wall_cells.push_back(w);
  }
  if (verbose)
    std::cout << "Wall cells: " << wall_cells.size() << std::endl;
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

}  // namespace near_wall
