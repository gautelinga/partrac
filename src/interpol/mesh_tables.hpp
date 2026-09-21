#ifndef __MESH_TABLES_HPP
#define __MESH_TABLES_HPP

// The tables a simplex mesh loader needs, built from plain arrays: one parallel
// key sort, the unique edges of the topology in dolfin's local order, the values
// by node from a stored per-cell dof table, and the facet neighbours with their
// periodic partners. Nothing here knows about dolfin.

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

#include "cell_walk.hpp"
#include "typedefs.hpp"

namespace mesh_tables {

// Facet k of a cell faces vertex k; across it: a cell, a wall, or a periodic
// image, encoded as the walk encodes it
using ::facet_wall;
using ::facet_periodic;

// Local edges in dolfin's order: the vertex pairs a < b, reverse lexicographic.
// A tet gives (2,3), (1,3), (1,2), (0,3), (0,2), (0,1), which is the order
// Tet::quadbasis pairs with N[4..9]; a triangle gives (1,2), (0,2), (0,1).
template<int NV>
constexpr std::array<std::array<int, 2>, NV*(NV-1)/2> local_edges(){
  constexpr int ne = NV*(NV-1)/2;
  std::array<std::array<int, 2>, ne> e{};
  int i = 0;
  for (int a = 0; a < NV; ++a)
    for (int b = a + 1; b < NV; ++b){ e[ne - 1 - i] = {a, b}; ++i; }
  return e;
}

template<int NV> constexpr int n_local_edges = NV*(NV-1)/2;

// Stable parallel radix sort by key; equal keys keep their input order
void sort_by_key(std::vector<std::uint64_t>& key, std::vector<std::uint32_t>& payload);

// The same, with each entry's input index as its payload: the order is then
// total, so it is the same at any thread count
void sort_by_key_indexed(std::vector<std::uint64_t>& key, std::vector<std::uint32_t>& payload);

// Numbers the unique edges of an ncells x NV topology; writes the
// ncells x n_local_edges<NV> table in dolfin's local edge order and returns the
// edge count
template<int NV>
std::size_t build_edge_table(const std::vector<std::uint32_t>& topo, const std::size_t ncells,
                             std::vector<std::uint32_t>& edges);

// Where each stored dof's value belongs, in CSR: one offset per stored dof into
// the node slots (node*ncomp + c) it feeds. A dof serves one slot in an
// unconstrained space and a node with its periodic images in a reduced one, so
// a later stamp read through this reaches every image. The slots are the node
// numbering itself, so there are n_total*ncomp of them and both arrays are
// uint32 under `check_entry_count`.
struct DofNodes {
  std::vector<std::uint32_t> start;   // ndofs + 1 offsets into slot
  std::vector<std::uint32_t> slot;    // node*ncomp + c, grouped by stored dof
};

// Values by node from the stored per-cell dof table, without a scatter whose
// last writer wins: each node's (cell, slot) entries are grouped by the sort,
// the lowest cell writes and every other is compared. `cell_dofs` is
// ncells x (n_nodes*ncomp) with the components blocked, n_nodes being NV for P1
// and NV + n_local_edges<NV> for P2 (pass an empty `edges` for P1). `values`
// comes back node-major, ncomp doubles a node, and `map` sends a stored dof to
// every node slot of it. Fails on a disagreement, naming the node and the two cells.
template<int NV>
void scatter_dofs_to_nodes(const std::vector<std::uint32_t>& topo,
                           const std::vector<std::uint32_t>& edges,
                           const std::size_t ncells, const std::size_t nverts,
                           const std::size_t nedges, const std::size_t ncomp,
                           const std::vector<std::uint32_t>& cell_dofs,
                           const std::vector<double>& vec,
                           std::vector<double>& values,
                           DofNodes& map);

// The neighbour across the facet facing each vertex of each cell: a cell id, or
// facet_wall where the facet is on the boundary. The exterior facets come back
// as (midpoint, slot) in cell-then-facet order, for the periodic match.
template<int NV>
void build_facet_table(const std::vector<std::uint32_t>& topo, const std::size_t ncells,
                       const std::vector<double>& coords, const Uint gdim,
                       std::vector<std::int32_t>& across,
                       std::vector<std::pair<Vector3d, std::size_t>>& exterior);

// Periodic partners of the exterior facets: the low side of each periodic axis
// shifted onto the high side and matched within a sorted window, written into
// `across` as facet_periodic ids. Unmatched facets stay walls, and are counted
// in a message.
void match_periodic_facets(std::vector<std::int32_t>& across,
                           const std::vector<std::pair<Vector3d, std::size_t>>& exterior,
                           const std::vector<bool>& periodic,
                           const Vector3d& x_min, const Vector3d& x_max,
                           const Uint dim, const std::size_t nv, const double tol);

// The master vertex of every vertex: across each periodic axis the image on the
// low face, followed to a vertex on no high face. A P1 field written from a
// constrained space gives a vertex and its images one value, so a loader that
// numbers its nodes by vertex reads them through this. The identity where no
// axis is periodic; unmatched high-face vertices stay their own and are counted
// in a message.
std::vector<std::uint32_t> match_periodic_vertices(const std::vector<double>& coords,
                                                   const std::size_t nverts, const Uint gdim,
                                                   const std::vector<bool>& periodic,
                                                   const Vector3d& x_min, const Vector3d& x_max,
                                                   const double tol);

// The facet table with its periodic partners, from the topology alone
template<int NV>
void build_facet_neighbours(const std::vector<std::uint32_t>& topo, const std::size_t ncells,
                            const std::vector<double>& coords, const Uint gdim,
                            const std::vector<bool>& periodic,
                            const Vector3d& x_min, const Vector3d& x_max,
                            const double tol,
                            std::vector<std::int32_t>& across);

}  // namespace mesh_tables

#endif
