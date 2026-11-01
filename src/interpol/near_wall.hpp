#ifndef __NEAR_WALL_HPP
#define __NEAR_WALL_HPP

// The near-wall rules, which raise the velocity's order next to walls at rest.
// A P1 velocity becomes P2 by the edge rule in the cells with a wall vertex:
// the midpoint of every edge rising from a wall vertex is a fixed linear map
// of the velocity at the edge's fluid end, so u.n grows as delta^2 above the
// facet. The map differs by dimension -- a matrix in 2D, u_v/2 + (u_v.n) q in
// 3D -- so the per-end record and the table are specialised. A P2 velocity
// gains, in the cells with a wall vertex at rest, the wall terms and the
// correction of wall_split.hpp. Nothing here needs dolfin: the wall facets are
// the facet table's walls or those wall_facets names (wall_facets.hpp), and
// their normals come from the coordinates.

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "Tet.hpp"
#include "Triangle.hpp"
#include "typedefs.hpp"
#include "wall_split.hpp"

class CellDofs;

namespace near_wall {

// Whether wall_upgrade reads the velocity as P2, else P1
template<typename Cell>
bool velocity_p2(const Uint ncoeffs){ return ncoeffs == Uint(Cell::n_dofs_max); }

// Warns, when there are wall cells, that their walls come from the velocity
void warn_velocity_walls(const std::vector<std::int32_t>& wall_index);
// The rest tolerance of a velocity field: round-off of the largest value a cell can read
double rest_tolerance(const std::vector<double>& u);
// The same for a field on the vertex table, where an image vertex's own value
// is never gathered, its master's is
double rest_tolerance_by_vertex(const std::vector<double>& u, const std::vector<std::uint32_t>& vclass,
                                const std::size_t nverts, const Uint dim);

// What an edge end stores, and how many ends a cell has
template<typename Cell> struct WallRule;

template<> struct WallRule<Triangle> {
  // Midpoint of a wall-vertex edge: M u_v
  struct End { double mxx, mxy, myx, myy; };
  // Ends by edge (01, 02, 12), then (first, second)
  static constexpr int n_ends = 6;
  static constexpr std::array<std::array<int, 2>, 3> edge_ends = {{{0, 1}, {0, 2}, {1, 2}}};
};

template<> struct WallRule<Tet> {
  // Midpoint of a wall-vertex edge: u_v/2 + (u_v.n) q
  struct End { double qx, qy, qz, nx, ny, nz; };
  // Ends by edge (01, 02, 03, 12, 13, 23), then (first, second)
  static constexpr int n_ends = 12;
  static constexpr std::array<std::array<int, 2>, 6> edge_ends =
    {{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
};

template<typename Cell>
struct WallEdges {
  std::array<typename WallRule<Cell>::End, WallRule<Cell>::n_ends> ends;
  std::uint8_t wall;   // bit k: vertex k lies on a wall (stamped: then at rest)
  WallEdges() {}       // unset: build_wall_edges fills the table in parallel
};

// A block's wall vertices at rest, by bit: those of wall whose values are all
// within tol, n values a direction
template<typename Cell>
inline unsigned rest_vertices(const double* u, const Uint n, const std::uint8_t wall, const double tol)
{
  unsigned rest = 0;
  for (int k = 0; k < Cell::n_verts; ++k){
    if (!(wall >> k & 1)) continue;
    bool r = true;
    for (int d = 0; d < Cell::n_verts - 1; ++d) r = r && std::abs(u[d*n + k]) <= tol;
    if (r) rest |= 1u << k;
  }
  return rest;
}

// Raise a P1 block, per dimension; rest by bit
template<typename Cell>
void raise_block(const double* u, double* u2, const WallEdges<Cell>& w, const unsigned rest);

template<>
inline void raise_block<Triangle>(const double* u, double* u2, const WallEdges<Triangle>& w,
                                  const unsigned rest)
{
  constexpr int n1 = 3;
  constexpr auto edge_ends = WallRule<Triangle>::edge_ends;
  for (int k = 0; k < 3; ++k){
    const bool r = rest >> k & 1;
    u2[k] = r ? 0. : u[k];
    u2[6 + k] = r ? 0. : u[n1 + k];
  }
  for (int e = 0; e < 3; ++e){
    const int a = edge_ends[e][0], c = edge_ends[e][1];
    const int m = Triangle::mid_[e];
    const bool wa = rest >> a & 1, wc = rest >> c & 1;
    const WallRule<Triangle>::End& we = w.ends[2*e + (wa ? 0 : 1)];
    if (wa != wc){
      const int v = wa ? c : a;
      const double ux = u[v], uy = u[n1 + v];
      u2[m] = we.mxx*ux + we.mxy*uy;
      u2[6 + m] = we.myx*ux + we.myy*uy;
    }
    else {
      u2[m] = 0.5*(u2[a] + u2[c]);
      u2[6 + m] = 0.5*(u2[6 + a] + u2[6 + c]);
    }
  }
}

template<>
inline void raise_block<Tet>(const double* u, double* u2, const WallEdges<Tet>& w, const unsigned rest)
{
  constexpr int n1 = 4;
  constexpr auto edge_ends = WallRule<Tet>::edge_ends;
  for (int k = 0; k < 4; ++k){
    const bool r = rest >> k & 1;
    for (int c = 0; c < 3; ++c)
      u2[10*c + k] = r ? 0. : u[c*n1 + k];
  }
  for (int e = 0; e < 6; ++e){
    const int a = edge_ends[e][0], b = edge_ends[e][1];
    const int m = Tet::mid_[e];
    const bool wa = rest >> a & 1, wb = rest >> b & 1;
    const WallRule<Tet>::End& we = w.ends[2*e + (wa ? 0 : 1)];
    if (wa != wb){
      const int v = wa ? b : a;
      const double ux = u[v], uy = u[n1 + v], uz = u[2*n1 + v];
      const double un = ux*we.nx + uy*we.ny + uz*we.nz;
      u2[m] = 0.5*ux + un*we.qx;
      u2[10 + m] = 0.5*uy + un*we.qy;
      u2[20 + m] = 0.5*uz + un*we.qz;
    }
    else {
      for (int c = 0; c < 3; ++c)
        u2[10*c + m] = 0.5*(u2[10*c + a] + u2[10*c + b]);
    }
  }
}

// The wall facets, here and below: wall_facets' bits by cell (bit k: the
// facet facing vertex k), else facet_neigh's facet_wall entries.
// The P2 wall cells, numbered in cell order, and the candidates of their
// rest bits (wall_split::Shape): the cell's wall facets, its edges (3D) or
// vertices (2D) on a wall facet of the mesh, and in 3D its wall vertices, by
// vertex class; verbose prints their count
template<typename Cell>
void build_wall_patterns(const std::vector<std::uint32_t>& topo, const std::size_t ncells,
                         const std::size_t nverts, const std::vector<std::int32_t>& facet_neigh,
                         const std::vector<std::uint8_t>& wall_facets, const std::vector<std::uint32_t>& vclass,
                         std::vector<std::int32_t>& wall_index, std::vector<wall_split::WallCell>& walls,
                         const bool verbose);

// Once the P2 wall cells' rest bits are fixed: those with none leave (wall
// index -1), the others are renumbered in cell order and get their wall terms
// and normals, and the tables of their patterns are built
template<typename Cell>
void build_wall_terms(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                      const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                      const std::vector<std::int32_t>& facet_neigh, const std::vector<std::uint8_t>& wall_facets,
                      const std::vector<std::uint32_t>& vclass, std::vector<std::int32_t>& wall_index,
                      wall_split::Tables<Cell::n_verts>& t, const bool verbose);

// The wall cells and their edge ends, from the mesh arrays and the wall
// facets. A wall facet's normal is the outward one the
// vertex it faces fixes the sign of; the facets are visited in the order of
// their sorted vertex ids, which is the order the entity numbering gives them,
// so no sum depends on how the cells are numbered. vclass sends a vertex to
// the master of its periodic images; verbose prints the wall cell count.
template<typename Cell>
void build_wall_edges(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                      const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                      const std::vector<std::int32_t>& facet_neigh,
                      const std::vector<std::uint8_t>& wall_facets,
                      const std::vector<std::uint32_t>& vclass,
                      std::vector<std::int32_t>& wall_index,
                      std::vector<WallEdges<Cell>>& wall_cells, const bool verbose);

// The wall tables: vclass where the format left it empty (P2: only with a
// periodic axis); for a P2 velocity the wall cells and their rest candidates, else
// each cell's wall index (-1: no wall vertex) and the P1 edge ends
template<typename Cell>
void build_tables(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                  const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                  const std::vector<bool>& periodic, const Vector3d& x_min, const Vector3d& x_max,
                  const double periodic_tol, const std::vector<std::int32_t>& facet_neigh,
                  const std::vector<std::uint8_t>& wall_facets, const bool p2, std::vector<std::uint32_t>& vclass,
                  std::vector<std::int32_t>& wall_index, std::vector<WallEdges<Cell>>& wall_cells,
                  std::vector<wall_split::WallCell>& walls, const bool verbose);

// The nodes of the walls wall_facets names against each field of us with
// flow, or with amps their sum weighted by |amps|: warns once on stderr of
// those whose |u| (bounded by sum |a_k||u_k|) exceeds 1e-6 of the largest
// (sum |a_k| max |u_k|), which stay no-slip walls; where names the fields.
// The nodes: a P2 wall cell's entities at rest (walls), else a P1 one's wall
// vertices (edges); the largest |u| by vertex where vclass is given. False
// where no field has flow
template<typename Cell>
bool check_moving_walls(const CellDofs& dofs, const std::vector<std::int32_t>& wall_index, const bool p2,
                        const std::vector<wall_split::WallCell>& walls, const std::vector<WallEdges<Cell>>& edges,
                        const std::vector<const std::vector<double>*>& us, const std::vector<std::uint32_t>& vclass,
                        const std::string& where, const std::vector<double>& amps = {});

// 0 bulk, 1 on a wall, 2 next to a cell on a wall; the walls are wall_facets'
// (walls, by cell) where given, else the facet table's
void label_cell_type(std::vector<int>& cell_type, const std::vector<std::int32_t>& across,
                     const Uint nv, const std::vector<std::uint8_t>& walls = {});

}  // namespace near_wall

#endif
