#ifndef __NEAR_WALL_HPP
#define __NEAR_WALL_HPP

// The near-wall P2 scheme of the XDMF loaders: where a cell's wall vertices are
// at rest, the midpoint of every edge rising from one of them is a fixed linear
// map of the velocity at the edge's fluid end, so u.n grows as delta^2 above
// the facet. The map differs by dimension -- a matrix in 2D, u_v/2 + (u_v.n) q
// in 3D -- so the per-end record and the table are specialised. Nothing here
// needs dolfin: the wall facets are the facet table's walls and their normals
// come from the coordinates.

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include "Tet.hpp"
#include "Triangle.hpp"
#include "typedefs.hpp"

namespace near_wall {

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
  std::uint8_t wall;   // bit k: vertex k lies on a wall
  WallEdges() {}       // unset: build_wall_edges fills the table in parallel
};

// The near-wall rule, per dimension: the P2 block of a P1 one whose wall
// vertices are at rest, false where none is
template<typename Cell>
bool wall_block(const double* u, double* u2, const WallEdges<Cell>& w, const double tol);

template<>
inline bool wall_block<Triangle>(const double* u, double* u2, const WallEdges<Triangle>& w,
                          const double tol)
{
  constexpr int n1 = 3;
  // Walls: listed vertices at rest
  unsigned rest = 0;
  for (int k = 0; k < 3; ++k)
    if ((w.wall >> k & 1) && std::abs(u[k]) <= tol && std::abs(u[n1 + k]) <= tol)
      rest |= 1u << k;
  if (rest == 0) return false;

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
  return true;
}

template<>
inline bool wall_block<Tet>(const double* u, double* u2, const WallEdges<Tet>& w,
                     const double tol)
{
  constexpr int n1 = 4;
  // Walls: listed vertices at rest
  unsigned rest = 0;
  for (int k = 0; k < 4; ++k)
    if ((w.wall >> k & 1) && std::abs(u[k]) <= tol && std::abs(u[n1 + k]) <= tol
        && std::abs(u[2*n1 + k]) <= tol)
      rest |= 1u << k;
  if (rest == 0) return false;

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
  return true;
}

// The wall cells and their edge ends, from the mesh arrays and the facet table.
// A wall facet is a facet_wall entry, and its normal is the outward one the
// vertex it faces fixes the sign of; the facets are visited in the order of
// their sorted vertex ids, which is the order the entity numbering gives them,
// so no sum depends on how the cells are numbered. vclass sends a vertex to
// the master of its periodic images; verbose prints the wall cell count.
template<typename Cell>
void build_wall_edges(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                      const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                      const std::vector<std::int32_t>& facet_neigh,
                      const std::vector<std::uint32_t>& vclass,
                      std::vector<std::int32_t>& wall_index,
                      std::vector<WallEdges<Cell>>& wall_cells, const bool verbose);

// 0 bulk, 1 on a wall, 2 next to a cell on a wall; from the facet table
void label_cell_type(std::vector<int>& cell_type, const std::vector<std::int32_t>& across,
                     const Uint nv);

}  // namespace near_wall

#endif
