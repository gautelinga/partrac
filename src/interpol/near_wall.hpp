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
};

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
