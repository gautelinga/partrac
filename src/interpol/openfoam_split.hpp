#ifndef __OPENFOAM_SPLIT_HPP
#define __OPENFOAM_SPLIT_HPP

// A conforming simplex split of an OpenFOAM polyMesh. Every face is fanned
// from findFaceBasePts' base point (the owner's, reused across a cyclic).
// W12: each cell the fan from its centre over its faces' triangles,
// cellPoint's own tets; W6: a hex Dompierre's 5 or 6 tets where its face
// diagonals meet at a corner, else the fan. A cell with a flat tet takes the
// other split. 2D: the front plane of an empty pair. Needs no OpenFOAM.

#include <cstdint>
#include <vector>

#include "openfoam_load.hpp"

namespace openfoam_split {

// Smallest simplex, relative to its cell
constexpr double valid_fraction = 1e-12;

struct SplitData {
  int nv = 4;                               // 4: tets, 3: triangles
  int gdim = 3;
  std::vector<std::uint32_t> cells;         // nv a simplex into the nodes, positively oriented, by cell
  std::vector<double> node_x;               // gdim a node
  std::vector<std::int8_t> node_kind;       // 0 a mesh point, 1 a cell centre
  std::vector<std::int32_t> node_point;     // the mesh point, -1 for a centre
  std::vector<std::int32_t> node_cell;      // the cell of a centre, -1 for a point
  std::vector<std::int32_t> cell_of;        // the OpenFOAM cell of each simplex
  std::vector<std::int32_t> facet_patch;    // nv a simplex: the facet facing each vertex, -1 inside
  std::vector<std::uint32_t> node_master;   // the lowest node among a node's cyclic images
  std::vector<std::int64_t> facet_partner;  // nv a simplex: the facet slot across a cyclic, else -1
  std::vector<int> inplane;                 // 2D: the two in-plane axes, ascending
  // What the split did
  std::size_t fan_cells = 0, five_tet_hexes = 0, fallback = 0;
  std::size_t invalid_cells = 0, invalid_simplices = 0, base_failures = 0;
  std::size_t nsimplices() const { return cell_of.size(); }
  std::size_t nnodes() const { return node_kind.size(); }
};

// The split of a case, W12 or W6; a case with an empty pair gives triangles
SplitData split(const openfoam_load::CaseData& c, int tets_per_hex);

// Each face's base point by findFaceBasePts' rule; where no base is good the
// best, counted in failures; only the faces marked in only, if given, the
// rest -1
std::vector<std::int32_t> face_bases(const openfoam_load::CaseData& c, std::size_t* failures = nullptr,
                                     const std::vector<char>* only = nullptr);

// The smallest cyclic image of every mesh point, the point itself without one
std::vector<std::int32_t> cyclic_masters(const openfoam_load::CaseData& c);

}  // namespace openfoam_split

#endif
