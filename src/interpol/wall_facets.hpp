#ifndef __WALL_FACETS_HPP
#define __WALL_FACETS_HPP

// The no-slip walls a facet file names (wall_facets, wall_markers): a
// meshio-style XDMF of facets, triangles in 3D and segments in 2D, each with
// a marker, matched to the mesh's exterior facets by the coordinates of their
// vertices. Nothing here needs dolfin.

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "typedefs.hpp"

namespace wall_facets {

// The markers a list selects: values and inclusive ranges, comma-separated,
// as 7:106 or 3,5,7:106; refuses an empty list
class Markers {
public:
  explicit Markers(const std::string& list);
  bool selects(const std::int64_t m) const;
private:
  std::vector<std::pair<std::int64_t, std::int64_t>> ranges_;
};

// A facet file: points (gdim a point), facets (gdim points a facet, as ids
// into points) and their markers, empty where the file holds none (refused by match)
struct File {
  std::vector<double> points;
  std::vector<std::int64_t> facets;
  std::vector<std::int64_t> markers;
};
File read(const std::string& xdmf, const Uint gdim);

// The selected facets of f as the mesh's, by cell: bit k the facet facing
// vertex k. A point of a selected facet maps to the boundary vertex within
// 1e-3 of the shortest boundary edge there, a facet to the exterior facet of
// those vertices. Refuses a point that matches no vertex, a facet that matches
// no exterior facet, and one on a periodic side; name is the file's in the
// messages; verbose prints the count
std::vector<std::uint8_t> match(const File& f, const Markers& m, const std::vector<std::uint32_t>& topo,
                                const std::vector<double>& coords, const std::size_t ncells,
                                const Uint gdim, const std::vector<std::int32_t>& facet_neigh,
                                const std::string& name, const bool verbose);

// read and match the file the parameter wall_facets names, in folder, with
// the markers of wall_markers
std::vector<std::uint8_t> load(const std::string& folder, const std::string& file, const std::string& markers,
                               const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                               const std::size_t ncells, const Uint gdim,
                               const std::vector<std::int32_t>& facet_neigh, const bool verbose);

}  // namespace wall_facets

#endif
