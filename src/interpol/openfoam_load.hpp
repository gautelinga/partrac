#ifndef __OPENFOAM_LOAD_HPP
#define __OPENFOAM_LOAD_HPP

// An OpenFOAM case as plain arrays: the polyMesh, its patches, the geometry
// OpenFOAM computes, the time directories, and one field of one time with
// each patch's condition. OpenFOAM's own libraries read the files, from one
// translation unit (openfoam_load.cpp) built as a plugin that is loaded only
// when a case is read (openfoam_host.cpp), so no other code and no app links
// OpenFOAM or sees one of its types. Every refusal comes back as a
// partrac::Error.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace openfoam_load {

// A boundary patch as constant/polyMesh/boundary has it
struct Patch {
  std::string name, type;
  std::int64_t start = 0, size = 0;       // its faces are start .. start + size
  std::vector<std::string> groups;
  // A cyclic's other side: face i here is face i there, its points reversed
  int neighbour = -1;
  bool owner = false;                     // the first of the pair
  std::array<double, 3> separation{};     // x there = x here + separation
};

// The polyMesh and what OpenFOAM derives from it; faces point out of their owner
struct CaseData {
  std::string dir;
  std::vector<double> points;                 // 3 a point
  std::vector<std::int64_t> face_start;       // nfaces + 1 offsets into face_points
  std::vector<std::int32_t> face_points;
  std::vector<std::int32_t> owner;            // nfaces
  std::vector<std::int32_t> neighbour;        // n_internal
  std::int64_t n_internal = 0, ncells = 0;
  std::vector<Patch> patches;
  std::vector<double> cell_centres;           // 3 a cell, primitiveMesh's
  std::vector<double> cell_volumes;
  std::vector<double> bface_centres;          // 3 a boundary face, face f at f - n_internal
  std::vector<double> bface_areas;            // area vectors, out of the owner
  // The time directories holding the velocity, by value; constant and 0
  // excluded unless 0 is the only one
  std::vector<std::string> times;
  std::vector<double> time_values;
  int empty_axis = -1;                        // the normal of a 2D case's empty pair, else -1
  std::size_t nfaces() const { return owner.size(); }
  std::size_t npoints() const { return points.size()/3; }
};

// The diagonal of a case's bounding box
inline double box_diagonal(const CaseData& c){
  std::array<double, 3> lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
  for (std::size_t i = 0; i < c.npoints(); ++i)
    for (std::size_t d = 0; d < 3; ++d){
      lo[d] = std::min(lo[d], c.points[3*i + d]);
      hi[d] = std::max(hi[d], c.points[3*i + d]);
    }
  double s = 0.;
  for (std::size_t d = 0; d < 3; ++d) s += (hi[d] - lo[d])*(hi[d] - lo[d]);
  return std::sqrt(s);
}

// One patch of a field: its condition, and the values where it holds them
struct FieldPatch {
  std::string condition;
  bool fixes_value = false;           // fvPatchField::fixesValue()
  bool has_value = false;
  bool unknown = false;               // OpenFOAM cannot build it: its library is not loaded
  std::vector<double> values;         // ncomp a face
};

// One field of one time directory
struct FieldData {
  std::string name, time;
  int ncomp = 1;
  std::vector<double> internal;       // ncomp a cell
  std::vector<FieldPatch> patches;    // as CaseData::patches
};

// The plugin's interface, its version and the sizes of what its entry points
// pass: the host refuses a plugin whose stamp is not its own
constexpr int layout_version = 2;
constexpr std::uint64_t layout_stamp =
  ((((std::uint64_t(layout_version)*1000 + sizeof(CaseData))*1000 + sizeof(FieldData))*1000
    + sizeof(Patch))*1000 + sizeof(FieldPatch))*1000 + sizeof(std::string);

// The mesh of a case; velocity_field names the field whose directories are the times
CaseData read_case(const std::string& dir, const std::string& velocity_field);

// A field of one time directory, the file alone read; the patches are the case's
FieldData read_field(const CaseData& c, const std::string& time, const std::string& name);

// Whether a case has empty patches, from its boundary file alone
bool has_empty_patches(const std::string& dir);

// Whether this build reads OpenFOAM cases
bool available();

}  // namespace openfoam_load

#endif
