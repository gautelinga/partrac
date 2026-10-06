#ifndef __OPENFOAMINTERPOL_HPP
#define __OPENFOAMINTERPOL_HPP

// An OpenFOAM case traced in place: the polyMesh split into tets (or, one cell
// thick between an empty pair, its front plane into triangles), each time
// directory a stamp, the cyclic pairs the periodic box. The velocity and the
// pressure are P1 on the split, their node values by least squares (exact for
// linear fields) or cellPoint's by inverse distance (openfoam_nodes.hpp), each
// field's operator built at load from its first stamp's conditions. A phase
// field is P1 on cellPoint's node values, each fanned cell's centre moved so
// its volume above 1/2 is its volume fraction's (openfoam_phase.hpp); its
// gradient is P1 on the least-squares fit's. Every boundary facet but a
// cyclic one's is a wall, whatever its patch. The files are read through
// openfoam_load.hpp; nothing here includes OpenFOAM.

#include <string>
#include <vector>

#include "Params.hpp"
#include "Timestamps.hpp"
#include "openfoam_load.hpp"
#include "openfoam_nodes.hpp"
#include "openfoam_phase.hpp"

template<typename Cell, typename Format> class StampedInterpol;

// A time directory per stamp
struct OpenFoamFormat {
  using Key = std::string;   // the time directory
  static constexpr bool vertex_fields = true;
  static partrac::Schema schema(int D);

  Timestamps ts;
  std::string u_field, p_field, phi_field;   // phi_field empty: no phase field
  // What a field's read needs of the case: its directory, cells and patches
  openfoam_load::CaseData layout;
  bool least_squares = true;
  openfoam_nodes::LeastSquares w_u, w_p;
  openfoam_nodes::InverseDistance idw;
  openfoam_nodes::Gradient g_phi;
  openfoam_phase::Cells phase_cells;   // by the split's rows and vertices
  openfoam_load::FieldData first_u, first_p, first_phi;   // the first stamp's, read at load, until its read
  std::vector<std::uint32_t> slot;   // a split node's vertex, after the renumbering
  std::vector<int> comps;            // the velocity components the cells take
  std::vector<openfoam_nodes::PatchClass> patch_class;    // from the first stamp's velocity

  // The case, its split, the tables and the first stamp
  template<typename I> void load(I& intp, const std::string& infilename);
  // The boundary facets' count by class, and the patches of each
  void report_boundary(const std::vector<openfoam_load::Patch>& patches,
                       const std::vector<std::int32_t>& facet_patch) const;
  // One time directory's fields by vertex
  template<typename I> void read(I& intp, const Key& time, typename I::Stamp& s);
  template<typename I> Key key(const I&, const Stamp& st) const { return st.filename; }
  std::string name(const Stamp& st) const { return layout.dir + "/" + st.filename; }
};

#include "StampedInterpol.hpp"

template<typename Cell>
using OpenFoamInterpol = StampedInterpol<Cell, OpenFoamFormat>;

#endif
