#ifndef __XDMFINTERPOL_HPP
#define __XDMFINTERPOL_HPP

// P1 velocity, pressure and phase fields written as XDMF, on triangles or
// tets. The dofs of a P1 field are the mesh vertices, so the values are read
// by vertex and every table comes from the mesh arrays the XDMF names; nothing
// here needs dolfin.

#include <string>
#include <vector>

#include "Params.hpp"
#include "Timestamps.hpp"
#include "typedefs.hpp"

template<typename Cell, typename Format> class StampedInterpol;

// An XDMF file per field, naming a dataset per stamp; the velocity's carries the mesh
struct XDMFFormat {
  using Key = Uint;   // the stamp's index in the timestamps
  static constexpr bool vertex_fields = true;
  static partrac::Schema schema(int D);

  MultiTimestamps ts;

  // The mesh, the vertex classes and the node table; the stamps come with update
  template<typename I> void load(I& intp, const std::string& infilename);
  // One stamp's fields by vertex
  template<typename I> void read(I& intp, const Key& it, typename I::Stamp& s);
  template<typename I> Key key(const I&, const MultiStamp& st) const { return st.it; }
  std::string name(const MultiStamp& st) const {
    const auto path = ts.get_path("u", st.it);
    return path[0] + ":" + path[1];
  }
};

#include "StampedInterpol.hpp"

template<typename Cell>
using XDMFInterpol = StampedInterpol<Cell, XDMFFormat>;

#endif
