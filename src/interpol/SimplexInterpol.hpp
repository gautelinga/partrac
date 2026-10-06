#ifndef __SIMPLEXINTERPOL_HPP
#define __SIMPLEXINTERPOL_HPP

// Velocity, pressure and phase stamps written as dolfin HDF5 checkpoints, on
// triangles or tets: P1 or P2 fields read from the stored dofmaps and the mesh
// arrays, and evaluated through this code's own cell walk, cell tree and
// basis. Nothing here needs dolfin.

#include <string>

#include "Params.hpp"
#include "mesh_tables.hpp"
#include "Timestamps.hpp"

template<typename Cell, typename Format> class StampedInterpol;

// A file per stamp, listed with its time in the timestamps file; the mesh in a
// file of its own
struct DolfinH5Format {
  using Key = std::string;   // the stamp's file
  static constexpr bool vertex_fields = false;
  static partrac::Schema schema(int D);

  Timestamps ts;
  std::string u_field, p_field, phi_field;
  // Stored dof -> every node slot it feeds, for the later stamps
  mesh_tables::DofNodes u_map, p_map, phi_map;

  // The mesh, the tables and the first stamp, or all of it from the mesh cache
  template<typename I> void load(I& intp, const std::string& infilename);
  // One stamp's vectors into a buffer, through the retained dof mappings
  template<typename I> void read(I& intp, const Key& file, typename I::Stamp& s);
  template<typename I> Key key(const I& intp, const Stamp& st) const {
    return intp.get_folder() + "/" + st.filename;
  }
  std::string name(const Stamp& st) const { return st.filename; }
};

#include "StampedInterpol.hpp"

template<typename Cell>
using SimplexInterpol = StampedInterpol<Cell, DolfinH5Format>;

#endif
