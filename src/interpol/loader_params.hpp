#ifndef __LOADER_PARAMS_HPP
#define __LOADER_PARAMS_HPP

// Keys of the parameter files the field loaders read

#include <string>
#include "Params.hpp"

// Every mesh loader
inline void add_mesh_params(partrac::Schema& s, const bool three_d){
  s.opt<bool>("periodic_x", false, "periodic along x");
  s.opt<bool>("periodic_y", false, "periodic along y");
  s.opt<bool>("periodic_z", false, three_d ? "periodic along z" : "periodic along z; not read in 2D");
  s.optional<double>("rho", "density");
  s.opt<bool>("ignore_pressure", false, "do not read the pressure");
}

// The datasets the field files hold
inline void add_field_names(partrac::Schema& s){
  s.opt<std::string>("velocity_field", "u", "velocity dataset in the field files");
  s.opt<std::string>("pressure_field", "p", "pressure dataset in the field files");
}

// Fields in dolfin HDF5 files: mode tet, triangle or fenics
inline partrac::Schema dolfin_h5_schema(const std::string& mode){
  partrac::Schema s("dolfin_params.dat, mode=" + mode, "");
  add_mesh_params(s, mode != "triangle");
  s.require<std::string>("timestamps", "file listing each stamp's time and field file");
  s.require<std::string>("mesh", "mesh file");
  s.require<std::string>("velocity_space", "velocity element, as P1 or P2");
  s.require<std::string>("pressure_space", "pressure element, as P1 or P2");
  if (mode == "fenics"){
    // Read by DolfInterpol only
    s.opt<std::string>("renumber_cells", "auto", "renumber cells by their dofs: auto (if poorly ordered), never, always");
    s.choices("renumber_cells", {"auto", "never", "always"});
  }
  else {
    // Read by SimplexInterpol only
    s.opt<bool>("mesh_cache", false, "keep the loaded tables beside the mesh and read them back");
  }
  add_field_names(s);
  return s;
}

// A time series as frequency components: mode trianglefreq
inline partrac::Schema triangle_freq_schema(){
  partrac::Schema s("dolfin_params.dat, mode=trianglefreq", "");
  add_mesh_params(s, false);
  s.require<std::string>("freqstamps", "file listing each component's time shift, amplitude and field file");
  s.require<std::string>("mesh", "mesh file");
  s.require<std::string>("velocity_space", "velocity element, as P1 or P2");
  s.require<std::string>("pressure_space", "pressure element, as P1 or P2");
  add_field_names(s);
  s.require<double>("tau", "base period; 0 or less for none");
  s.require<double>("t_min", "start of the time interval");
  s.require<double>("t_max", "end of the time interval");
  return s;
}

// P1 fields in XDMF files: mode xdmftriangle or xdmftet
inline partrac::Schema xdmf_schema(const std::string& mode){
  partrac::Schema s("dolfin_params.dat, mode=" + mode, "");
  add_mesh_params(s, mode == "xdmftet");
  s.require<std::string>("u", "velocity XDMF file");
  s.require_if<std::string>("p", [](const partrac::Params& p){ return !p.get<bool>("ignore_pressure"); },
                            "ignore_pressure is false", "pressure XDMF file");
  s.opt<bool>("include_phi", false, "read a phase field");
  s.require_if<std::string>("phi", [](const partrac::Params& p){ return p.get<bool>("include_phi"); },
                            "include_phi is true", "phase field XDMF file");
  s.opt<bool>("include_pf", false, "not read");
  s.opt<std::string>("wall_p2", "edge", "velocity next to walls at rest: edge (quadratic) or none (P1)");
  s.choices("wall_p2", {"edge", "none"});
  return s;
}

// Lattice Boltzmann fields: mode felbm (structured, lbm)
inline partrac::Schema felbm_schema(){
  partrac::Schema s("felbm_params.dat", "");
  s.opt<std::string>("timestamps", "timestamps.dat", "file listing each stamp's time and field file");
  s.opt<std::string>("is_solid_file", "output_is_solid.h5", "solid mask");
  s.opt<bool>("ignore_pressure", false, "do not read the pressure");
  s.opt<bool>("ignore_density", false, "do not read the density");
  s.opt<bool>("ignore_uz", false, "do not read u_z");
  s.opt<std::string>("interpolation", "linear",
                     "in space: linear (trilinear, no slip at the walls) or constant (the nearest node, no gradient)");
  s.choices("interpolation", {"linear", "constant"});
  return s;
}

#endif
