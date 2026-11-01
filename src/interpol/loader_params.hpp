#ifndef __LOADER_PARAMS_HPP
#define __LOADER_PARAMS_HPP

// Keys of the parameter files the field loaders read

#include <string>
#include "Params.hpp"

// The near-wall upgrade; wall_p2 is refused, naming it
inline void add_wall_upgrade(partrac::Schema& s, const bool def, const std::string& doc){
  s.opt<bool>("wall_upgrade", def, doc);
  s.replaced("wall_p2", "wall_upgrade = true (wall_p2 = edge) or false (wall_p2 = none)");
}

// The no-slip walls as a facet file names them, for wall_upgrade
inline void add_wall_facets(partrac::Schema& s){
  s.optional<std::string>("wall_facets", "XDMF facet file (meshio's, markers per facet) naming the no-slip walls "
                                         "wall_upgrade raises the velocity at, all at rest; without it the "
                                         "exterior facets off the periodic sides at rest in the velocity");
  s.optional<std::string>("wall_markers", "markers of wall_facets' no-slip walls, values and inclusive ranges as "
                                          "3,5,7:106; required with wall_facets");
  s.check([](const partrac::Params& p){ return !p.has("wall_markers") || p.has("wall_facets"); },
          "wall_markers selects facets of wall_facets, which is not given");
  s.check([](const partrac::Params& p){ return !p.has("wall_facets") || p.has("wall_markers"); },
          "wall_facets needs wall_markers, the markers of its no-slip walls");
  s.check([](const partrac::Params& p){
            return !p.has("wall_markers")
              || p.get<std::string>("wall_markers").find_first_not_of(" \t") != std::string::npos; },
          "wall_markers is empty; give the markers of the no-slip walls, as 3,5,7:106");
  s.check([](const partrac::Params& p){ return !p.has("wall_facets") || p.get<bool>("wall_upgrade"); },
          "wall_facets names the walls of wall_upgrade, which is false");
}

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
  s.require<std::string>("velocity_space", mode == "fenics" ? "velocity element, P1 to P4"
                                                            : "velocity element, P1 or P2");
  s.require<std::string>("pressure_space", "pressure element, as P1 or P2");
  s.optional<std::string>("phase_field_space", "phase field element, as P1 or P2; checked against the field "
                                               "file or the mesh cache where the phase field is read");
  if (mode == "fenics"){
    // Read by DolfInterpol only
    s.opt<std::string>("renumber_cells", "auto", "renumber cells by their dofs: auto (if poorly ordered), never, always");
    s.choices("renumber_cells", {"auto", "never", "always"});
  }
  else {
    // Read by SimplexInterpol and, where divfree is true, by SplitInterpol
    s.opt<bool>("mesh_cache", false, "keep the loaded tables beside the mesh and read them back");
    s.opt<bool>("include_phi", false, "read a phase field");
    s.opt<std::string>("phase_field", "phi", "phase field dataset in the field files");
    add_wall_upgrade(s, false, "raise the velocity's order next to walls at rest: a P1 velocity to P2 "
                               "(the edge rule), a P2 one by quartic terms and a correction in the cells with "
                               "a wall vertex at rest; the walls are wall_facets', else those at rest in the "
                               "first stamp with flow, which every later stamp must keep at rest; modes xdmftet, "
                               "xdmftriangle and openfoam default to true");
    add_wall_facets(s);
    s.opt<bool>("divfree", false, "evaluate a P2 velocity on each cell's barycentric split, where it is "
                                  "pointwise divergence-free; every cell's net flux must already be zero, "
                                  "which python/divfree/divfree_clean.py prepares; refuses wall_upgrade = true "
                                  "and mesh_cache = true beside it");
  }
  add_field_names(s);
  return s;
}

// A time series as frequency components: mode trianglefreq or tetfreq
inline partrac::Schema simplex_freq_schema(const std::string& mode){
  partrac::Schema s("dolfin_params.dat, mode=" + mode, "");
  add_mesh_params(s, mode != "trianglefreq");
  s.require<std::string>("freqstamps", "file of lines `t_k a_k file`, each component a_k cos(2 pi/tau (k t + t_k)) "
                                       "times its field with k the line number, or `omega_k phi_k a_k file`, "
                                       "a_k cos(omega_k t + phi_k)");
  s.require<std::string>("mesh", "mesh file");
  s.require<std::string>("velocity_space", "velocity element, as P1 or P2");
  s.require<std::string>("pressure_space", "pressure element, as P1 or P2");
  add_field_names(s);
  s.optional<double>("tau", "base period of `t_k a_k file` lines, 0 or less for none; not read with omega_k");
  add_wall_upgrade(s, false, "raise the velocity's order next to walls at rest in every component: a P1 "
                             "velocity to P2 (the edge rule), a P2 one by quartic terms and a correction in the "
                             "cells with a wall vertex at rest");
  add_wall_facets(s);
  s.require<double>("t_min", "start of the time interval");
  s.require<double>("t_max", "end of the time interval");
  s.check([](const partrac::Params& p){ return p.get<double>("t_max") >= p.get<double>("t_min"); },
          "t_max is before t_min");
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
  add_wall_upgrade(s, true, "velocity next to walls at rest: true (quadratic) or false (P1); the walls are "
                            "wall_facets', else those at rest in the first stamp with flow, which every later "
                            "stamp must keep at rest; the dolfin modes tet and triangle default to false");
  add_wall_facets(s);
  return s;
}

// An OpenFOAM case read in place: mode openfoam. The mesh, the times and the
// periodicity are the case's own.
inline partrac::Schema openfoam_schema(){
  partrac::Schema s("partrac_params.dat, mode=openfoam", "");
  add_mesh_params(s, true);
  s.opt<std::string>("velocity_field", "U", "velocity field in the time directories");
  s.opt<std::string>("pressure_field", "p", "pressure field in the time directories");
  s.opt<std::string>("split", "12", "tets a hex: 12, the fan from the cell centre (cellPoint's own), or 6, "
                                    "Dompierre's where the faces allow it; a 2D case 4 or 2 triangles a quad");
  s.choices("split", {"12", "6"});
  s.opt<std::string>("nodes", "least_squares", "values of U and p at the mesh points: least_squares (exact for "
                                               "linear fields) or inverse_distance (cellPoint's)");
  s.choices("nodes", {"least_squares", "inverse_distance"});
  add_wall_upgrade(s, true, "velocity next to walls at rest: true (quadratic) or false (P1); the walls are "
                            "those at rest in each stamp");
  s.optional<std::string>("wall_facets", "not read in mode openfoam");
  s.optional<std::string>("wall_markers", "not read in mode openfoam");
  s.check([](const partrac::Params& p){ return !p.has("wall_facets") && !p.has("wall_markers"); },
          "wall_facets and wall_markers are not read in mode openfoam: its walls are those at rest in the velocity");
  s.optional<std::string>("phase_field", "phase field in the time directories, a volScalarField such as "
                                          "alpha.water; none read without it");
  s.opt<bool>("include_phi", false, "not read: phase_field names the phase field");
  s.check([](const partrac::Params& p){ return !p.get<bool>("include_phi"); },
          "include_phi = true: mode openfoam reads the phase field phase_field names, as phase_field=alpha.water");
  s.check([](const partrac::Params& p){
            return !p.was_set("periodic_x") && !p.was_set("periodic_y") && !p.was_set("periodic_z"); },
          "periodic_x, periodic_y and periodic_z are the case's: its cyclic patches say which axes are periodic");
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
