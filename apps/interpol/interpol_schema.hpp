#ifndef __INTERPOL_SCHEMA_HPP
#define __INTERPOL_SCHEMA_HPP

#include "typedefs.hpp"
#include "Params.hpp"

// Parameters accepted by this app
inline partrac::Schema interpol_schema(){
  partrac::Schema s("interpol");
  s.require<std::string>("mode", "interpolator type");
  s.require<Uint>("Nrw", "number of probe points");
  s.require<int>("int_order", "interpolation order");
  s.opt<double>("Dm", 0.0, "diffusivity, enters the folder name only");
  s.opt<double>("dt", 1.0, "timestep, enters the folder name only");
  s.opt<double>("t0", 0.0, "time to probe the field at");
  s.opt<double>("U", 1.0, "velocity scale");
  s.opt<int>("seed", 0, "random seed");
  s.opt<int>("num_threads", 0, "OpenMP threads, 0 = leave alone");
  s.opt<bool>("random", true, "draw the seed randomly");
  s.opt<std::string>("tag", "", "appended to the folder name");
  s.opt<std::string>("restart_folder", "", "folder to restart from");
  s.runtime<std::string>("folder", "", "output folder");
  s.choices("mode", {"analytic", "structured", "lbm", "felbm", "fenics",
                     "tet", "triangle", "trianglefreq", "tetfreq", "xdmftriangle", "xdmftet", "openfoam"});
  return s;
}

#endif
