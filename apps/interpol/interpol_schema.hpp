#ifndef __INTERPOL_SCHEMA_HPP
#define __INTERPOL_SCHEMA_HPP

#include "typedefs.hpp"
#include "Params.hpp"
#include "interpol_factory.hpp"
#include "AppParams.hpp"

// Parameters accepted by this app
inline partrac::Schema interpol_schema(){
  partrac::Schema s("interpol");
  s.require<std::string>("mode", "interpolator type");
  s.require<Uint>("Nrw", "number of probe points");
  s.require<int>("int_order", "interpolation order");
  s.opt<double>("Dm", 0.0, "diffusivity, enters the folder name only");
  s.opt<double>("dt", 1.0, "timestep, enters the folder name only");
  s.opt<double>("t0", 0.0, "time to probe the field at");
  add_app_params(s);
  s.choices("mode", interpol_modes());
  return s;
}

#endif
