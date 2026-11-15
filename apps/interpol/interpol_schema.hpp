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
  s.opt<std::string>("points", "", "HDF5 file of points to probe in place of random ones: a row a point in its order, with inside 1 or 0 and the ids beside them if any");
  s.opt<std::string>("points_dataset", "points", "N x 3 positions in points, as a dump's <t>/points; probed at t0, not at the dump's time");
  s.require_if<Uint>("Nrw", 0,
                     [](const partrac::Params& p){ return p.get<std::string>("points").empty(); },
                     "points is not given", "number of probe points; with points, must be the file's count");
  add_int_order(s, "interpolation order");
  s.opt<double>("Dm", 0.0, "diffusivity, enters the folder name only");
  s.opt<double>("dt", 1.0, "timestep, enters the folder name only");
  s.opt<double>("t0", 0.0, "time to probe the field at");
  add_app_params(s);
  s.choices("mode", interpol_modes());
  return s;
}

#endif
