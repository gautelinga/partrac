#ifndef __TRACERVECTORS_TRIANGLE_SPATIAL_SCHEMA_HPP
#define __TRACERVECTORS_TRIANGLE_SPATIAL_SCHEMA_HPP

#include "typedefs.hpp"
#include "Params.hpp"
#include "experimental/initializer.hpp"

// Parameters accepted by this app
inline partrac::Schema tracervectors_triangle_spatial_schema(){
  partrac::Schema s("tracervectors_triangle_spatial");
  add_common_app_params(s);
  add_experimental_initializer_params(s);
  s.token_choices("init_mode", "_", {"points"});
  add_restart_params(s);
  s.opt<int>("num_threads", 0, "OpenMP threads, 0 = leave alone");
  s.opt<int>("sort_every", 0, "reorder particles by cell every this many steps, 0 = never");
  s.require<double>("ds_max", "max edge length");
  s.opt<double>("Lt", 0.0, "tangential extent of the exit plane");
  s.opt<double>("u_eps", 1e-7, "velocity cutoff");
  return s;
}

#endif
