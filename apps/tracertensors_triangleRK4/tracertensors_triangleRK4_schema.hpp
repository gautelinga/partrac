#ifndef __TRACERTENSORS_TRIANGLERK4_SCHEMA_HPP
#define __TRACERTENSORS_TRIANGLERK4_SCHEMA_HPP

#include "typedefs.hpp"
#include "Params.hpp"
#include "experimental/initializer.hpp"

// Parameters accepted by this app
inline partrac::Schema tracertensors_triangleRK4_schema(){
  partrac::Schema s("tracertensors_triangleRK4");
  add_common_app_params(s);
  add_experimental_initializer_params(s);
  s.token_choices("init_mode", "_", {"points"});
  add_restart_params(s);
  s.opt<int>("num_threads", 0, "OpenMP threads, 0 = leave alone");
  s.opt<int>("sort_every", 0, "reorder particles by cell every this many steps, 0 = never");
  return s;
}

#endif
