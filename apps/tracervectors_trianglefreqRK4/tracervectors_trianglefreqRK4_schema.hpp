#ifndef __TRACERVECTORS_TRIANGLEFREQRK4_SCHEMA_HPP
#define __TRACERVECTORS_TRIANGLEFREQRK4_SCHEMA_HPP

#include "typedefs.hpp"
#include "Params.hpp"
#include "experimental/initializer.hpp"

// Parameters accepted by this app
inline partrac::Schema tracervectors_trianglefreqRK4_schema(){
  partrac::Schema s("tracervectors_trianglefreqRK4");
  add_common_app_params(s);
  add_experimental_initializer_params(s);
  add_restart_params(s);
  s.opt<int>("num_threads", 0, "OpenMP threads, 0 = leave alone");
  return s;
}

#endif
