#ifndef __TRACERS_TRIANGLERK4_SCHEMA_HPP
#define __TRACERS_TRIANGLERK4_SCHEMA_HPP

#include "typedefs.hpp"
#include "Params.hpp"
#include "experimental/initializer.hpp"

// Parameters accepted by this app
inline partrac::Schema tracers_triangleRK4_schema(){
  partrac::Schema s("tracers_triangleRK4");
  add_common_app_params(s);
  add_experimental_initializer_params(s);
  add_restart_params(s);
  s.require<int>("int_order", "integration order");
  return s;
}

#endif
