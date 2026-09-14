#ifndef __TRACERS_SCHEMA_HPP
#define __TRACERS_SCHEMA_HPP

#include "TracerParams.hpp"

// Parameters accepted by this app
inline partrac::Schema tracers_schema(){
  partrac::Schema s("tracers");
  add_tracer_params(s, {"reinject", 0, false, false, false});
  return s;
}

#endif
