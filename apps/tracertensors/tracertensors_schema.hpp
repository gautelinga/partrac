#ifndef __TRACERTENSORS_SCHEMA_HPP
#define __TRACERTENSORS_SCHEMA_HPP

#include "TracerParams.hpp"

// Parameters accepted by this app
inline partrac::Schema tracertensors_schema(){
  partrac::Schema s("tracertensors");
  add_tracer_params(s, {"reinject", 2, false, false, false});
  return s;
}

#endif
