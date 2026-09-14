#ifndef __TRACERVECTORS_SPATIAL_SCHEMA_HPP
#define __TRACERVECTORS_SPATIAL_SCHEMA_HPP

#include "TracerParams.hpp"

// Parameters accepted by this app
inline partrac::Schema tracervectors_spatial_schema(){
  partrac::Schema s("tracervectors_spatial");
  add_march_tracer_params(s, {"remove", 2, false, false, true});
  return s;
}

#endif
