#ifndef __TRACERVECTORS_SCHEMA_HPP
#define __TRACERVECTORS_SCHEMA_HPP

#include "TracerParams.hpp"

// Parameters accepted by this app
inline partrac::Schema tracervectors_schema(){
  partrac::Schema s("tracervectors");
  add_tracer_params(s, {"mark", 2, true, true, false});
  return s;
}

#endif
