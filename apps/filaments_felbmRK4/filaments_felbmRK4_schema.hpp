#ifndef __FILAMENTS_FELBMRK4_SCHEMA_HPP
#define __FILAMENTS_FELBMRK4_SCHEMA_HPP

#include "typedefs.hpp"
#include "Params.hpp"
#include "experimental/initializer.hpp"

// Parameters accepted by this app
inline partrac::Schema filaments_felbmRK4_schema(){
  partrac::Schema s("filaments_felbmRK4");
  add_common_app_params(s);
  add_experimental_initializer_params(s);
  add_restart_params(s);
  s.require<int>("int_order", "integration order");
  s.require<double>("ds_init", "initial edge length");
  s.opt<double>("resize_intv", 0.0, "resize interval");
  return s;
}

#endif
