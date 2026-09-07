#ifndef __OMP_TEST_SCHEMA_HPP
#define __OMP_TEST_SCHEMA_HPP

#include "typedefs.hpp"
#include "Params.hpp"
#include "experimental/initializer.hpp"

// Parameters accepted by this app
inline partrac::Schema omp_test_schema(){
  partrac::Schema s("omp_test");
  add_common_app_params(s);
  add_experimental_initializer_params(s);
  s.require<int>("int_order", "integration order");
  s.opt<int>("num_threads", 0, "OpenMP threads, 0 = leave alone");
  return s;
}

#endif
