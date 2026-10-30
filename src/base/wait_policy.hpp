#ifndef __WAIT_POLICY_HPP
#define __WAIT_POLICY_HPP

// Passive OpenMP waiting unless OMP_WAIT_POLICY or GOMP_SPINCOUNT is set;
// libgomp reads them before main, hence the re-exec (set one under valgrind)

#include <cstdlib>
#include <unistd.h>

namespace partrac {

inline void wait_passively(char* argv[]){
  if (std::getenv("OMP_WAIT_POLICY") || std::getenv("GOMP_SPINCOUNT")) return;
  setenv("OMP_WAIT_POLICY", "passive", 1);
  execv("/proc/self/exe", argv);   // returns only if it failed: run as is
}

}  // namespace partrac

#endif
