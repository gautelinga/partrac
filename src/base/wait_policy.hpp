#ifndef __WAIT_POLICY_HPP
#define __WAIT_POLICY_HPP

// Passive OpenMP waiting unless OMP_WAIT_POLICY or GOMP_SPINCOUNT is set;
// libgomp reads them before main, hence the re-exec (set one under valgrind)

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <sched.h>
#include <unistd.h>
#include <vector>

namespace partrac {

// The places' union as the affinity mask: libgomp binds the main thread to
// its first place, and an exec keeps that mask
inline void bind_to_all_places(){
  const int n = omp_get_num_places();
  if (n <= 0)
    return;
  std::vector<int> ids;
  for (int p = 0; p < n; ++p){
    std::vector<int> place(omp_get_place_num_procs(p));
    omp_get_place_proc_ids(p, place.data());
    ids.insert(ids.end(), place.begin(), place.end());
  }
  int top = -1;
  for (const int id : ids) top = std::max(top, id);
  if (top < 0)
    return;
  cpu_set_t* set = CPU_ALLOC(top + 1);
  const std::size_t size = CPU_ALLOC_SIZE(top + 1);
  CPU_ZERO_S(size, set);
  for (const int id : ids)
    if (id >= 0)
      CPU_SET_S(id, size, set);
  if (sched_setaffinity(0, size, set) != 0)
    std::cerr << "Warning: the OpenMP places could not be restored as the affinity mask, "
              << "so the threads may share one place; set OMP_WAIT_POLICY to skip the restart" << std::endl;
  CPU_FREE(set);
}

inline void wait_passively(char* argv[]){
  if (std::getenv("OMP_WAIT_POLICY") || std::getenv("GOMP_SPINCOUNT")) return;
  setenv("OMP_WAIT_POLICY", "passive", 1);
  bind_to_all_places();
  execv("/proc/self/exe", argv);   // returns only if it failed: run as is
}

}  // namespace partrac

#endif
