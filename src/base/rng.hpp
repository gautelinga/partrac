#ifndef __RNG_HPP
#define __RNG_HPP

#include <iostream>
#include <random>
#include <vector>
#include <omp.h>

#include "Params.hpp"

// One generator per thread; random = true draws the seed and records it in prm
inline std::vector<std::mt19937> make_generators(partrac::Params& prm){
  if (prm.get<bool>("random")){
    std::random_device rd;
    const int seed = static_cast<int>(rd() >> 1);
    prm.set<int>("seed", seed);
    std::cout << "Seed " << seed << " drawn (random=true); random=false seed=" << seed
              << " repeats this run" << std::endl;
  }
  std::vector<std::mt19937> gens;
  const int N = omp_get_max_threads();
  gens.reserve(N);
  for (int i = 0; i < N; ++i){
    std::seed_seq rd{prm.get<int>("seed"), i};
    gens.emplace_back(rd);
  }
  return gens;
}

#endif
