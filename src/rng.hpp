#ifndef __RNG_HPP
#define __RNG_HPP

#include <random>
#include <vector>
#include <omp.h>

#include "Params.hpp"

// One generator per thread. seed_seq mixes the pair, so the streams are independent.
inline std::vector<std::mt19937> make_generators(const partrac::Params& prm){
  std::vector<std::mt19937> gens;
  const int N = omp_get_max_threads();
  gens.reserve(N);
  for (int i = 0; i < N; ++i){
    std::mt19937 gen;
    if (prm.get<bool>("random")){
      std::random_device rd;
      gen.seed(rd());
    }
    else {
      std::seed_seq rd{prm.get<int>("seed"), i};
      gen.seed(rd);
    }
    gens.emplace_back(gen);
  }
  return gens;
}

#endif
