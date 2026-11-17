#ifndef __RNG_HPP
#define __RNG_HPP

#include <cstdint>
#include <cstring>
#include <iostream>
#include <random>
#include <vector>
#include <omp.h>

#include "Params.hpp"

// One generator per thread; a restart adds its time to the seed
inline std::vector<std::mt19937> make_generators(partrac::Params& prm){
  if (prm.get<bool>("random")){
    std::random_device rd;
    const int seed = static_cast<int>(rd() >> 1);
    prm.set<int>("seed", seed);
    std::cout << "Seed " << seed << " drawn (random=true); random=false seed=" << seed
              << " repeats this run" << std::endl;
  }
  std::vector<std::uint32_t> words{std::uint32_t(prm.get<int>("seed")), 0};
  if (prm.has("restart_folder") && prm.get<std::string>("restart_folder") != ""){
    const double t = prm.get<double>("t");
    std::uint64_t bits;
    std::memcpy(&bits, &t, sizeof bits);
    words.push_back(std::uint32_t(bits));
    words.push_back(std::uint32_t(bits >> 32));
  }
  std::vector<std::mt19937> gens;
  const int N = omp_get_max_threads();
  gens.reserve(N);
  for (int i = 0; i < N; ++i){
    words[1] = std::uint32_t(i);
    std::seed_seq rd(words.begin(), words.end());
    gens.emplace_back(rd);
  }
  return gens;
}

// A generator of the run's own for one sequence of draws (stream), as the
// strip filter's: the run's seed and words no thread's generator has
inline std::mt19937 stream_generator(const partrac::Params& prm, const std::uint32_t stream){
  const std::vector<std::uint32_t> words{std::uint32_t(prm.get<int>("seed")), stream, 0, 0, 0};
  std::seed_seq rd(words.begin(), words.end());
  return std::mt19937(rd);
}

#endif
