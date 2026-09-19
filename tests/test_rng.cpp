#include <catch2/catch.hpp>

#include <random>
#include <set>
#include <vector>

// One generator per thread, seeded seed_seq{seed, i}. An ensemble over
// consecutive seeds is only independent if no stream is shared.

static std::vector<unsigned int> first_draws(int seed, int i, int n = 4) {
  std::seed_seq s{seed, i};
  std::mt19937 gen;
  gen.seed(s);
  std::vector<unsigned int> out;
  for (int k = 0; k < n; ++k) out.push_back(gen());
  return out;
}

TEST_CASE("every (seed, thread) pair gets its own stream", "[rng]") {
  std::set<std::vector<unsigned int>> seen;
  for (int seed = 1; seed <= 20; ++seed)
    for (int i = 0; i < 16; ++i)
      REQUIRE(seen.insert(first_draws(seed, i)).second);
}

TEST_CASE("consecutive seeds share no stream", "[rng]") {
  // thread i of seed k must not match thread i-1 of seed k+1
  for (int seed = 1; seed <= 20; ++seed)
    for (int i = 1; i < 16; ++i)
      REQUIRE(first_draws(seed, i) != first_draws(seed + 1, i - 1));
}
