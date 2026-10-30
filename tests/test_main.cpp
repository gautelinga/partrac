#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>

#include "wait_policy.hpp"

int main(int argc, char* argv[]){
  partrac::wait_passively(argv);
  return Catch::Session().run(argc, argv);
}
