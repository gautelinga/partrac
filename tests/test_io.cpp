#include <catch2/catch.hpp>
#include "io.hpp"

TEST_CASE("Passing case", "[pass]") {
  REQUIRE ( 1 < 2 );
}

TEST_CASE("Failing case", "[fail]") {
  REQUIRE ( 1 > 2 );
}

TEST_CASE("print_param", "[print_param]") {
  // not implemented
  print_param("test", 2.0);

  REQUIRE ( 1 == 1 );
}
