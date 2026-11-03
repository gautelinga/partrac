// The step counts of the dump, statistics and checkpoint intervals: an
// interval that is a whole number of timesteps is that number of steps, even
// where intv/dt falls just short of it in floating point, as 0.3/0.1 does;
// a fractional ratio is cut down; anything up to one step is one step.
#include <catch2/catch.hpp>

#include <limits>

#include "intervals.hpp"

TEST_CASE("An interval of whole timesteps is that many steps", "[intervals]"){
  CHECK(0.3/0.1 < 3.);                    // the round-off this is about
  CHECK(steps_per(0.3, 0.1) == 3);
  CHECK(steps_per(0.7, 0.1) == 7);
  CHECK(steps_per(0.6, 0.2) == 3);
  CHECK(steps_per(1.0, 0.1) == 10);
  CHECK(steps_per(0.125, 0.0625) == 2);
  CHECK(steps_per(3e-3, 1e-3) == 3);
  for (int k = 1; k <= 1000; ++k)
    CHECK(steps_per(k*0.1, 0.1) == Uint(k));
}

TEST_CASE("A fractional interval is cut down to whole steps", "[intervals]"){
  CHECK(steps_per(0.25, 0.1) == 2);
  CHECK(steps_per(0.39, 0.1) == 3);
  CHECK(steps_per(1.0 - 1e-6, 0.1) == 9);
}

TEST_CASE("Up to one timestep, or no number, is one step", "[intervals]"){
  CHECK(steps_per(0.1, 0.1) == 1);
  CHECK(steps_per(0.05, 0.1) == 1);
  CHECK(steps_per(0., 0.1) == 1);
  CHECK(steps_per(std::numeric_limits<double>::quiet_NaN(), 0.1) == 1);
  CHECK(steps_per(1e300, 1e-300) == std::numeric_limits<Uint>::max());
}

TEST_CASE("at_interval dumps every whole interval", "[intervals]"){
  CHECK_FALSE(at_interval(2, 0.3, 0.1));
  CHECK(at_interval(3, 0.3, 0.1));
  CHECK(at_interval(6, 0.3, 0.1));
  CHECK_FALSE(at_interval(3, 0., 0.1));
}
