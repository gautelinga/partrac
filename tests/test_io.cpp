#include <catch2/catch.hpp>
#include "io.hpp"
#include "Interpol.hpp"
#include "Timestamps.hpp"

TEST_CASE("Timestamps span every key, sorted or not", "[timestamps]") {
  // t_max was tested only for a key that did not lower t_min, so the first key
  // never counted: one stamp, or a file in descending order, left it at -1e14
  auto span = [](std::vector<std::pair<double, std::string>> items){
    Timestamps ts;
    ts.initialize(items);
    return std::make_pair(ts.get_t_min(), ts.get_t_max());
  };
  REQUIRE(span({{0., "a"}}) == std::make_pair(0., 0.));
  REQUIRE(span({{0., "a"}, {10., "b"}}) == std::make_pair(0., 10.));
  REQUIRE(span({{10., "b"}, {0., "a"}}) == std::make_pair(0., 10.));
  REQUIRE(span({{20., "c"}, {0., "a"}, {10., "b"}}) == std::make_pair(0., 20.));

  MultiTimestamps mts;
  mts.initialize({{5., {"u.h5", "u"}}});
  REQUIRE(mts.get_t_min() == 5.);
  REQUIRE(mts.get_t_max() == 5.);
}

TEST_CASE("A single stamp brackets every time with itself", "[timestamps]") {
  std::vector<std::pair<double, std::string>> items = {{0., "up_0.h5"}};
  Timestamps ts;
  ts.initialize(items);
  StampPair sp = ts.get(1.);
  REQUIRE(sp.prev.t == 0.);
  REQUIRE(sp.next.t == 0.);
  REQUIRE(sp.prev.filename == sp.next.filename);
}

TEST_CASE("A degenerate stamp bracket holds the field", "[timestamps]") {
  // the plain expressions, where the bracket has width
  REQUIRE(stamp_weight(0.25, 0., 1.) == 0.25);
  REQUIRE(stamp_rate(3., 1., 0., 2.) == 1.);
  // and no 0/0 where it has none
  REQUIRE(stamp_weight(2., 2., 2.) == 0.);
  REQUIRE(stamp_rate(3., 3., 2., 2.) == 0.);
  REQUIRE(stamp_rate(Vector3d(1., 2., 3.), Vector3d(1., 2., 3.), 2., 2.).isZero(0.));
  REQUIRE(stamp_rate(Matrix3d::Identity().eval(), Matrix3d::Identity().eval(), 2., 2.).isZero(0.));
}

TEST_CASE("Passing case", "[pass]") {
  REQUIRE ( 1 < 2 );
}

TEST_CASE("print_param", "[print_param]") {
  // not implemented
  print_param("test", 2.0);

  REQUIRE ( 1 == 1 );
}
