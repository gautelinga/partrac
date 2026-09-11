#include <catch2/catch.hpp>

#include <cmath>
#include <vector>

#include "typedefs.hpp"
#include "ParticleSet.hpp"
#include "mesh.hpp"

// This file exists partly for what it tests and partly for the fact that it
// compiles at all: the core headers define their functions in the header, so
// until they were marked inline a second translation unit including them was a
// multiple-definition link error, and every app avoided it only by being one
// TU. If that regresses, this file stops linking.
//
// ParticleSet's constructor only stores the interpolator, so a null one is
// enough for anything that reads positions.

TEST_CASE("a triangle's area is half the cross product of two of its edges",
          "[mesh]") {
  ParticleSet ps(nullptr, 8);
  ps.add({{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.},   // a right triangle
          {2., 0., 0.}, {3., 0., 0.}, {4., 0., 0.}},  // three points on a line
         0);
  EdgesType edges;
  edges.push_back({{0, 1}, 1.0});
  edges.push_back({{1, 2}, 1.0});
  edges.push_back({{3, 4}, 1.0});
  edges.push_back({{4, 5}, 1.0});

  REQUIRE( ps.triangle_area(0, 1, edges) == Approx(0.5) );
  REQUIRE( ps.triangle_area(2, 3, edges) == 0.0 );   // collinear, exactly none
}

TEST_CASE("the distance between two nodes is symmetric and zero on itself",
          "[mesh]") {
  ParticleSet ps(nullptr, 4);
  ps.add({{0., 0., 0.}, {3., 4., 0.}}, 0);
  REQUIRE( ps.dist(0, 1) == Approx(5.0) );
  REQUIRE( ps.dist(1, 0) == Approx(5.0) );
  REQUIRE( ps.dist(0, 0) == 0.0 );
}

TEST_CASE("a node inserted between two lands on the midpoint", "[mesh]") {
  // check_if_inside is off: locating a point needs an interpolator, and the
  // placement is what is being checked
  ParticleSet ps(nullptr, 4);
  ps.add({{0., 0., 0.}, {1., 2., 4.}}, 0);
  REQUIRE( ps.N() == 2 );
  REQUIRE( ps.insert_node_between(0, 1, false) );
  REQUIRE( ps.N() == 3 );
  REQUIRE( ps.x(2)[0] == Approx(0.5) );
  REQUIRE( ps.x(2)[1] == Approx(1.0) );
  REQUIRE( ps.x(2)[2] == Approx(2.0) );
}

TEST_CASE("get_common_entry finds the node two edges share", "[mesh]") {
  EdgesType edges;
  edges.push_back({{4, 7}, 1.0});
  edges.push_back({{7, 9}, 1.0});
  REQUIRE( get_common_entry(0, 1, edges) == 7 );
  REQUIRE( get_common_entry(1, 0, edges) == 7 );
}
