#include <catch2/catch.hpp>

#include <cmath>
#include <memory>
#include <vector>

#include "typedefs.hpp"
#include "Error.hpp"
#include "Interpol.hpp"
#include "ParticleSet.hpp"
#include "mesh.hpp"

// This file exists partly for what it tests and partly for the fact that it
// compiles at all: the core headers define functions in the header, and this
// second translation unit must link.
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

// ---------------------------------------------------------------------------
// Remeshing on meshes built in memory. The app tests reach mesh.cpp only
// through whole runs, where a failure says that something moved, not what.

namespace {

// Fluid everywhere, at rest: enough for the inside check of an insertion
struct EverywhereInterpol final : public Interpol {
  EverywhereInterpol() : Interpol("") {}
  double get_t_min() override { return 0.; }
  double get_t_max() override { return 1.; }
  void update(const double) override {}
  bool locate(const Vector3d&, const double, CellPos&) override { return true; }
  void evaluate(const Vector3d&, const double, const CellPos&, PointValues&) override {}
  using Interpol::locate;
  using Interpol::evaluate;
};

double total_length(const EdgesType& edges, const ParticleSet& ps){
  double s = 0.;
  for (const auto& e : edges) s += ps.dist(e.first[0], e.first[1]);
  return s;
}

}  // namespace

TEST_CASE("refining a strip splits every edge longer than ds_max at its midpoint", "[mesh]") {
  ParticleSet ps(std::make_shared<EverywhereInterpol>(), 64);
  ps.add({{0., 0., 0.}, {1., 0., 0.}, {2., 0., 0.}}, 0);
  FacesType faces;
  EdgesType edges = {{{0, 1}, 1.0}, {{1, 2}, 1.0}};
  Edge2FacesType edge2faces;
  Node2EdgesType node2edges;
  compute_edge2faces(edge2faces, faces, edges);
  compute_node2edges(node2edges, edges, ps.N());
  EdgesListType edges_inlet;
  NodesListType nodes_inlet;
  std::vector<Vector3d> pos_inj;
  EdgesType edges_inj;

  const Uint n_add = refinement(faces, edges, edge2faces, node2edges, edges_inlet, nodes_inlet,
                                pos_inj, edges_inj, ps, 0.6, 0., false);
  REQUIRE( n_add == 2 );
  REQUIRE( ps.N() == 5 );
  REQUIRE( edges.size() == 4 );
  for (const auto& e : edges)
    REQUIRE( ps.dist(e.first[0], e.first[1]) == Approx(0.5) );
  REQUIRE( total_length(edges, ps) == Approx(2.0) );
  // Every node of a strip has one or two edges, and the table says so
  std::size_t ends = 0;
  for (Uint i = 0; i < ps.N(); ++i){
    REQUIRE( (node2edges[i].size() == 1 || node2edges[i].size() == 2) );
    if (node2edges[i].size() == 1) ++ends;
  }
  REQUIRE( ends == 2 );
  // Short enough already: nothing to do
  REQUIRE( refinement(faces, edges, edge2faces, node2edges, edges_inlet, nodes_inlet,
                      pos_inj, edges_inj, ps, 0.6, 0., false) == 0 );

  SECTION("and coarsening a straight strip keeps its length"){
    const Uint n_rem = coarsening(faces, edges, edge2faces, node2edges, edges_inlet, nodes_inlet,
                                  ps, 0.75, 0.);
    REQUIRE( n_rem > 0 );
    REQUIRE( total_length(edges, ps) == Approx(2.0) );
  }
}

TEST_CASE("refining a sheet keeps its area and stays a disk", "[mesh]") {
  // The unit square as two triangles across the diagonal 0-2
  ParticleSet ps(nullptr, 64);
  ps.add({{0., 0., 0.}, {1., 0., 0.}, {1., 1., 0.}, {0., 1., 0.}}, 0);
  EdgesType edges = {{{0, 1}, 1.0}, {{1, 2}, 1.0}, {{0, 2}, std::sqrt(2.)}, {{2, 3}, 1.0}, {{0, 3}, 1.0}};
  FacesType faces = {{{0, 1, 2}, 0.5}, {{2, 3, 4}, 0.5}};
  Edge2FacesType edge2faces;
  Node2EdgesType node2edges;
  compute_edge2faces(edge2faces, faces, edges);
  compute_node2edges(node2edges, edges, ps.N());
  EdgesListType edges_inlet;
  NodesListType nodes_inlet;
  std::vector<Vector3d> pos_inj;
  EdgesType edges_inj;

  // Only the diagonal is longer than 1.2
  const Uint n_add = sheet_refinement(faces, edges, edge2faces, node2edges, edges_inlet, nodes_inlet,
                                      pos_inj, edges_inj, ps, 1.2, 0., false, false);
  REQUIRE( n_add == 1 );
  REQUIRE( ps.N() == 5 );
  REQUIRE( faces.size() == 4 );
  REQUIRE( edges.size() == 8 );
  double area = 0.;
  for (Uint f = 0; f < faces.size(); ++f) area += ps.triangle_area(f, faces, edges);
  REQUIRE( area == Approx(1.0) );
  // Euler's formula for a disk
  REQUIRE( int(ps.N()) - int(edges.size()) + int(faces.size()) == 1 );
  // The new node is the diagonal's midpoint
  REQUIRE( (ps.x(4) - Vector3d(0.5, 0.5, 0.)).norm() < 1e-14 );
}

TEST_CASE("a broken mesh invariant throws an internal error", "[mesh][errors]") {
  EdgesType edges = {{{0, 1}, 1.0}, {{2, 3}, 1.0}};
  REQUIRE_THROWS_WITH( get_common_entry(0, 1, edges), Catch::Contains("internal:") );
}

TEST_CASE("a field a particle set does not have is an error, by name", "[mesh][errors]") {
  ParticleSet ps(nullptr, 4);
  ps.add({{0., 0., 0.}}, 0);
  for (const char* what : {"scalar", "vector", "tensor"}){
    INFO(what);
    const std::string w = what;
    const auto dump = [&]{
      if (w == "scalar") ps.dump_scalar("unused.h5", "no_such_field");
      if (w == "vector") ps.dump_vector("unused.h5", "no_such_field");
      if (w == "tensor") ps.dump_tensor("unused.h5", "no_such_field");
    };
    const auto load = [&]{
      if (w == "scalar") ps.load_scalar("unused.h5", "no_such_field");
      if (w == "vector") ps.load_vector("unused.h5", "no_such_field");
      if (w == "tensor") ps.load_tensor("unused.h5", "no_such_field");
    };
    REQUIRE_THROWS_WITH( dump(), Catch::Contains("no field 'no_such_field'") );
    REQUIRE_THROWS_WITH( load(), Catch::Contains("no field 'no_such_field'") );
  }
  REQUIRE_THROWS_AS( ps.dump_as("no_such_field", "name"), partrac::Error );
}
