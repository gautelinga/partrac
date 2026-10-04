#include <catch2/catch.hpp>
#include <limits>
#include <sstream>
#include "Error.hpp"
#include "param_print.hpp"
#include "io.hpp"
#include "Interpol.hpp"
#include "Timestamps.hpp"

TEST_CASE("Timestamps span every key, sorted or not", "[timestamps]") {
  // t_max counts every key: one stamp, or a file in descending order
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

TEST_CASE("No stamps give no bracket", "[timestamps]") {
  Timestamps ts;
  REQUIRE_THROWS_AS(ts.get(0.), partrac::Error);
}

TEST_CASE("A time is bracketed by the stamps at and after it", "[timestamps]") {
  // stamps given out of order; before the first and from the last on, one stamp twice
  std::vector<std::pair<double, std::string>> items = {{2., "b"}, {0., "a"}, {5., "c"}};
  Timestamps ts;
  ts.initialize(items);
  auto bracket = [&](const double t){
    StampPair sp = ts.get(t);
    return std::make_pair(sp.prev.filename, sp.next.filename);
  };
  REQUIRE(bracket(-1.) == std::make_pair(std::string("a"), std::string("a")));
  REQUIRE(bracket(0.) == std::make_pair(std::string("a"), std::string("b")));
  REQUIRE(bracket(1.) == std::make_pair(std::string("a"), std::string("b")));
  REQUIRE(bracket(2.) == std::make_pair(std::string("b"), std::string("c")));
  REQUIRE(bracket(4.9) == std::make_pair(std::string("b"), std::string("c")));
  REQUIRE(bracket(5.) == std::make_pair(std::string("c"), std::string("c")));
  REQUIRE(bracket(7.) == std::make_pair(std::string("c"), std::string("c")));
  StampPair sp = ts.get(3.);
  REQUIRE(sp.prev.t == 2.);
  REQUIRE(sp.next.t == 5.);
}

TEST_CASE("The next stamp is the first strictly after a time", "[timestamps]") {
  // where the run loop cuts a step; none past the last
  std::vector<std::pair<double, std::string>> items = {{2., "b"}, {0., "a"}, {5., "c"}};
  Timestamps ts;
  ts.initialize(items);
  const double none = std::numeric_limits<double>::infinity();
  REQUIRE(ts.next_after(-1.) == 0.);
  REQUIRE(ts.next_after(0.) == 2.);
  REQUIRE(ts.next_after(1.) == 2.);
  REQUIRE(ts.next_after(2.) == 5.);
  REQUIRE(ts.next_after(5.) == none);
  REQUIRE(ts.next_after(7.) == none);

  MultiTimestamps mts;
  mts.initialize({{0., {"u.h5", "u0"}}, {2., {"u.h5", "u1"}}, {5., {"u.h5", "u2"}}});
  REQUIRE(mts.next_after(-1.) == 0.);
  REQUIRE(mts.next_after(0.) == 2.);
  REQUIRE(mts.next_after(4.9) == 5.);
  REQUIRE(mts.next_after(5.) == none);
}

TEST_CASE("The XDMF stamps bracket a time as the others do", "[timestamps]") {
  MultiTimestamps mts;
  mts.initialize({{0., {"u.h5", "u0"}}, {2., {"u.h5", "u1"}}, {5., {"u.h5", "u2"}}});
  auto bracket = [&](const double t){
    MultiStampPair sp = mts.get(t);
    return std::make_pair(sp.prev.it, sp.next.it);
  };
  REQUIRE(bracket(-1.) == std::make_pair(Uint(0), Uint(0)));
  REQUIRE(bracket(0.) == std::make_pair(Uint(0), Uint(1)));
  REQUIRE(bracket(1.) == std::make_pair(Uint(0), Uint(1)));
  REQUIRE(bracket(2.) == std::make_pair(Uint(1), Uint(2)));
  REQUIRE(bracket(5.) == std::make_pair(Uint(2), Uint(2)));
  REQUIRE(bracket(7.) == std::make_pair(Uint(2), Uint(2)));
  REQUIRE(mts.get(3.).prev.t == 2.);
  REQUIRE(mts.get(3.).next.t == 5.);

  // searched by bisection: out of order is refused
  MultiTimestamps unsorted;
  REQUIRE_THROWS_AS(unsorted.initialize({{2., {"u.h5", "u1"}}, {0., {"u.h5", "u0"}}}), partrac::Error);
}

TEST_CASE("An XDMF field with another number of steps than u is refused", "[timestamps]") {
  // p and phi are indexed by u's steps: more would read past u's keys, fewer past their own
  MultiTimestamps mts;
  mts.initialize({{0., {"u.h5", "u0"}}, {2., {"u.h5", "u1"}}});
  REQUIRE_THROWS_AS(mts.add("p", {{0., {"p.h5", "p0"}}, {2., {"p.h5", "p1"}}, {5., {"p.h5", "p2"}}}),
                    partrac::Error);
  REQUIRE_THROWS_AS(mts.add("p", {{0., {"p.h5", "p0"}}}), partrac::Error);
  mts.add("p", {{0., {"p.h5", "p0"}}, {2., {"p.h5", "p1"}}});
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

TEST_CASE("print_param writes a key and its value", "[print_param]") {
  std::ostringstream buf;
  std::streambuf* old = std::cout.rdbuf(buf.rdbuf());
  print_param("dt", 0.5);
  print_param("Nrw", 3);
  print_param("mode", std::string("tet"));
  std::cout.rdbuf(old);
  REQUIRE(buf.str() == "dt = 0.5\nNrw = 3\nmode = tet\n");
}

// ---------------------------------------------------------------------------
// The checkpoint fields and the HDF5 writers. A whole run writes and reads
// only what its element carries, so every writer is exercised here.

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <H5Cpp.h>
#include "ParticleSet.hpp"
#include "case_dir.hpp"

namespace {

// A file name of its own per case, removed when it goes out of scope
struct TempFile {
  std::string name;
  explicit TempFile(const std::string& tag) : name(temp_path(tag).string()) {}
  ~TempFile(){ std::remove(name.c_str()); }
};

// carry() is what allocates the element's arrays
ParticleSet three_particles(const TransportElement e = TransportElement::Point){
  ParticleSet ps(nullptr, 8);
  ps.add({{0., 0., 0.}, {1., 2., 3.}, {-1., 0.5, 7.}}, 0);
  ps.carry(e);
  return ps;
}

// A line element's w and a walker's generation are both dumped as w: no app
// asks for both, so the set refuses it whichever comes first
TEST_CASE("a line element and a walker generation cannot share a set", "[io]"){
  ParticleSet vec = three_particles(TransportElement::Vector);
  REQUIRE_THROWS_AS(vec.record_generation(), partrac::Error);
  ParticleSet walker = three_particles();
  walker.record_generation();
  REQUIRE_THROWS_AS(walker.carry(TransportElement::Vector), partrac::Error);
}

// One dataset of a written file, as doubles
std::vector<double> read_dataset(const std::string& file, const std::string& name, const std::size_t n){
  H5::H5File h5(file, H5F_ACC_RDONLY);
  H5::DataSet d = h5.openDataSet(name);
  std::vector<double> v(n);
  d.read(v.data(), H5::PredType::NATIVE_DOUBLE);
  return v;
}

}  // namespace

TEST_CASE("a checkpoint survives a write and a read, bit for bit", "[io]") {
  // Values a float would round and a careless format would truncate
  auto write_read = [](const ParticleSet& ps, ParticleSet& back, const bool with_t_loc){
    TempFile f("checkpoint.h5");
    {
      H5::H5File h5(f.name, H5F_ACC_TRUNC);
      ps.write_checkpoint(h5, with_t_loc);
    }
    H5::H5File h5(f.name, H5F_ACC_RDONLY);
    back.read_checkpoint(h5, with_t_loc);
  };
  SECTION("positions, ids and t_loc"){
    ParticleSet ps = three_particles();
    for (Uint i = 0; i < ps.N(); ++i){
      ps.set_x(i, Vector3d(1. / 3. + i, -1e-300, 1e8 / 7.));
      ps.set_t_loc(i, 0.5 * i + 1. / 3.);
    }
    ParticleSet back(nullptr, 8);
    write_read(ps, back, true);
    REQUIRE(back.N() == ps.N());
    for (Uint i = 0; i < back.N(); ++i){
      REQUIRE(back.x(i) == ps.x(i));
      REQUIRE(back.t_loc(i) == ps.t_loc(i));
      REQUIRE(back.id(i) == ps.id(i));
    }
  }
  SECTION("a vector element"){
    ParticleSet ps = three_particles(TransportElement::Vector);
    for (Uint i = 0; i < ps.N(); ++i){
      ps.set_rhohat(i, Vector3d(0.1 + i, 1. / 3., -1e-8 * (i + 1)).normalized());
      ps.set_w(i, 1. / (3. + i));
      ps.set_S(i, -2. / 7. * i);
    }
    ParticleSet back(nullptr, 8);
    back.carry(TransportElement::Vector);
    write_read(ps, back, false);
    for (Uint i = 0; i < back.N(); ++i){
      REQUIRE(back.rhohat(i) == ps.rhohat(i));
      REQUIRE(back.w(i) == ps.w(i));
      REQUIRE(back.S(i) == ps.S(i));
    }
  }
  SECTION("a tensor element, its frame as carried"){
    ParticleSet ps = three_particles(TransportElement::Tensor);
    for (Uint i = 0; i < ps.N(); ++i){
      Matrix3d F;
      F << 1. + i, 2.5, -3., 0.25, 1e-9, 7., -0.5, 1. / 7., 1e8;
      ps.set_F(i, F);
      Matrix3d M;
      M << 1.1, 0.2, 0., -0.3, 0.9, 0.1, 0., 0.05, 1.3;
      ps.advance_frame(i, M * ps.frame(i));
    }
    ParticleSet back(nullptr, 8);
    back.carry(TransportElement::Tensor);
    write_read(ps, back, false);
    for (Uint i = 0; i < back.N(); ++i){
      REQUIRE(back.frame(i) == ps.frame(i));
      REQUIRE(back.F(i) == ps.F(i));
    }
  }
  SECTION("a walker generation"){
    ParticleSet ps = three_particles();
    ps.record_generation();
    for (Uint i = 0; i < ps.N(); ++i) ps.set_generation(i, 3. * i);
    ParticleSet back(nullptr, 8);
    back.record_generation();
    write_read(ps, back, false);
    for (Uint i = 0; i < back.N(); ++i) REQUIRE(back.generation(i) == ps.generation(i));
  }
}

TEST_CASE("a checkpoint keeps the next id past removed particles; one without it numbers on from the largest id", "[io]") {
  TempFile f("checkpoint.h5");
  ParticleSet ps = three_particles();
  ps.set_N(2);   // id 2 removed
  {
    H5::H5File h5(f.name, H5F_ACC_TRUNC);
    ps.write_checkpoint(h5, false);
  }
  auto next_after_read = [&](){
    ParticleSet back(nullptr, 8);
    H5::H5File h5(f.name, H5F_ACC_RDONLY);
    back.read_checkpoint(h5, false);
    back.add({{0.5, 0.5, 0.5}}, back.N());
    return back.id(back.N() - 1);
  };
  REQUIRE(next_after_read() == 3);
  {
    H5::H5File h5(f.name, H5F_ACC_RDWR);
    h5.removeAttr("next_id");
  }
  REQUIRE(next_after_read() == 2);
}

TEST_CASE("a checkpoint without a field, or with one of the wrong shape, stops the run", "[io]") {
  TempFile f("checkpoint.h5");
  ParticleSet ps = three_particles();
  {
    H5::H5File h5(f.name, H5F_ACC_TRUNC);
    ps.write_checkpoint(h5, false);
    std::vector<double> short_t(2, 0.);
    scalar2hdf5(h5, "t_loc", short_t, 2);
  }
  H5::H5File h5(f.name, H5F_ACC_RDONLY);
  SECTION("missing"){
    ParticleSet back(nullptr, 8);
    back.carry(TransportElement::Vector);
    REQUIRE_THROWS_WITH(back.read_checkpoint(h5, false), Catch::Contains("no dataset 'rhohat'"));
  }
  SECTION("one row short"){
    ParticleSet back(nullptr, 8);
    REQUIRE_THROWS_WITH(back.read_checkpoint(h5, true), Catch::Contains("'t_loc' is 2 x 1, expected 3 x 1"));
  }
  SECTION("more particles than the set has room for"){
    ParticleSet back(nullptr, 2);
    REQUIRE_THROWS_WITH(back.read_checkpoint(h5, false), Catch::Contains("more than Nrw_max"));
  }
}

TEST_CASE("an old text checkpoint field is read to the digit", "[io]") {
  // As the text checkpoints were written: 17 significant digits, one particle a line
  auto write = [](const std::string& name, const std::vector<std::vector<double>>& rows){
    std::ofstream out(name);
    out << std::setprecision(17);
    for (const auto& r : rows){
      for (Uint k = 0; k < r.size(); ++k) out << r[k] << (k + 1 == r.size() ? "\n" : " ");
    }
  };
  SECTION("a vector field"){
    ParticleSet ps = three_particles(TransportElement::Vector);
    std::vector<std::vector<double>> rows;
    for (Uint i = 0; i < ps.N(); ++i){
      const Vector3d n = Vector3d(0.1 + i, 1. / 3., -1e-8 * (i + 1)).normalized();
      rows.push_back({n[0], n[1], n[2]});
    }
    TempFile f("rhohat.vec");
    write(f.name, rows);
    ps.load_vector(f.name, "rhohat");
    for (Uint i = 0; i < ps.N(); ++i)
      REQUIRE(ps.rhohat(i) == Vector3d(rows[i][0], rows[i][1], rows[i][2]));
  }
  SECTION("a tensor field held whole, factored on load"){
    ParticleSet ps = three_particles(TransportElement::Tensor);
    std::vector<Matrix3d> F0;
    std::vector<std::vector<double>> rows;
    for (Uint i = 0; i < ps.N(); ++i){
      Matrix3d F;
      F << 1. + i, 2.5, -3., 0.25, 1e-9, 7., -0.5, 1. / 7., 1e8;
      F0.push_back(F);
      rows.push_back({});
      for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) rows.back().push_back(F(r, c));
    }
    TempFile f("F.ten");
    write(f.name, rows);
    ps.load_tensor(f.name, "F");
    for (Uint i = 0; i < ps.N(); ++i)
      REQUIRE((ps.F(i) - F0[i]).norm() / F0[i].norm() < 1e-14);
  }
  SECTION("a scalar field"){
    ParticleSet ps = three_particles();
    std::vector<std::vector<double>> rows;
    for (Uint i = 0; i < ps.N(); ++i) rows.push_back({0.5 * i + 1. / 3.});
    TempFile f("t_loc.dat");
    write(f.name, rows);
    ps.load_scalar(f.name, "t_loc");
    for (Uint i = 0; i < ps.N(); ++i)
      REQUIRE(ps.t_loc(i) == rows[i][0]);
  }
}

TEST_CASE("the HDF5 writers lay a field out one particle after another", "[io]") {
  TempFile f("fields.h5");
  const Uint n = 3;
  std::vector<Vector3d> a = {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}};
  std::vector<Matrix3d> M(n);
  for (Uint i = 0; i < n; ++i)
    M[i] << 1. + 9 * i, 2. + 9 * i, 3. + 9 * i, 4. + 9 * i, 5. + 9 * i,
            6. + 9 * i, 7. + 9 * i, 8. + 9 * i, 9. + 9 * i;
  std::vector<double> s = {0.5, 1.5, 2.5};
  {
    H5::H5File h5(f.name, H5F_ACC_TRUNC);
    vector2hdf5(h5, "/v", a, n);
    tensor2hdf5(h5, "/t", M, n);
    scalar2hdf5(h5, "/s", s, n);
  }
  const std::vector<double> v = read_dataset(f.name, "/v", 3 * n);
  for (Uint i = 0; i < n; ++i)
    for (int k = 0; k < 3; ++k)
      REQUIRE(v[3 * i + k] == a[i][k]);
  const std::vector<double> t = read_dataset(f.name, "/t", 9 * n);
  for (Uint i = 0; i < n; ++i)
    for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 3; ++c)
        REQUIRE(t[9 * i + 3 * r + c] == M[i](r, c));
  const std::vector<double> sc = read_dataset(f.name, "/s", n);
  for (Uint i = 0; i < n; ++i) REQUIRE(sc[i] == s[i]);
}
