#include <catch2/catch.hpp>
#include "param_print.hpp"
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

// ---------------------------------------------------------------------------
// The checkpoint fields and the HDF5 writers. A whole run writes and reads
// these, but only for what its element carries, so the tensor and vector
// writers were never exercised by a test.

#include <cstdio>
#include <filesystem>
#include <H5Cpp.h>
#include "ParticleSet.hpp"

namespace {

// A file name of its own per case, removed when it goes out of scope
struct TempFile {
  std::string name;
  explicit TempFile(const std::string& tag)
    : name((std::filesystem::temp_directory_path() / ("partrac_test_" + tag)).string()) {}
  ~TempFile(){ std::remove(name.c_str()); }
};

// carry() is what allocates the element's arrays
ParticleSet three_particles(const TransportElement e = TransportElement::Point){
  ParticleSet ps(nullptr, 8);
  ps.add({{0., 0., 0.}, {1., 2., 3.}, {-1., 0.5, 7.}}, 0);
  ps.carry(e);
  return ps;
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

TEST_CASE("a checkpoint field survives a write and a read, to the digit", "[io]") {
  // Values a float would round and a careless format would truncate
  SECTION("a vector field"){
    ParticleSet ps = three_particles(TransportElement::Vector);
    std::vector<Vector3d> n0;
    for (Uint i = 0; i < ps.N(); ++i){
      ps.set_rhohat(i, Vector3d(0.1 + i, 1. / 3., -1e-8 * (i + 1)).normalized());
      n0.push_back(ps.rhohat(i));
    }
    TempFile f("rhohat.dat");
    ps.dump_vector(f.name, "rhohat");
    ParticleSet back = three_particles(TransportElement::Vector);
    back.load_vector(f.name, "rhohat");
    for (Uint i = 0; i < back.N(); ++i)
      REQUIRE((back.rhohat(i) - n0[i]).norm() < 1e-14);
  }
  SECTION("a tensor field"){
    ParticleSet ps = three_particles(TransportElement::Tensor);
    std::vector<Matrix3d> F0;
    for (Uint i = 0; i < ps.N(); ++i){
      Matrix3d F;
      F << 1. + i, 2.5, -3., 0.25, 1e-9, 7., -0.5, 1. / 7., 1e8;
      ps.set_F(i, F);
      F0.push_back(F);
    }
    TempFile f("F.dat");
    ps.dump_tensor(f.name, "F");
    ParticleSet back = three_particles(TransportElement::Tensor);
    back.load_tensor(f.name, "F");
    for (Uint i = 0; i < back.N(); ++i)
      REQUIRE((back.F(i) - F0[i]).norm() / F0[i].norm() < 1e-14);
  }
  SECTION("a scalar field"){
    ParticleSet ps = three_particles();
    TempFile f("t_loc.dat");
    for (Uint i = 0; i < ps.N(); ++i) ps.set_t_loc(i, 0.5 * i + 1. / 3.);
    ps.dump_scalar(f.name, "t_loc");
    ParticleSet back = three_particles();
    back.load_scalar(f.name, "t_loc");
    for (Uint i = 0; i < back.N(); ++i)
      REQUIRE(back.t_loc(i) == Approx(0.5 * i + 1. / 3.).epsilon(1e-14));
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
