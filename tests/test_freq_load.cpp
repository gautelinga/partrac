// The dolfin-free loader of SimplexFreqInterpol, on fixtures dolfin writes
// here the way data_example/sine_trianglefreq_p2/generate_up.py does: one
// steady Taylor-Hood pair per frequency component, written from a constrained
// space, and a freqstamps file giving each component its shift and amplitude.
// The components share one function space, so only the first is read with its
// dof table and the rest go through the mapping it built -- which is where a
// periodic-reduced space bites: a stored dof serves a node and its images, and
// a later component that reaches only one of them leaves the first component's
// values on a fifth of the cells. The tolerances here are the interpolation
// error; that failure is of order one. Triangles and tets are both covered,
// since the two modes share every table and the evaluation.
#ifdef USE_DOLFIN
#include <catch2/catch.hpp>

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <dolfin.h>

#include "Error.hpp"
#include "h5direct.hpp"
#include "TriangleFreqInterpol.hpp"
#include "TetFreqInterpol.hpp"
#include "dolfin_ref.hpp"
#include "taylor_hood.hpp"

namespace {

// A case directory of its own per test, removed with its contents
struct CaseDir {
  std::filesystem::path path;
  explicit CaseDir(const std::string& tag)
    : path(std::filesystem::temp_directory_path() / ("partrac_freq_" + tag)) {
    std::filesystem::remove_all(path);
    std::filesystem::create_directories(path);
  }
  ~CaseDir(){ std::filesystem::remove_all(path); }
  std::string params() const { return (path / "dolfin_params.dat").string(); }
};

// The amplitude and time shift of each component, as freqstamps.dat gives them
constexpr int n_freq = 2;
constexpr double amp[n_freq] = {1.0, 0.5};
constexpr double shift[n_freq] = {0.0, 0.125};
constexpr double tau = 1.0;

// Component k of the velocity, of period one along every axis
template<int D>
Vector3d comp_u(const Vector3d& x, const int k){
  if constexpr (D == 2)
    return Vector3d(std::sin(2.*M_PI*(x[1] + 0.1*double(k))),
                    std::sin(2.*M_PI*(x[0] + 0.2*double(k))), 0.);
  else
    return Vector3d(std::sin(2.*M_PI*(x[1] + 0.1*double(k))),
                    std::sin(2.*M_PI*(x[2] + 0.2*double(k))),
                    std::sin(2.*M_PI*(x[0] + 0.3*double(k))));
}

double comp_p(const Vector3d& x, const int k){
  return std::cos(2.*M_PI*(x[0] + 0.3*double(k)));
}

// The field the loader has to give back: the components summed with a cosine
// of the base frequency, as evaluate weights them
template<int D>
Vector3d want_u(const Vector3d& x, const double t){
  Vector3d u = Vector3d::Zero();
  for (int k = 0; k < n_freq; ++k)
    u += amp[k]*std::cos(2.*M_PI/tau*(double(k)*t + shift[k]))*comp_u<D>(x, k);
  return u;
}

double want_p(const Vector3d& x, const double t){
  double p = 0.;
  for (int k = 0; k < n_freq; ++k)
    p += amp[k]*std::cos(2.*M_PI/tau*(double(k)*t + shift[k]))*comp_p(x, k);
  return p;
}

template<int D>
class UExpr : public dolfin::Expression {
public:
  explicit UExpr(const int k) : dolfin::Expression(D), k_(k) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    const Vector3d u = comp_u<D>(Vector3d(x[0], x[1], D == 3 ? x[2] : 0.), k_);
    for (int d = 0; d < D; ++d) v[d] = u[d];
  }
private:
  int k_;
};

template<int D>
class PExpr : public dolfin::Expression {
public:
  explicit PExpr(const int k) : k_(k) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    v[0] = comp_p(Vector3d(x[0], x[1], D == 3 ? x[2] : 0.), k_);
  }
private:
  int k_;
};

// A periodic case with one file per component, written from a constrained
// space so that a vertex and its image share a stored dof. u_spaces names the
// element each component's velocity is written in, so a case whose components
// disagree can be written too.
template<typename Cell>
void write_case(const CaseDir& c, const std::size_t n, const std::vector<std::string>& u_spaces){
  constexpr int D = Cell::n_verts - 1;
  std::shared_ptr<dolfin::Mesh> mesh;
  if constexpr (D == 2) mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitSquareMesh(n, n));
  else                  mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitCubeMesh(n, n, n));
  const std::vector<bool> per = {true, true, D == 3};
  const Vector3d x_max(1., 1., D == 3 ? 1. : 0.);
  const std::shared_ptr<const dolfin::SubDomain> pbc =
    std::make_shared<PeriodicBC>(per, Vector3d::Zero(), x_max, D);
  {
    dolfin::HDF5File f(MPI_COMM_WORLD, (c.path / "mesh.h5").string(), "w");
    f.write(*mesh, "mesh");
  }
  std::ofstream stamps(c.path / "freqstamps.dat");
  for (int k = 0; k < int(u_spaces.size()); ++k){
    std::shared_ptr<dolfin::FunctionSpace> V, P;
    Uint ncoeffs_u = 0, ncoeffs_p = 0;
    taylor_hood_spaces<Cell>(u_spaces[std::size_t(k)], "P1", true, mesh, pbc, V, P,
                             ncoeffs_u, ncoeffs_p);
    dolfin::Function u(V), p(P);
    const UExpr<D> ue(k);
    const PExpr<D> pe(k);
    u.interpolate(ue);
    p.interpolate(pe);
    const std::string name = "up_" + std::to_string(k) + ".h5";
    dolfin::HDF5File f(MPI_COMM_WORLD, (c.path / name).string(), "w");
    f.write(u, "u");
    f.write(p, "p");
    stamps << shift[k] << " " << amp[k] << " " << name << "\n";
  }
  stamps.close();
  std::ofstream(c.params())
    << "velocity_space=" << u_spaces[0] << "\npressure_space=P1\nfreqstamps=freqstamps.dat\n"
    << "mesh=mesh.h5\nperiodic_x=true\nperiodic_y=true\nperiodic_z="
    << (D == 3 ? "true" : "false") << "\n"
    << "tau=" << tau << "\nt_min=0\nt_max=1e8\n";
}

// The cell centroids of the case's mesh, and how many sit in a cell with a
// vertex on a max face: those cells read image nodes
template<typename Cell>
std::vector<Vector3d> centroids(const CaseDir& c, std::size_t& touching){
  constexpr int nv = Cell::n_verts;
  constexpr Uint gdim = Uint(nv - 1);
  const partrac::H5Id file = partrac::h5_open_read((c.path / "mesh.h5").string());
  std::vector<std::uint32_t> topo;
  std::vector<double> coords;
  partrac::h5_read(file, "mesh/topology", topo);
  partrac::h5_read(file, "mesh/coordinates", coords, gdim);
  const std::size_t ncells = topo.size()/nv;
  std::vector<Vector3d> out(ncells, Vector3d::Zero());
  touching = 0;
  for (std::size_t i = 0; i < ncells; ++i){
    bool on_face = false;
    for (int k = 0; k < nv; ++k){
      const double* v = coords.data() + std::size_t(topo[i*nv + k])*gdim;
      for (Uint d = 0; d < gdim; ++d){
        out[i][d] += v[d]/double(nv);
        if (v[d] > 1. - 1e-12) on_face = true;
      }
    }
    if (on_face) ++touching;
  }
  return out;
}

// Every component reaches every cell, at times where the two carry different
// weights, one of them (t = 0.25) where the second component alone moves the field
template<typename Interpol, typename Cell>
void check_components_reach_every_cell(const std::string& tag, const std::size_t n){
  CaseDir c(tag);
  write_case<Cell>(c, n, {"P2", "P2"});
  std::size_t touching = 0;
  const std::vector<Vector3d> pts = centroids<Cell>(c, touching);
  REQUIRE(touching > 0);
  constexpr int D = Cell::n_verts - 1;
  Interpol intp(c.params());
  for (const double t : {0., 0.25, 0.7}){
    intp.update(t);
    std::size_t inside = 0;
    for (const Vector3d& x : pts){
      CellPos pos;
      if (!intp.locate(x, t, pos)) continue;
      PointValues got(1.);
      intp.evaluate(x, t, pos, got);
      const Vector3d u = want_u<D>(x, t);
      for (int d = 0; d < D; ++d)
        REQUIRE(std::abs(got.get_u()[d] - u[d]) < 1e-2);
      REQUIRE(std::abs(got.get_p() - want_p(x, t)) < 6e-2);
      ++inside;
    }
    REQUIRE(inside == pts.size());
  }
}

// The velocity at the centroids and times a check samples
template<typename Interpol>
std::vector<double> sample_u(Interpol& intp, const std::vector<Vector3d>& pts, const int D){
  std::vector<double> out;
  for (const double t : {0., 0.25, 0.7}){
    intp.update(t);
    for (const Vector3d& x : pts){
      CellPos pos;
      REQUIRE(intp.locate(x, t, pos));
      PointValues got(1.);
      intp.evaluate(x, t, pos, got);
      for (int d = 0; d < D; ++d) out.push_back(got.get_u()[d]);
      out.push_back(got.get_a()[0]);
    }
  }
  return out;
}

void write_stamps(const CaseDir& c, const std::string& text){
  std::ofstream(c.path / "freqstamps.dat") << text;
}

// The case's parameter file without its base period, for `omega phi a file` lines
void drop_tau(const CaseDir& c){
  std::ifstream in(c.params());
  std::string text, line;
  while (std::getline(in, line))
    if (line.rfind("tau=", 0) != 0) text += line + "\n";
  in.close();
  std::ofstream(c.params()) << text;
}

}  // namespace

TEST_CASE("Every frequency component reaches the cells on a periodic face", "[freq_load]") {
  SECTION("triangle") {
    check_components_reach_every_cell<TriangleFreqInterpol, Triangle>("per_tri_p2", 12);
  }
  SECTION("tet") {
    check_components_reach_every_cell<TetFreqInterpol, Tet>("per_tet_p2", 10);
  }
}

TEST_CASE("A frequency component in another element is refused", "[freq_load]") {
  // The components are one function space written many times; a file holding
  // another element cannot be read through the first one's mapping
  SECTION("triangle") {
    CaseDir c("mixed_elements_tri");
    write_case<Triangle>(c, 4, {"P2", "P1"});
    REQUIRE_THROWS_AS(TriangleFreqInterpol(c.params()), partrac::Error);
  }
  SECTION("tet") {
    CaseDir c("mixed_elements_tet");
    write_case<Tet>(c, 3, {"P2", "P1"});
    REQUIRE_THROWS_AS(TetFreqInterpol(c.params()), partrac::Error);
  }
}

// What a step reads must be what evaluate gives: the frequency weights are
// kept per thread between calls, so a stale one would show here first
TEST_CASE("evaluate_motion gives evaluate's velocity, acceleration and gradient", "[freq_load]") {
  CaseDir c("motion_tet_p2");
  write_case<Tet>(c, 4, {"P2", "P2"});
  std::size_t touching = 0;
  const std::vector<Vector3d> pts = centroids<Tet>(c, touching);
  TetFreqInterpol intp(c.params());
  intp.set_int_order(2);   // the gradient is read too
  // the times of one RK4 step, asked in the order a step asks them
  for (const double t : {0.3, 0.35, 0.35, 0.4, 0.3}){
    intp.update(t);
    for (const Vector3d& x : pts){
      CellPos pos;
      REQUIRE(intp.locate(x, t, pos));
      PointValues all(1.), motion(1.);
      intp.evaluate(x, t, pos, all);
      intp.evaluate_motion(x, t, pos, motion);
      REQUIRE((all.U - motion.U).norm() == 0.);
      REQUIRE((all.A - motion.A).norm() == 0.);
      REQUIRE((all.gradU - motion.gradU).norm() == 0.);
      REQUIRE((all.gradA - motion.gradA).norm() == 0.);
      // the weights are a cosine of the time; a stale one would freeze them
      const Vector3d u = want_u<3>(x, t);
      for (int d = 0; d < 3; ++d)
        REQUIRE(std::abs(all.get_u()[d] - u[d]) < 4e-2);
    }
  }
}

TEST_CASE("A freqstamps line may give its angular frequency and phase", "[freq_load]") {
  // `omega phi a file` is a cos(omega t + phi); the case's components as
  // harmonics of tau = 1 read the same field that way, in any order
  CaseDir c("omega");
  write_case<Triangle>(c, 6, {"P2", "P2"});
  std::size_t touching = 0;
  const std::vector<Vector3d> pts = centroids<Triangle>(c, touching);
  std::vector<double> harmonic;
  { TriangleFreqInterpol intp(c.params()); harmonic = sample_u(intp, pts, 2); }
  drop_tau(c);
  for (const char* text : {"0 0 1 up_0.h5\n6.283185307179586 0.7853981633974483 0.5 up_1.h5\n",
                           "6.283185307179586 0.7853981633974483 0.5 up_1.h5\n\n0 0 1 up_0.h5\n"}){
    write_stamps(c, text);
    TriangleFreqInterpol intp(c.params());
    const std::vector<double> got = sample_u(intp, pts, 2);
    REQUIRE(got.size() == harmonic.size());
    for (std::size_t i = 0; i < got.size(); ++i)
      REQUIRE(std::abs(got[i] - harmonic[i]) < 1e-12);
  }
}

TEST_CASE("One frequency on two lines is a Fourier mode's cosine and sine parts", "[freq_load]") {
  // A travelling mode needs two fields a frequency: A cos(w t) + B sin(w t),
  // the sine part a cosine of phase -pi/2
  CaseDir c("fourier_pair");
  write_case<Triangle>(c, 12, {"P2", "P2"});
  drop_tau(c);
  const double w = 2.5;
  write_stamps(c, "2.5 0 1 up_0.h5\n2.5 -1.5707963267948966 1 up_1.h5\n");
  std::size_t touching = 0;
  const std::vector<Vector3d> pts = centroids<Triangle>(c, touching);
  TriangleFreqInterpol intp(c.params());
  for (const double t : {0., 0.1, 0.25, 0.7}){
    intp.update(t);
    for (const Vector3d& x : pts){
      CellPos pos;
      REQUIRE(intp.locate(x, t, pos));
      PointValues got(1.);
      intp.evaluate(x, t, pos, got);
      const Vector3d u = std::cos(w*t)*comp_u<2>(x, 0) + std::sin(w*t)*comp_u<2>(x, 1);
      const Vector3d a = w*(-std::sin(w*t)*comp_u<2>(x, 0) + std::cos(w*t)*comp_u<2>(x, 1));
      for (int d = 0; d < 2; ++d){
        REQUIRE(std::abs(got.get_u()[d] - u[d]) < 1e-2);
        REQUIRE(std::abs(got.get_a()[d] - a[d]) < 3e-2);
      }
    }
  }
}

TEST_CASE("Components all of frequency 0 are a steady field", "[freq_load]") {
  // Several means, each a_k times its field whatever the time; no acceleration
  CaseDir c("steady");
  write_case<Triangle>(c, 6, {"P2", "P2"});
  drop_tau(c);
  write_stamps(c, "0 0 0.707 up_0.h5\n 0 0 0.103 up_1.h5\n 0 0 0.405 up_0.h5\n");
  std::size_t touching = 0;
  const std::vector<Vector3d> pts = centroids<Triangle>(c, touching);
  TriangleFreqInterpol intp(c.params());
  intp.set_int_order(2);
  std::vector<Vector3d> first;
  for (const double t : {0., 0.25, 0.7, 13.1}){
    intp.update(t);
    for (std::size_t i = 0; i < pts.size(); ++i){
      const Vector3d& x = pts[i];
      CellPos pos;
      REQUIRE(intp.locate(x, t, pos));
      PointValues got(1.);
      intp.evaluate(x, t, pos, got);
      REQUIRE(got.get_a().norm() == 0.);
      REQUIRE(got.gradA.norm() == 0.);
      if (t == 0.){
        first.push_back(got.get_u());
        const Vector3d u = 1.112*comp_u<2>(x, 0) + 0.103*comp_u<2>(x, 1);
        for (int d = 0; d < 2; ++d) REQUIRE(std::abs(got.get_u()[d] - u[d]) < 1e-2);
      }
      else {
        REQUIRE(got.get_u() == first[i]);
      }
    }
  }
}

TEST_CASE("A freqstamps file may have comments and Windows line ends", "[freq_load]") {
  // A line whose first word starts with '#' is skipped, and a trailing CR is
  // not part of the file name
  CaseDir c("stamps_text");
  write_case<Triangle>(c, 6, {"P2", "P2"});
  std::size_t touching = 0;
  const std::vector<Vector3d> pts = centroids<Triangle>(c, touching);
  std::vector<double> plain, annotated;
  { TriangleFreqInterpol intp(c.params()); plain = sample_u(intp, pts, 2); }
  write_stamps(c, "# t_k a_k file\r\n0 1 up_0.h5\r\n  #the first harmonic\r\n0.125 0.5 up_1.h5\r\n");
  { TriangleFreqInterpol intp(c.params()); annotated = sample_u(intp, pts, 2); }
  REQUIRE(annotated == plain);
}

TEST_CASE("A freqstamps file it cannot read as one form is refused", "[freq_load]") {
  CaseDir c("bad_stamps");
  write_case<Triangle>(c, 3, {"P2", "P2"});
  const auto refused = [&](const std::string& text){
    write_stamps(c, text);
    REQUIRE_THROWS_AS(TriangleFreqInterpol(c.params()), partrac::Error);
  };
  SECTION("a phase on the mean"){ refused("0.1 1 up_0.h5\n0.125 0.5 up_1.h5\n"); }
  SECTION("the two forms mixed"){ refused("0.125 0.5 up_1.h5\n0 0 1 up_0.h5\n"); }
  SECTION("a number that is not one"){ refused("0 x up_0.h5\n"); }
  SECTION("a line of two columns"){ refused("0 up_0.h5\n"); }
  SECTION("tau beside omega"){ refused("0 0 1 up_0.h5\n"); }
  SECTION("omega without tau, a phase on the mean"){
    drop_tau(c);
    refused("0 0.3 1 up_0.h5\n");
  }
  SECTION("harmonics without tau"){
    drop_tau(c);
    refused("0 1 up_0.h5\n");
  }
}

#endif
