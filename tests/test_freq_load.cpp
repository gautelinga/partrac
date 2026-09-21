// The dolfin-free loader of TriangleFreqInterpol, on fixtures dolfin writes
// here the way data_example/sine_trianglefreq_p2/generate_up.py does: one
// steady Taylor-Hood pair per frequency component, written from a constrained
// space, and a freqstamps file giving each component its shift and amplitude.
// The components share one function space, so only the first is read with its
// dof table and the rest go through the mapping it built -- which is where a
// periodic-reduced space bites: a stored dof serves a node and its images, and
// a later component that reaches only one of them leaves the first component's
// values on a fifth of the cells. The tolerances here are the interpolation
// error; that failure is of order one.
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

// Component k of the velocity, of period one in x and y
Vector3d comp_u(const Vector3d& x, const int k){
  return Vector3d(std::sin(2.*M_PI*(x[1] + 0.1*double(k))),
                  std::sin(2.*M_PI*(x[0] + 0.2*double(k))), 0.);
}

double comp_p(const Vector3d& x, const int k){
  return std::cos(2.*M_PI*(x[0] + 0.3*double(k)));
}

// The field the loader has to give back: the components summed with a cosine
// of the base frequency, as evaluate weights them
Vector3d want_u(const Vector3d& x, const double t){
  Vector3d u = Vector3d::Zero();
  for (int k = 0; k < n_freq; ++k)
    u += amp[k]*std::cos(2.*M_PI/tau*(double(k)*t + shift[k]))*comp_u(x, k);
  return u;
}

double want_p(const Vector3d& x, const double t){
  double p = 0.;
  for (int k = 0; k < n_freq; ++k)
    p += amp[k]*std::cos(2.*M_PI/tau*(double(k)*t + shift[k]))*comp_p(x, k);
  return p;
}

class UExpr : public dolfin::Expression {
public:
  explicit UExpr(const int k) : dolfin::Expression(2), k_(k) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    const Vector3d u = comp_u(Vector3d(x[0], x[1], 0.), k_);
    v[0] = u[0];
    v[1] = u[1];
  }
private:
  int k_;
};

class PExpr : public dolfin::Expression {
public:
  explicit PExpr(const int k) : k_(k) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    v[0] = comp_p(Vector3d(x[0], x[1], 0.), k_);
  }
private:
  int k_;
};

// A periodic case with one file per component, written from a constrained
// space so that a vertex and its image share a stored dof. u_spaces names the
// element each component's velocity is written in, so a case whose components
// disagree can be written too.
void write_case(const CaseDir& c, const std::size_t n, const std::vector<std::string>& u_spaces){
  auto mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitSquareMesh(n, n));
  const std::vector<bool> per = {true, true, false};
  const std::shared_ptr<const dolfin::SubDomain> pbc =
    std::make_shared<PeriodicBC>(per, Vector3d::Zero(), Vector3d(1., 1., 0.), 2);
  {
    dolfin::HDF5File f(MPI_COMM_WORLD, (c.path / "mesh.h5").string(), "w");
    f.write(*mesh, "mesh");
  }
  std::ofstream stamps(c.path / "freqstamps.dat");
  for (int k = 0; k < int(u_spaces.size()); ++k){
    std::shared_ptr<dolfin::FunctionSpace> V, P;
    Uint ncoeffs_u = 0, ncoeffs_p = 0;
    taylor_hood_spaces<Triangle>(u_spaces[std::size_t(k)], "P1", true, mesh, pbc, V, P,
                                 ncoeffs_u, ncoeffs_p);
    dolfin::Function u(V), p(P);
    const UExpr ue(k);
    const PExpr pe(k);
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
    << "mesh=mesh.h5\nperiodic_x=true\nperiodic_y=true\nperiodic_z=false\n"
    << "tau=" << tau << "\nt_min=0\nt_max=1e8\n";
}

// The cell centroids of the case's mesh, and how many sit in a cell with a
// vertex on a max face: those cells read image nodes
std::vector<Vector3d> centroids(const CaseDir& c, std::size_t& touching){
  const partrac::H5Id file = partrac::h5_open_read((c.path / "mesh.h5").string());
  std::vector<std::uint32_t> topo;
  std::vector<double> coords;
  partrac::h5_read(file, "mesh/topology", topo);
  partrac::h5_read(file, "mesh/coordinates", coords, 2);
  const std::size_t ncells = topo.size()/3;
  std::vector<Vector3d> out(ncells, Vector3d::Zero());
  touching = 0;
  for (std::size_t i = 0; i < ncells; ++i){
    bool on_face = false;
    for (int k = 0; k < 3; ++k){
      const double* v = coords.data() + std::size_t(topo[i*3 + k])*2;
      for (Uint d = 0; d < 2; ++d){
        out[i][d] += v[d]/3.;
        if (v[d] > 1. - 1e-12) on_face = true;
      }
    }
    if (on_face) ++touching;
  }
  return out;
}

}  // namespace

TEST_CASE("Every frequency component reaches the cells on a periodic face", "[freq_load]") {
  CaseDir c("per_tri_p2");
  write_case(c, 12, {"P2", "P2"});
  std::size_t touching = 0;
  const std::vector<Vector3d> pts = centroids(c, touching);
  REQUIRE(touching > 0);
  TriangleFreqInterpol intp(c.params());
  // three times at which the two components carry different weights, one of
  // them (t = 0.25) where the second component alone moves the field
  for (const double t : {0., 0.25, 0.7}){
    intp.update(t);
    std::size_t inside = 0;
    for (const Vector3d& x : pts){
      CellPos pos;
      if (!intp.locate(x, t, pos)) continue;
      PointValues got(1.);
      intp.evaluate(x, t, pos, got);
      const Vector3d u = want_u(x, t);
      for (Uint d = 0; d < 2; ++d)
        REQUIRE(std::abs(got.get_u()[d] - u[d]) < 1e-2);
      REQUIRE(std::abs(got.get_p() - want_p(x, t)) < 6e-2);
      ++inside;
    }
    REQUIRE(inside == pts.size());
  }
}

TEST_CASE("A frequency component in another element is refused", "[freq_load]") {
  // The components are one function space written many times; a file holding
  // another element cannot be read through the first one's mapping
  CaseDir c("mixed_elements");
  write_case(c, 4, {"P2", "P1"});
  REQUIRE_THROWS_AS(TriangleFreqInterpol(c.params()), partrac::Error);
}

#endif
