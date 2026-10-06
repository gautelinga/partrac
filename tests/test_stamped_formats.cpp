// The two stamp formats behind one evaluation: a P1 field written by dolfin as
// a checkpoint and as XDMF has to read the same through either loader, with
// the near-wall rule, the phase field and the cell types, since only the
// reading of the files differs. Around that, what the checkpoint side gained
// with it: wall_p2 refused on a velocity that is P2 already, a phase field in
// an element of its own, and a mesh cache that carries the phase field.
#ifdef USE_DOLFIN
#include <catch2/catch.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include <dolfin.h>

#include "Error.hpp"
#include "Tet.hpp"
#include "Triangle.hpp"
#include "SimplexInterpol.hpp"
#include "XDMFInterpol.hpp"
#include "case_dir.hpp"
#include "dolfin_ref.hpp"
#include "taylor_hood.hpp"

namespace {

// Stamp k of a velocity at rest on the boundary of the unit box, so the
// near-wall rule applies at every wall; the factors keep it off every symmetry
class UExpr : public dolfin::Expression {
public:
  UExpr(const Uint gdim, const int k) : dolfin::Expression(gdim), gdim_(gdim), k_(k) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    double b = 1.;
    for (Uint d = 0; d < gdim_; ++d) b *= x[d]*(1. - x[d]);
    for (Uint d = 0; d < gdim_; ++d)
      v[d] = b*(1. + double(k_))*(d == 0 ? 4.*(1. + x[1]) : 1. + x[0] + 0.5*double(d));
  }
private:
  Uint gdim_;
  int k_;
};

class PExpr : public dolfin::Expression {
public:
  PExpr(const Uint gdim, const int k) : gdim_(gdim), k_(k) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    v[0] = std::cos(2.*x[0] + 0.2) - 0.5*x[1] + (gdim_ == 3 ? 0.3*x[2] : 0.) + double(k_);
  }
private:
  Uint gdim_;
  int k_;
};

// A quadratic, so a P2 space holds it exactly
double phi_exact(const Vector3d& x, const int k, const Uint gdim){
  return x[0]*x[0] - x[0]*x[1] + 0.5*x[1]*x[1] + (gdim == 3 ? x[2]*x[2] - 0.3*x[2] : 0.)
       + 0.5*double(k);
}

class PhiExpr : public dolfin::Expression {
public:
  PhiExpr(const Uint gdim, const int k) : gdim_(gdim), k_(k) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    Vector3d xx = Vector3d::Zero();
    for (Uint d = 0; d < gdim_; ++d) xx[d] = x[d];
    v[0] = phi_exact(xx, k_, gdim_);
  }
private:
  Uint gdim_;
  int k_;
};

// Two stamps, at t = 0 and 1, as a checkpoint (mesh.h5, up_k.h5) and, where
// every field is P1, as XDMF (u, p, phi .xdmf); the parameter files name
// neither wall_p2 nor include_phi, which the tests append
template<typename Cell>
void write_case(const CaseDir& c, const std::size_t n, const std::string& u_el,
                const std::string& phi_el){
  constexpr int nv = Cell::n_verts;
  constexpr Uint gdim = nv - 1;
  std::shared_ptr<dolfin::Mesh> mesh;
  if constexpr (nv == 3) mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitSquareMesh(n, n));
  else                   mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitCubeMesh(n, n, n));
  std::shared_ptr<dolfin::FunctionSpace> V, P;
  Uint ncoeffs_u = 0, ncoeffs_p = 0;
  taylor_hood_spaces<Cell>(u_el, "P1", true, mesh, nullptr, V, P, ncoeffs_u, ncoeffs_p);
  const auto F = lagrange_space<gdim, false>(phi_el, mesh, nullptr, "phase field");
  const bool xdmf = u_el == "P1" && phi_el == "P1";

  {
    dolfin::HDF5File f(MPI_COMM_WORLD, (c.path / "mesh.h5").string(), "w");
    f.write(*mesh, "mesh");
  }
  std::unique_ptr<dolfin::XDMFFile> xf[3];
  if (xdmf){
    for (int i = 0; i < 3; ++i){
      const char* name[3] = {"u.xdmf", "p.xdmf", "phi.xdmf"};
      xf[i] = std::make_unique<dolfin::XDMFFile>(mesh->mpi_comm(), (c.path / name[i]).string());
      xf[i]->parameters["functions_share_mesh"] = true;
      xf[i]->parameters["rewrite_function_mesh"] = false;
    }
  }
  for (int k = 0; k < 2; ++k){
    dolfin::Function u(V), p(P), phi(F);
    const UExpr ue(gdim, k);
    const PExpr pe(gdim, k);
    const PhiExpr fe(gdim, k);
    u.interpolate(ue);
    p.interpolate(pe);
    phi.interpolate(fe);
    dolfin::HDF5File f(MPI_COMM_WORLD, (c.path / ("up_" + std::to_string(k) + ".h5")).string(), "w");
    f.write(u, "u");
    f.write(p, "p");
    f.write(phi, "phi");
    if (xdmf){
      xf[0]->write(u, double(k));
      xf[1]->write(p, double(k));
      xf[2]->write(phi, double(k));
    }
  }
  for (auto& f : xf) if (f) f->close();

  std::ofstream(c.path / "timestamps.dat") << "0\tup_0.h5\n1\tup_1.h5\n";
  std::ofstream(c.file("h5_params.dat"))
    << "velocity_space=" << u_el << "\npressure_space=P1\ntimestamps=timestamps.dat\nmesh=mesh.h5\n";
  if (xdmf)
    std::ofstream(c.file("xdmf_params.dat")) << "u=u.xdmf\np=p.xdmf\nphi=phi.xdmf\n";
}

void append(const std::string& params, const std::string& lines){
  std::ofstream(params, std::ios::app) << lines;
}

// Points spread over the box, a quarter of them in the first layer of cells
// above the walls, where the near-wall rule acts
std::vector<Vector3d> points(const Uint gdim, const std::size_t n, const double layer){
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> any(0.01, 0.99), low(0.001, layer);
  std::vector<Vector3d> pts;
  for (std::size_t i = 0; i < n; ++i){
    Vector3d x = Vector3d::Zero();
    for (Uint d = 0; d < gdim; ++d) x[d] = any(rng);
    if (i % 4 == 0) x[i/4 % gdim] = low(rng);
    pts.push_back(x);
  }
  return pts;
}

// Everything an evaluation gives, at these points and times
struct Sample {
  std::vector<double> u, gradu, p, phi;
  std::vector<int> cell_type;
};

template<typename I>
Sample sample(I& intp, const std::vector<Vector3d>& pts, const Uint gdim){
  intp.set_int_order(2);   // the gradient too
  Sample s;
  for (const double t : {0., 0.3, 1.}){
    intp.update(t);
    for (const Vector3d& x : pts){
      CellPos pos;
      REQUIRE(intp.locate(x, t, pos));
      PointValues got(1.);
      intp.evaluate(x, t, pos, got);
      for (Uint d = 0; d < gdim; ++d) s.u.push_back(got.get_u()[d]);
      for (Uint i = 0; i < gdim; ++i)
        for (Uint j = 0; j < gdim; ++j) s.gradu.push_back(got.gradU(i, j));
      s.p.push_back(got.get_p());
      s.phi.push_back(got.get_phi());
      s.cell_type.push_back(got.get_cell_type());
    }
  }
  return s;
}

template<typename Cell>
void check_formats_agree(const std::string& tag, const std::size_t n){
  constexpr Uint gdim = Cell::n_verts - 1;
  CaseDir c(tag);
  write_case<Cell>(c, n, "P1", "P1");
  // has_phase_field decides interpol's phi column: include_phi's, whatever the file holds
  REQUIRE_FALSE(SimplexInterpol<Cell>(c.file("h5_params.dat")).has_phase_field());
  REQUIRE_FALSE(XDMFInterpol<Cell>(c.file("xdmf_params.dat")).has_phase_field());
  append(c.file("h5_params.dat"), "include_phi=true\nwall_p2=edge\n");
  append(c.file("xdmf_params.dat"), "include_phi=true\nwall_p2=edge\n");
  const std::vector<Vector3d> pts = points(gdim, 200, 0.5/double(n));

  SimplexInterpol<Cell> h5(c.file("h5_params.dat"));
  XDMFInterpol<Cell> xdmf(c.file("xdmf_params.dat"));
  REQUIRE(h5.has_phase_field());
  REQUIRE(xdmf.has_phase_field());
  REQUIRE_FALSE(h5.has_phase_gradient());
  REQUIRE_FALSE(xdmf.has_phase_gradient());
  const Sample a = sample(h5, pts, gdim), b = sample(xdmf, pts, gdim);
  REQUIRE(a.u == b.u);
  REQUIRE(a.gradu == b.gradu);
  REQUIRE(a.p == b.p);
  REQUIRE(a.phi == b.phi);
  REQUIRE(a.cell_type == b.cell_type);

  // the rule acted, and the cells next to the walls were labelled
  CaseDir d(tag + "_none");
  write_case<Cell>(d, n, "P1", "P1");
  append(d.file("h5_params.dat"), "include_phi=true\n");
  SimplexInterpol<Cell> p1(d.file("h5_params.dat"));
  const Sample l = sample(p1, pts, gdim);
  REQUIRE(l.u != a.u);
  REQUIRE(l.p == a.p);
  std::size_t labelled = 0;
  for (const int t : a.cell_type) if (t > 0) ++labelled;
  REQUIRE(labelled > 0);
}

}  // namespace

TEST_CASE("A P1 field reads the same as a checkpoint and as XDMF", "[stamped]") {
  // One evaluation behind both formats: velocity, gradient, pressure, phase
  // field and cell type agree to the last bit, near-wall P2 included
  SECTION("triangles"){ check_formats_agree<Triangle>("tri_agree", 8); }
  SECTION("tets"){ check_formats_agree<Tet>("tet_agree", 4); }
}

TEST_CASE("wall_p2 = edge on a P2 checkpoint is refused", "[stamped]") {
  // The rule makes a P1 velocity quadratic next to the walls; a P2 one is
  // already, so the key would do nothing
  CaseDir c("p2_edge");
  write_case<Triangle>(c, 4, "P2", "P1");
  append(c.file("h5_params.dat"), "wall_p2=edge\n");
  REQUIRE_THROWS_AS(SimplexInterpol<Triangle>(c.file("h5_params.dat")), partrac::Error);
}

TEST_CASE("A checkpoint's phase field is read in its own element", "[stamped]") {
  // A P2 phase field beside a P1 pressure, and one read with the pressure
  // ignored: either way its own basis, so the quadratic is exact
  SECTION("P2 phase field, P1 pressure"){
    CaseDir c("phi_p2");
    write_case<Tet>(c, 3, "P2", "P2");
    append(c.file("h5_params.dat"), "include_phi=true\n");
    SimplexInterpol<Tet> intp(c.file("h5_params.dat"));
    for (const double t : {0., 0.3, 1.}){
      intp.update(t);
      for (const Vector3d& x : points(3, 40, 0.1)){
        CellPos pos;
        REQUIRE(intp.locate(x, t, pos));
        PointValues got(1.);
        intp.evaluate(x, t, pos, got);
        const double want = (1. - t)*phi_exact(x, 0, 3) + t*phi_exact(x, 1, 3);
        REQUIRE(std::abs(got.get_phi() - want) < 1e-12);
      }
    }
  }
  SECTION("P1 phase field, pressure ignored"){
    CaseDir c("phi_nop");
    write_case<Triangle>(c, 6, "P1", "P1");
    append(c.file("h5_params.dat"), "include_phi=true\n");
    SimplexInterpol<Triangle> with_p(c.file("h5_params.dat"));
    const Sample a = sample(with_p, points(2, 40, 0.1), 2);
    append(c.file("h5_params.dat"), "ignore_pressure=true\n");
    SimplexInterpol<Triangle> without_p(c.file("h5_params.dat"));
    const Sample b = sample(without_p, points(2, 40, 0.1), 2);
    REQUIRE(a.phi == b.phi);
    REQUIRE(a.u == b.u);
  }
}

TEST_CASE("An XDMF phase field is read with the pressure ignored", "[stamped]") {
  // The phase field has its own basis; the pressure's is not computed
  CaseDir c("xdmf_phi_nop");
  write_case<Triangle>(c, 6, "P1", "P1");
  append(c.file("xdmf_params.dat"), "include_phi=true\n");
  XDMFInterpol<Triangle> with_p(c.file("xdmf_params.dat"));
  const Sample a = sample(with_p, points(2, 40, 0.1), 2);
  append(c.file("xdmf_params.dat"), "ignore_pressure=true\n");
  XDMFInterpol<Triangle> without_p(c.file("xdmf_params.dat"));
  const Sample b = sample(without_p, points(2, 40, 0.1), 2);
  REQUIRE(a.phi == b.phi);
}

TEST_CASE("The mesh cache carries the phase field", "[stamped]") {
  // A cached load holds the first stamp's phase field and its node table, and
  // reads the later stamp's through the cached mapping
  CaseDir c("phi_cache");
  write_case<Tet>(c, 3, "P2", "P2");
  append(c.file("h5_params.dat"), "include_phi=true\nmesh_cache=true\n");
  const std::vector<Vector3d> pts = points(3, 40, 0.1);
  Sample fresh, cached;
  { SimplexInterpol<Tet> intp(c.file("h5_params.dat")); fresh = sample(intp, pts, 3); }
  REQUIRE(std::filesystem::exists(c.path / "mesh_partrac_tet.h5"));
  { SimplexInterpol<Tet> intp(c.file("h5_params.dat")); cached = sample(intp, pts, 3); }
  REQUIRE(cached.phi == fresh.phi);
  REQUIRE(cached.u == fresh.u);
  REQUIRE(cached.p == fresh.p);
}

#endif
