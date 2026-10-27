// The two stamp formats behind one evaluation: a P1 field written by dolfin as
// a checkpoint and as XDMF has to read the same through either loader, with
// the near-wall rule, the phase field and the cell types, since only the
// reading of the files differs. Around that, what the checkpoint side gained
// with it: wall_upgrade on a P2 velocity acting in the wall cells only, a phase field in
// an element of its own, and a mesh cache that carries the phase field.
#ifdef USE_DOLFIN
#include <catch2/catch.hpp>

#include <array>
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
#include "held_eval.hpp"
#include "taylor_hood.hpp"

namespace {

// Stamp k of a velocity at rest on the boundary of the unit box, so the
// near-wall rule applies at every wall; the factors keep it off every symmetry.
// A drift added to every component moves the walls
class UExpr : public dolfin::Expression {
public:
  UExpr(const Uint gdim, const int k, const double drift = 0.)
    : dolfin::Expression(gdim), gdim_(gdim), k_(k), drift_(drift) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    double b = 1.;
    for (Uint d = 0; d < gdim_; ++d) b *= x[d]*(1. - x[d]);
    for (Uint d = 0; d < gdim_; ++d)
      v[d] = b*(1. + double(k_))*(d == 0 ? 4.*(1. + x[1]) : 1. + x[0] + 0.5*double(d)) + drift_;
  }
private:
  Uint gdim_;
  int k_;
  double drift_;
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

// Stamps at t = 0, 1, ..., two by default, as a checkpoint (mesh.h5, up_k.h5)
// and, where every field is P1, as XDMF (u, p, phi .xdmf); the parameter files
// name neither wall_upgrade nor include_phi, which the tests append. Stamp k
// holds UExpr(ks[k], drift[k])
template<typename Cell>
void write_case(const CaseDir& c, const std::size_t n, const std::string& u_el,
                const std::string& phi_el, const std::vector<int>& ks = {0, 1},
                const std::vector<double>& drift = {0., 0.}){
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
  std::ofstream ts(c.path / "timestamps.dat");
  for (int k = 0; k < int(ks.size()); ++k){
    ts << k << "\tup_" << k << ".h5\n";
    dolfin::Function u(V), p(P), phi(F);
    const UExpr ue(gdim, ks[std::size_t(k)], drift[std::size_t(k)]);
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
  append(c.file("h5_params.dat"), "include_phi=true\nwall_upgrade=true\n");
  append(c.file("xdmf_params.dat"), "include_phi=true\nwall_upgrade=true\n");
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

// The upgraded and the plain P2 field of one case, at the same points
template<typename Cell>
void check_upgrade_p2(const std::string& tag, const std::size_t n){
  constexpr Uint gdim = Cell::n_verts - 1;
  CaseDir c(tag);
  write_case<Cell>(c, n, "P2", "P1");
  CaseDir d(tag + "_plain");
  write_case<Cell>(d, n, "P2", "P1");
  append(c.file("h5_params.dat"), "wall_upgrade=true\n");
  const std::vector<Vector3d> pts = points(gdim, 400, 0.5/double(n));
  SimplexInterpol<Cell> up(c.file("h5_params.dat")), plain(d.file("h5_params.dat"));
  const Sample a = sample(up, pts, gdim), b = sample(plain, pts, gdim);
  // Off the boundary's layer of lattice boxes, no cell has a wall vertex
  const double reach = 1./double(n) + 1e-9;
  std::size_t far = 0, changed = 0;
  for (std::size_t q = 0; q < a.cell_type.size(); ++q){
    const Vector3d& x = pts[q % pts.size()];
    double dist = 1.;
    for (Uint k = 0; k < gdim; ++k) dist = std::min({dist, x[k], 1. - x[k]});
    bool same = true;
    for (Uint k = 0; k < gdim; ++k) same = same && a.u[q*gdim + k] == b.u[q*gdim + k];
    for (Uint k = 0; k < gdim*gdim; ++k) same = same && a.gradu[q*gdim*gdim + k] == b.gradu[q*gdim*gdim + k];
    if (dist > reach){
      REQUIRE(same);
      ++far;
    }
    else if (!same) ++changed;
  }
  REQUIRE(far > 0);
  REQUIRE(changed > 0);
  // No slip on the walls: every term vanishes on a facet at rest
  up.update(0.3);
  std::mt19937 rng(3);
  std::uniform_real_distribution<double> any(0.05, 0.95);
  for (int i = 0; i < 50; ++i){
    Vector3d x = Vector3d::Zero();
    for (Uint k = 0; k < gdim; ++k) x[k] = any(rng);
    x[Uint(i) % gdim] = i % 2 ? 0. : 1.;
    CellPos pos;
    REQUIRE(up.locate(x, 0.3, pos));
    PointValues got(1.);
    up.evaluate(x, 0.3, pos, got);
    REQUIRE(got.get_u().norm() < 1e-14);
  }
}

TEST_CASE("wall_upgrade = true on a P2 checkpoint changes the wall cells only", "[stamped]") {
  // The quartic terms live in the cells with a wall vertex at rest; every
  // other cell reads the plain P2 field bit for bit, and the wall stays at rest
  SECTION("triangles"){ check_upgrade_p2<Triangle>("p2_up_tri", 6); }
  SECTION("tets"){ check_upgrade_p2<Tet>("p2_up_tet", 4); }
}

// The fields of a located point
template<typename Cell>
PointValues located(SimplexInterpol<Cell>& intp, const Vector3d& x, const double t){
  CellPos pos;
  REQUIRE(intp.locate(x, t, pos));
  PointValues got(1.);
  intp.evaluate(x, t, pos, got);
  return got;
}

// The same held in its cell, as RK4cells evaluates
template<typename Cell>
PointValues held(SimplexInterpol<Cell>& intp, const Vector3d& x, const double t){
  CellPos pos;
  REQUIRE(intp.locate(x, t, pos));
  PointValues got(1.);
  HeldEval<SimplexInterpol<Cell>> ev(intp, got);
  ev.bind(intp.region_in(pos.id, x));
  ev.at(x, t);
  return got;
}

template<typename M>
double largest(const Eigen::MatrixBase<M>& m){ return m.cwiseAbs().maxCoeff(); }

// intp at t between its stamps at tp and tn, whose fields prev and next hold at
// t_ref: the blend of the two, their rate, and the held evaluation the same
template<typename Cell>
void check_blend(SimplexInterpol<Cell>& intp, SimplexInterpol<Cell>& prev, SimplexInterpol<Cell>& next,
                 const double tp, const double tn, const double t, const double t_ref,
                 const std::vector<Vector3d>& pts){
  const double a = (t - tp)/(tn - tp);
  double scale = 0., err = 0., err_held = 0.;
  for (const Vector3d& x : pts){
    const PointValues b = located(intp, x, t), h = held(intp, x, t);
    const PointValues p = located(prev, x, t_ref), n = located(next, x, t_ref);
    const Vector3d u = (1. - a)*p.U + a*n.U, A = (n.U - p.U)/(tn - tp);
    const Matrix3d g = (1. - a)*p.gradU + a*n.gradU, gA = (n.gradU - p.gradU)/(tn - tp);
    scale = std::max({scale, largest(u), largest(g), largest(A), largest(gA)});
    err = std::max({err, largest(b.U - u), largest(b.gradU - g), largest(b.A - A), largest(b.gradA - gA)});
    err_held = std::max({err_held, largest(h.U - b.U), largest(h.gradU - b.gradU), largest(h.A - b.A),
                         largest(h.gradA - b.gradA)});
  }
  INFO("t = " << t << ", scale " << scale);
  REQUIRE(scale > 0.);
  REQUIRE(err <= 1e-13*scale);
  REQUIRE(err_held <= 1e-13*scale);
}

// Stamps at rest on the walls blended with ones whose walls move: a moving
// stamp reads its plain field (P1, or P2 without terms), and between two
// stamps the field is the blend of the one at rest's raised field (the edge
// rule on P1, the terms on P2) and the moving one's plain field, and its rate
// their difference over the interval, in either order, located or held; past
// a stamp, the buffer's swap keeps each stamp's raised reading with its values
template<typename Cell>
void check_upgrade_moving(const std::string& tag, const std::size_t n, const std::string& u_el){
  constexpr Uint gdim = Cell::n_verts - 1;
  constexpr double drift = 0.05, t = 0.3, t_ref = 0.5;
  const auto params = [&](const CaseDir& c, const std::vector<int>& ks, const std::vector<double>& d,
                          const bool upgrade){
    write_case<Cell>(c, n, u_el, "P1", ks, d);
    append(c.file("h5_params.dat"), upgrade ? "wall_upgrade=true\n" : "wall_upgrade=false\n");
    return c.file("h5_params.dat");
  };
  CaseDir cb(tag), cv(tag + "_reverse"), c3(tag + "_three"), cr(tag + "_rest"), cr2(tag + "_rest2"),
          cm(tag + "_moving"), cp(tag + "_moving_plain"), cq(tag + "_rest_plain");
  SimplexInterpol<Cell> both(params(cb, {0, 1}, {0., drift}, true));
  SimplexInterpol<Cell> reverse(params(cv, {1, 0}, {drift, 0.}, true));
  SimplexInterpol<Cell> three(params(c3, {0, 1, 2}, {0., drift, 0.}, true));
  SimplexInterpol<Cell> rest(params(cr, {0, 0}, {0., 0.}, true));
  SimplexInterpol<Cell> rest2(params(cr2, {2, 2}, {0., 0.}, true));
  SimplexInterpol<Cell> moving(params(cm, {1, 1}, {drift, drift}, true));
  SimplexInterpol<Cell> moving_plain(params(cp, {1, 1}, {drift, drift}, false));
  SimplexInterpol<Cell> rest_plain(params(cq, {0, 0}, {0., 0.}, false));
  const std::vector<Vector3d> pts = points(gdim, 400, 0.5/double(n));
  for (auto* i : {&both, &reverse, &three, &rest, &rest2, &moving, &moving_plain, &rest_plain}){
    i->set_int_order(2);
    i->set_needs_gradient(true);
    i->update(t_ref);
  }
  std::size_t raised = 0;
  for (const Vector3d& x : pts){
    const PointValues r = located(rest, x, t_ref), m = located(moving, x, t_ref);
    const PointValues mp = located(moving_plain, x, t_ref), rp = located(rest_plain, x, t_ref);
    // Where the wall moves: the plain field, bit for bit
    REQUIRE(m.U == mp.U);
    REQUIRE(m.gradU == mp.gradU);
    if (r.U != rp.U || r.gradU != rp.gradU) ++raised;
  }
  REQUIRE(raised > 0);
  both.update(t);
  check_blend(both, rest, moving, 0., 1., t, t_ref, pts);
  reverse.update(t);
  check_blend(reverse, moving, rest, 0., 1., t, t_ref, pts);
  three.update(t);
  check_blend(three, rest, moving, 0., 1., t, t_ref, pts);
  // Past stamp 1: the moving stamp becomes prev
  three.update(1. + t);
  check_blend(three, moving, rest2, 1., 2., 1. + t, t_ref, pts);
}

TEST_CASE("wall_upgrade = true blends stamps at rest with ones whose walls move", "[stamped]") {
  // Whether the walls are at rest is a stamp's own: the raised field of the
  // one at rest is blended with the plain field of the moving one
  SECTION("P2 triangles"){ check_upgrade_moving<Triangle>("p2_move_tri", 6, "P2"); }
  SECTION("P2 tets"){ check_upgrade_moving<Tet>("p2_move_tet", 4, "P2"); }
  SECTION("P1 triangles"){ check_upgrade_moving<Triangle>("p1_move_tri", 6, "P1"); }
  SECTION("P1 tets"){ check_upgrade_moving<Tet>("p1_move_tet", 4, "P1"); }
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
