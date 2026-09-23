// SplitInterpol on dolfin-written checkpoints: the loader around
// split_eval.hpp -- the element check, the net-flux check every stamp must
// pass, the interior values held with the stamp, and the blend in time.
//
// The fixture needs no cleaner: a quadratic divergence-free velocity
// interpolated to P2 has exactly zero net flux in every cell, and the split
// field is then the polynomial itself, so the loader must return it to
// round-off at every point, value and gradient. Beside that, what a loader owes
// its parameter file: a P1 velocity, an unprepared file, a non-finite value,
// wall_p2 and mesh_cache refused; a phase field in an element of its own; a
// periodic mesh.
#ifdef USE_DOLFIN
#include <catch2/catch.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <dolfin.h>

#include "Error.hpp"
#include "SimplexInterpol.hpp"
#include "SplitInterpol.hpp"
#include "Tet.hpp"
#include "Triangle.hpp"
#include "case_dir.hpp"
#include "divfree_poly.hpp"
#include "dolfin_ref.hpp"
#include "taylor_hood.hpp"

namespace {

// The same field, independent of x, so a mesh periodic along x carries it: the
// transverse components still balance, since only y and z enter the divergence
Vector3d per_u(const Vector3d& p, const Uint gdim){
  const double y = p[1], z = p[2];
  if (gdim == 2)
    return {0.8*y*y - 0.4*y + 0.35, 0.6, 0.};
  return {0.3*y*y - 0.5*z*z + 0.7*y*z + 0.2*y + 0.15, 0.9*y*z + 0.25*y,
          -0.45*z*z - 0.25*z + 0.4};
}

Matrix3d per_grad(const Vector3d& p, const Uint gdim){
  const double y = p[1], z = p[2];
  Matrix3d g = Matrix3d::Zero();
  if (gdim == 2){
    g(0, 1) = 1.6*y - 0.4;
    return g;
  }
  g(0, 1) = 0.6*y + 0.7*z + 0.2;   g(0, 2) = -1.0*z + 0.7*y;
  g(1, 1) = 0.9*z + 0.25;          g(1, 2) = 0.9*y;
  g(2, 2) = -0.9*z - 0.25;
  return g;
}

// Stamp k of the velocity: the polynomial times 1 + k, so the blend and the
// rate in time are the polynomial again
class UExpr : public dolfin::Expression {
public:
  UExpr(const Uint gdim, const int k, const bool periodic)
    : dolfin::Expression(gdim), gdim_(gdim), k_(k), periodic_(periodic) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    Vector3d xx = Vector3d::Zero();
    for (Uint d = 0; d < gdim_; ++d) xx[d] = x[d];
    const Vector3d u = periodic_ ? per_u(xx, gdim_) : poly_u(xx, gdim_);
    for (Uint d = 0; d < gdim_; ++d) v[d] = (1. + double(k_))*u[d];
  }
private:
  Uint gdim_;
  int k_;
  bool periodic_;
};

// A velocity with a divergence, so its cells do not balance
class DirtyExpr : public dolfin::Expression {
public:
  explicit DirtyExpr(const Uint gdim) : dolfin::Expression(gdim), gdim_(gdim) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    for (Uint d = 0; d < gdim_; ++d)
      v[d] = x[d]*x[d] + 0.3*x[(d + 1) % gdim_] + 0.1*double(d);
  }
private:
  Uint gdim_;
};

class PExpr : public dolfin::Expression {
public:
  explicit PExpr(const Uint gdim) : gdim_(gdim) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    v[0] = 0.7*x[0] - 0.4*x[1] + (gdim_ == 3 ? 0.3*x[2] : 0.) + 0.2;
  }
private:
  Uint gdim_;
};

// A quadratic, so a P2 space holds it exactly
double phi_exact(const Vector3d& x, const Uint gdim){
  return x[0]*x[0] - x[0]*x[1] + 0.5*x[1]*x[1] + (gdim == 3 ? x[2]*x[2] - 0.3*x[2] : 0.);
}

class PhiExpr : public dolfin::Expression {
public:
  explicit PhiExpr(const Uint gdim) : gdim_(gdim) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    Vector3d xx = Vector3d::Zero();
    for (Uint d = 0; d < gdim_; ++d) xx[d] = x[d];
    v[0] = phi_exact(xx, gdim_);
  }
private:
  Uint gdim_;
};

// Every edge midpoint moved along its own edge. A facet's flux is a Simpson
// rule over the facet's nodes, and a move along the facet changes no flux, so
// every cell still balances; but the field is no longer the trace of a smooth
// divergence-free one, which is what a squashed cell amplifies. The move is a
// function of the edge alone, so the two cells sharing it write the same value.
std::size_t nudge_midpoints(dolfin::Function& u, const dolfin::FunctionSpace& V,
                            const Uint gdim, const double amp)
{
  const std::vector<double> dc = V.tabulate_dof_coordinates();
  std::vector<double> vals;
  u.vector()->get_local(vals);
  const dolfin::Mesh& mesh = *V.mesh();
  std::size_t moved = 0;
  for (dolfin::CellIterator c(mesh); !c.end(); ++c){
    std::vector<Vector3d> vx;
    for (dolfin::VertexIterator v(*c); !v.end(); ++v){
      Vector3d p = Vector3d::Zero();
      for (Uint d = 0; d < gdim; ++d) p[d] = v->x(d);
      vx.push_back(p);
    }
    const auto dofs = V.dofmap()->cell_dofs(c->index());
    const auto at = [&](const int k){
      Vector3d p = Vector3d::Zero();
      for (Uint d = 0; d < gdim; ++d) p[d] = dc[std::size_t(dofs[k])*gdim + d];
      return p;
    };
    // A cell's dofs are one component after another, as block_value reads them
    const int nn = int(dofs.size())/int(gdim);
    for (int k = 0; k < nn; ++k){
      for (Uint d = 1; d < gdim; ++d)
        REQUIRE((at(k + int(d)*nn) - at(k)).norm() == 0.);
      const Vector3d x = at(k);
      // The vertex pair this node is the midpoint of, ordered by position so
      // that both of the edge's cells build the same tangent
      for (std::size_t a = 0; a < vx.size(); ++a)
        for (std::size_t b = a + 1; b < vx.size(); ++b){
          const bool a_first = std::lexicographical_compare(vx[a].data(), vx[a].data() + 3,
                                                            vx[b].data(), vx[b].data() + 3);
          const Vector3d& lo = a_first ? vx[a] : vx[b];
          const Vector3d& hi = a_first ? vx[b] : vx[a];
          const Vector3d mid = 0.5*(lo + hi);
          if ((mid - x).norm() > 1e-12*(1. + x.norm())) continue;
          const Vector3d tang = (hi - lo).normalized();
          const double w = amp*std::sin(13.*mid[0] + 7.*mid[1]);
          for (Uint d = 0; d < gdim; ++d) vals[std::size_t(dofs[k + int(d)*nn])] += w*tang[d];
          ++moved;
        }
    }
  }
  u.vector()->set_local(vals);
  u.vector()->apply("insert");
  return moved;
}

// A velocity (g(y), 0, ...) whose nodal value depends on the node's y alone. A
// cell's three y levels carry one quadratic in y, so on any mesh the P2
// interpolant of such data is that quadratic: div u = 0 pointwise and every net
// flux is exactly zero, whatever g. Here g is 1 above y = 0.5 and `quiet`
// below, so the lower cells are numerically dead -- the bead pack's dead-end
// pores, where the Stokes velocity is at the solver's noise floor and a cell's
// own facet fluxes are round-off. `kick` is then added to the y component of
// one midpoint down there, and is the stamp's only net flux.
void layered_velocity(dolfin::Function& u, const dolfin::FunctionSpace& V,
                      const Uint gdim, const std::size_t n,
                      const double quiet, const double kick)
{
  const std::vector<double> dc = V.tabulate_dof_coordinates();
  std::vector<double> vals;
  u.vector()->get_local(vals);
  const dolfin::Mesh& mesh = *V.mesh();
  const double h = 1./double(n);
  std::ptrdiff_t kicked = -1;
  std::vector<std::ptrdiff_t> kicked_dofs(gdim, 0);
  for (dolfin::CellIterator c(mesh); !c.end(); ++c){
    const auto dofs = V.dofmap()->cell_dofs(c->index());
    const int nn = int(dofs.size())/int(gdim);
    for (int k = 0; k < nn; ++k){
      Vector3d x = Vector3d::Zero();
      for (Uint d = 0; d < gdim; ++d) x[d] = dc[std::size_t(dofs[k])*gdim + d];
      for (Uint d = 0; d < gdim; ++d)
        vals[std::size_t(dofs[k + int(d)*nn])] = d == 0 ? (x[1] < 0.5 ? quiet : 1.) : 0.;
      // a midpoint of the bottom row of cells, so every cell holding it is dead
      if (kicked < 0 && x[1] > 0.2*h && x[1] < 0.8*h && x[0] > 0.2 && x[0] < 0.8){
        kicked = k;
        for (Uint d = 0; d < gdim; ++d) kicked_dofs[d] = dofs[k + int(d)*nn];
      }
    }
  }
  REQUIRE(kicked >= 0);
  // a facet's flux picks up the moved node through its own normal, so the move
  // is off every facet normal of a structured mesh
  for (Uint d = 0; d < gdim; ++d)
    vals[std::size_t(kicked_dofs[d])] += kick*double(1 << (2*d));
  u.vector()->set_local(vals);
  u.vector()->apply("insert");
}

// What a case writes: the element of the velocity, how many stamps, whether the
// field is the periodic one and whether a phase field goes beside it
struct CaseOpts {
  std::string u_el = "P2";
  std::string phi_el = "";     // empty: no phase field
  int nstamps = 2;
  bool periodic = false;       // the mesh and the field periodic along x
  bool dirty = false;          // a velocity whose cells do not balance
  double nudge = 0.;           // every edge midpoint moved along its own edge
  double squash = 0.;          // the mesh flattened along y by this factor
  double quiet = 0.;           // the velocity below y = 0.5 scaled to this
  double kick = 0.;            // added to one midpoint's y component down there
  std::optional<double> poison;  // one velocity value of the last stamp replaced by this
  bool poison_first = false;     // ... of the first stamp instead
};

template<typename Cell>
void write_case(const CaseDir& c, const std::size_t n, const CaseOpts& o){
  constexpr int nv = Cell::n_verts;
  constexpr Uint gdim = nv - 1;
  std::shared_ptr<dolfin::Mesh> mesh;
  if constexpr (nv == 3) mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitSquareMesh(n, n));
  else                   mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitCubeMesh(n, n, n));
  if (o.squash > 0.){
    std::vector<double>& xs = mesh->coordinates();
    for (std::size_t k = 1; k < xs.size(); k += gdim) xs[k] *= o.squash;
  }
  std::shared_ptr<const dolfin::SubDomain> pbc;
  if (o.periodic){
    std::vector<bool> per(3, false);
    per[0] = true;
    pbc = std::make_shared<PeriodicBC>(per, Vector3d::Zero(), Vector3d(1., 1., 1.), gdim);
  }
  std::shared_ptr<dolfin::FunctionSpace> V, P;
  Uint ncoeffs_u = 0, ncoeffs_p = 0;
  taylor_hood_spaces<Cell>(o.u_el, "P1", true, mesh, pbc, V, P, ncoeffs_u, ncoeffs_p);
  std::shared_ptr<dolfin::FunctionSpace> F;
  if (!o.phi_el.empty())
    F = lagrange_space<gdim, false>(o.phi_el, mesh, pbc, "phase field");

  {
    dolfin::HDF5File f(MPI_COMM_WORLD, (c.path / "mesh.h5").string(), "w");
    f.write(*mesh, "mesh");
  }
  std::ofstream stamps(c.path / "timestamps.dat");
  for (int k = 0; k < o.nstamps; ++k){
    dolfin::Function u(V), p(P);
    const UExpr ue(gdim, k, o.periodic);
    const DirtyExpr de(gdim);
    const PExpr pe(gdim);
    if (o.dirty) u.interpolate(de); else u.interpolate(ue);
    if (o.quiet > 0.) layered_velocity(u, *V, gdim, n, o.quiet, o.kick);
    if (o.nudge > 0.) REQUIRE(nudge_midpoints(u, *V, gdim, o.nudge) > 0);
    if (o.poison && k == (o.poison_first ? 0 : o.nstamps - 1)){
      std::vector<double> vals;
      u.vector()->get_local(vals);
      vals[vals.size()/2] = *o.poison;
      u.vector()->set_local(vals);
      u.vector()->apply("insert");
    }
    p.interpolate(pe);
    dolfin::HDF5File f(MPI_COMM_WORLD, (c.path / ("up_" + std::to_string(k) + ".h5")).string(), "w");
    f.write(u, "u");
    f.write(p, "p");
    if (F){
      dolfin::Function phi(F);
      const PhiExpr fe(gdim);
      phi.interpolate(fe);
      f.write(phi, "phi");
    }
    stamps << k << "\tup_" << k << ".h5\n";
  }
  stamps.close();
  std::ofstream(c.params())
    << "velocity_space=" << o.u_el << "\npressure_space=P1\ntimestamps=timestamps.dat\n"
    << "mesh=mesh.h5\ndivfree=true\n"
    << "periodic_x=" << (o.periodic ? "true" : "false") << "\n";
}

void append(const std::string& params, const std::string& lines){
  std::ofstream(params, std::ios::app) << lines;
}

// Points spread over the box, a quarter of them in the first layer of cells
// above the walls, where no slip is what the construction buys
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

template<typename Cell>
void check_returns_polynomial(const std::string& tag, const std::size_t n){
  constexpr Uint gdim = Cell::n_verts - 1;
  CaseDir c(tag);
  write_case<Cell>(c, n, CaseOpts{});
  SplitInterpol<Cell> intp(c.params());
  intp.set_int_order(2);
  const std::vector<Vector3d> pts = points(gdim, 300, 0.5/double(n));
  for (const double t : {0., 0.4, 1.}){
    intp.update(t);
    const double s = 1. + t;   // stamp k is the polynomial times 1 + k
    for (const Vector3d& x : pts){
      CellPos pos;
      REQUIRE(intp.locate(x, t, pos));
      PointValues got(1.);
      intp.evaluate(x, t, pos, got);
      const Vector3d u = poly_u(x, gdim);
      const Matrix3d g = poly_grad(x, gdim);
      for (Uint d = 0; d < gdim; ++d){
        REQUIRE(got.get_u()[d] == Approx(s*u[d]).margin(1e-11));
        // Stamp k is (1 + k) times one field, so the rate in time is the field
        REQUIRE(got.get_a()[d] == Approx(u[d]).margin(1e-11));
        for (Uint j = 0; j < gdim; ++j){
          REQUIRE(got.gradU(d, j) == Approx(s*g(d, j)).margin(1e-9));
          REQUIRE(got.gradA(d, j) == Approx(g(d, j)).margin(1e-9));
        }
      }
      // div u = 0 is what the construction is for
      double div = 0.;
      for (Uint d = 0; d < gdim; ++d) div += got.gradU(d, d);
      REQUIRE(std::abs(div) < 1e-9);
    }
  }
}

template<typename Cell>
void check_motion_and_scalars(const std::string& tag, const std::size_t n){
  constexpr Uint gdim = Cell::n_verts - 1;
  CaseDir c(tag);
  CaseOpts o;
  o.phi_el = "P2";   // an element of its own, beside a P1 pressure
  write_case<Cell>(c, n, o);
  append(c.params(), "include_phi=true\n");
  SplitInterpol<Cell> intp(c.params());
  intp.set_int_order(2);
  const std::vector<Vector3d> pts = points(gdim, 120, 0.5/double(n));
  const double t = 0.3;
  intp.update(t);
  for (const Vector3d& x : pts){
    CellPos pos;
    REQUIRE(intp.locate(x, t, pos));
    PointValues full(1.), motion(1.);
    intp.evaluate(x, t, pos, full);
    intp.evaluate_motion(x, t, pos, motion);
    for (Uint d = 0; d < gdim; ++d){
      REQUIRE(motion.get_u()[d] == full.get_u()[d]);
      REQUIRE(motion.get_a()[d] == full.get_a()[d]);
      for (Uint j = 0; j < gdim; ++j){
        REQUIRE(motion.gradU(d, j) == full.gradU(d, j));
        REQUIRE(motion.gradA(d, j) == full.gradA(d, j));
      }
    }
    // evaluate_motion leaves the scalars alone
    REQUIRE(motion.get_phi() == 0.);
    REQUIRE(full.get_phi() == Approx(phi_exact(x, gdim)).margin(1e-11));
    REQUIRE(full.get_p() == Approx(0.7*x[0] - 0.4*x[1] + (gdim == 3 ? 0.3*x[2] : 0.) + 0.2).margin(1e-12));
  }
}

template<typename Cell>
void check_one_stamp(const std::string& tag, const std::size_t n){
  constexpr Uint gdim = Cell::n_verts - 1;
  CaseDir c(tag);
  CaseOpts o;
  o.nstamps = 1;
  write_case<Cell>(c, n, o);
  SplitInterpol<Cell> intp(c.params());
  intp.update(0.);
  REQUIRE(intp.stamps_aliased());
  const std::vector<Vector3d> pts = points(gdim, 60, 0.5/double(n));
  for (const Vector3d& x : pts){
    CellPos pos;
    REQUIRE(intp.locate(x, 0., pos));
    PointValues got(1.);
    intp.evaluate(x, 0., pos, got);
    const Vector3d u = poly_u(x, gdim);
    for (Uint d = 0; d < gdim; ++d){
      REQUIRE(got.get_u()[d] == Approx(u[d]).margin(1e-11));
      REQUIRE(got.get_a()[d] == 0.);
    }
  }
}

template<typename Cell>
void check_periodic(const std::string& tag, const std::size_t n){
  constexpr Uint gdim = Cell::n_verts - 1;
  CaseDir c(tag);
  CaseOpts o;
  o.periodic = true;
  write_case<Cell>(c, n, o);
  SplitInterpol<Cell> intp(c.params());
  intp.set_int_order(2);
  intp.update(0.);
  std::mt19937 rng(3);
  std::uniform_real_distribution<double> any(0.01, 0.99), seam(0., 1e-3);
  for (int i = 0; i < 300; ++i){
    Vector3d x = Vector3d::Zero();
    for (Uint d = 0; d < gdim; ++d) x[d] = any(rng);
    // A quarter of them within a hair of the seam, on either side
    if (i % 4 == 0) x[0] = (i % 8 == 0) ? seam(rng) : 1. - seam(rng);
    CellPos pos;
    REQUIRE(intp.locate(x, 0., pos));
    PointValues got(1.);
    intp.evaluate(x, 0., pos, got);
    const Vector3d u = per_u(x, gdim);
    const Matrix3d g = per_grad(x, gdim);
    for (Uint d = 0; d < gdim; ++d){
      REQUIRE(got.get_u()[d] == Approx(u[d]).margin(1e-11));
      for (Uint j = 0; j < gdim; ++j)
        REQUIRE(got.gradU(d, j) == Approx(g(d, j)).margin(1e-9));
    }
  }
}

template<typename Cell>
void check_refusals(const std::string& tag, const std::size_t n){
  {
    // A file no cleaner has been over
    CaseDir c(tag + "_dirty");
    CaseOpts o;
    o.dirty = true;
    write_case<Cell>(c, n, o);
    REQUIRE_THROWS_AS(SplitInterpol<Cell>(c.params()), partrac::Error);
  }
  {
    // A P1 velocity: the construction reads a cell's P2 boundary data
    CaseDir c(tag + "_p1");
    CaseOpts o;
    o.u_el = "P1";
    write_case<Cell>(c, n, o);
    REQUIRE_THROWS_AS(SplitInterpol<Cell>(c.params()), partrac::Error);
  }
  {
    // The near-wall rule is what div u = 0 replaces
    CaseDir c(tag + "_wall");
    write_case<Cell>(c, n, CaseOpts{});
    append(c.params(), "wall_p2=edge\n");
    REQUIRE_THROWS_AS(SplitInterpol<Cell>(c.params()), partrac::Error);
  }
  {
    CaseDir c(tag + "_cache");
    write_case<Cell>(c, n, CaseOpts{});
    append(c.params(), "mesh_cache=true\n");
    REQUIRE_THROWS_AS(SplitInterpol<Cell>(c.params()), partrac::Error);
  }
}

// What was written to cerr while f ran
template<typename F>
std::string captured_cerr(F&& f){
  std::ostringstream buf;
  std::streambuf* old = std::cerr.rdbuf(buf.rdbuf());
  f();
  std::cerr.rdbuf(old);
  return buf.str();
}

}  // namespace

TEST_CASE("SplitInterpol returns a quadratic divergence-free velocity", "[split]"){
  check_returns_polynomial<Triangle>("poly2d", 6);
  check_returns_polynomial<Tet>("poly3d", 3);
}

TEST_CASE("SplitInterpol's evaluate_motion is its evaluate without the scalars", "[split]"){
  check_motion_and_scalars<Triangle>("motion2d", 6);
  check_motion_and_scalars<Tet>("motion3d", 3);
}

TEST_CASE("SplitInterpol on a single stamp aliases it", "[split]"){
  check_one_stamp<Triangle>("one2d", 6);
  check_one_stamp<Tet>("one3d", 3);
}

TEST_CASE("SplitInterpol reads a periodic mesh", "[split]"){
  check_periodic<Triangle>("per2d", 6);
  check_periodic<Tet>("per3d", 3);
}

TEST_CASE("SplitInterpol refuses what it cannot read", "[split]"){
  check_refusals<Triangle>("bad2d", 4);
  check_refusals<Tet>("bad3d", 2);
}

// The factory picks SplitInterpol from this key; the plain loader would read
// the same file and evaluate it as the trapping plain P2 field
TEST_CASE("SimplexInterpol refuses a file marked divfree", "[split]"){
  CaseDir c("plain2d");
  write_case<Triangle>(c, 4, CaseOpts{});
  REQUIRE_THROWS_AS(SimplexInterpol<Triangle>(c.params()), partrac::Error);
}

// A cell whose own facet fluxes are round-off is judged against the stamp's
// scale, not against its own: the ratio of its net flux to its own largest
// facet flux is a ratio of round-off to round-off. The tool refines against the
// same criterion, so what it writes is what loads.
TEST_CASE("SplitInterpol measures a dead cell against the stamp", "[split]"){
  CaseOpts o;
  o.nstamps = 1;
  o.quiet = 1e-13;           // the lower half of the box at the solver's noise floor
  {
    // a net flux far above 1e-9 of that cell's own fluxes and far below 1e-9 of
    // the floor, which is 1e-6 of the stamp's largest facet flux
    CaseDir c("dead2d");
    o.kick = 1e-17;
    write_case<Triangle>(c, 4, o);
    REQUIRE_NOTHROW(SplitInterpol<Triangle>(c.params()));
  }
  {
    CaseDir c("dead3d");
    o.kick = 1e-17;
    write_case<Tet>(c, 2, o);
    REQUIRE_NOTHROW(SplitInterpol<Tet>(c.params()));
  }
  {
    // the floor is a floor and not a pardon: four decades up it refuses again,
    // for the net flux, which shows the kick registers
    CaseDir c("deadbad2d");
    o.kick = 1e-13;
    write_case<Triangle>(c, 4, o);
    REQUIRE_THROWS_WITH(SplitInterpol<Triangle>(c.params()), Catch::Contains("has a net flux"));
  }
  {
    CaseDir c("deadbad3d");
    o.kick = 1e-13;
    write_case<Tet>(c, 2, o);
    REQUIRE_THROWS_WITH(SplitInterpol<Tet>(c.params()), Catch::Contains("has a net flux"));
  }
}

// A NaN compares false with every tolerance, so the net-flux refusal alone
// would pass it: a diverged solver's output must be refused by file and cell.
// A later stamp is refused by the split loader's own check, the first where
// its field is read, before the dof table would call it a disagreement.
template<typename Cell>
void check_non_finite(const std::string& tag, const std::size_t n, const double bad){
  CaseDir c(tag);
  CaseOpts o;
  o.poison = bad;
  write_case<Cell>(c, n, o);
  SplitInterpol<Cell> intp(c.params());
  REQUIRE_THROWS_WITH(intp.update(0.5),
                      Catch::Contains("up_1.h5: cell") && Catch::Contains("non-finite"));
}

TEST_CASE("SplitInterpol refuses a non-finite velocity", "[split]"){
  check_non_finite<Triangle>("nan2d", 4, std::numeric_limits<double>::quiet_NaN());
  check_non_finite<Tet>("nan3d", 2, std::numeric_limits<double>::quiet_NaN());
  check_non_finite<Triangle>("inf2d", 4, std::numeric_limits<double>::infinity());
  CaseDir c("nanfirst");
  CaseOpts o;
  o.poison = std::numeric_limits<double>::quiet_NaN();
  o.poison_first = true;
  write_case<Triangle>(c, 4, o);
  REQUIRE_THROWS_WITH(SplitInterpol<Triangle>(c.params()),
                      Catch::Contains("up_0.h5: u holds a non-finite value"));
}

// A squashed cell returns interior values in proportion to its aspect, which
// is the construction's known weakness; the loader says so once a stamp
TEST_CASE("SplitInterpol warns about a squashed cell", "[split]"){
  CaseOpts o;
  o.nudge = 0.5;
  o.squash = 1e-4;
  o.nstamps = 1;
  {
    CaseDir c("sliver2d");
    write_case<Triangle>(c, 4, o);
    const std::string said = captured_cerr([&]{ SplitInterpol<Triangle> intp(c.params()); });
    REQUIRE(said.find("Warning") != std::string::npos);
    REQUIRE(said.find("cond(J)") != std::string::npos);
    // One line, not one a cell
    REQUIRE(std::count(said.begin(), said.end(), '\n') == 1);
  }
  {
    CaseDir c("sliver3d");
    write_case<Tet>(c, 2, o);
    const std::string said = captured_cerr([&]{ SplitInterpol<Tet> intp(c.params()); });
    REQUIRE(said.find("Warning") != std::string::npos);
    REQUIRE(said.find("cond(J)") != std::string::npos);
    REQUIRE(std::count(said.begin(), said.end(), '\n') == 1);
  }
  {
    // The same field on a mesh of well-shaped cells says nothing
    CaseOpts fine = o;
    fine.squash = 0.;
    CaseDir c2("nosliver2d"), c3("nosliver3d");
    write_case<Triangle>(c2, 4, fine);
    write_case<Tet>(c3, 2, fine);
    REQUIRE(captured_cerr([&]{ SplitInterpol<Triangle> intp(c2.params()); }).empty());
    REQUIRE(captured_cerr([&]{ SplitInterpol<Tet> intp(c3.params()); }).empty());
  }
}

#endif
