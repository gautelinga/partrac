#include <catch2/catch.hpp>

#include <cmath>
#include <map>
#include <string>

#include "typedefs.hpp"
#include "utils.hpp"
#include "expressions/Expr_ABCFlow.hpp"
#include "expressions/Expr_BatchelorVortex.hpp"
#include "expressions/Expr_BrinkmanCylinder.hpp"
#include "expressions/Expr_HagenPoiseuille.hpp"
#include "expressions/Expr_InfinitePlate.hpp"
#include "expressions/Expr_PlanePoiseuille.hpp"
#include "expressions/Expr_SineFlow.hpp"
#include "expressions/Expr_StokesSphere.hpp"
#include "expressions/Expr_TaylorCouette.hpp"

namespace {

using Prm = std::map<std::string, std::string>;

// gradU(i, j) must be dU_i/dx_j
void check_gradient(Expr& e, const Vector3d& x, const double t){
  const double h = 1e-6;
  PointValues pv(1.0);
  e.eval(x, t, pv);
  for (int j = 0; j < 3; ++j){
    Vector3d xp = x, xm = x;
    xp[j] += h;
    xm[j] -= h;
    PointValues pp(1.0), pm(1.0);
    e.eval(xp, t, pp);
    e.eval(xm, t, pm);
    for (int i = 0; i < 3; ++i){
      const double fd = (pp.U[i] - pm.U[i])/(2*h);
      INFO("gradU(" << i << ", " << j << ")");
      REQUIRE(std::abs(fd - pv.gradU(i, j))
              <= 1e-5*std::max(1.0, std::abs(pv.gradU(i, j))));
    }
  }
}

// The integrator reads the PointValues overload, the initializers the accessors
void check_light_matches_heavy(Expr& e, const Vector3d& x, const double t){
  PointValues pv(1.0);
  e.eval(x, t, pv);
  e.eval(x, t);
  const double tol = 1e-12;
  REQUIRE(e.ux() == Approx(pv.U[0]).margin(tol));
  REQUIRE(e.uy() == Approx(pv.U[1]).margin(tol));
  REQUIRE(e.uz() == Approx(pv.U[2]).margin(tol));
  REQUIRE(e.p() == Approx(pv.get_p()).margin(tol));
  REQUIRE(e.rho() == Approx(pv.get_rho()).margin(tol));
  const double g[9] = {e.uxx(), e.uxy(), e.uxz(),
                       e.uyx(), e.uyy(), e.uyz(),
                       e.uzx(), e.uzy(), e.uzz()};
  for (int i = 0; i < 3; ++i){
    for (int j = 0; j < 3; ++j){
      INFO("gradU(" << i << ", " << j << ")");
      REQUIRE(g[3*i + j] == Approx(pv.gradU(i, j)).margin(tol));
    }
  }
  REQUIRE(e.inside(x, t) == e.inside());
}

void check_incompressible(Expr& e, const Vector3d& x, const double t){
  PointValues pv(1.0);
  e.eval(x, t, pv);
  REQUIRE(pv.gradU.trace() == Approx(0.).margin(1e-9));
}

// The PointValues overload runs inside an omp for, so it must be stateless
void check_stateless(Expr& e, const Vector3d& x1, const Vector3d& x2, const double t){
  PointValues first(1.0), other(1.0), again(1.0);
  e.eval(x1, t, first);
  e.eval(x2, t, other);
  e.eval(x1, t, again);
  e.eval(x2, t);  // the member-writing overload must not disturb it either
  PointValues after_light(1.0);
  e.eval(x1, t, after_light);
  for (int i = 0; i < 3; ++i){
    REQUIRE(again.U[i] == first.U[i]);
    REQUIRE(after_light.U[i] == first.U[i]);
    for (int j = 0; j < 3; ++j){
      REQUIRE(again.gradU(i, j) == first.gradU(i, j));
      REQUIRE(after_light.gradU(i, j) == first.gradU(i, j));
    }
  }
}

void check_all(Expr& e, const Vector3d& x, const Vector3d& x2, const double t){
  check_gradient(e, x, t);
  check_light_matches_heavy(e, x, t);
  check_incompressible(e, x, t);
  check_stateless(e, x, x2, t);
}

Prm abc_params(){
  return {{"A", "1.0"}, {"B", "0.7"}, {"C", "1.3"}, {"L", "10"}, {"rho", "1.2"},
          {"x0", "5.0"}, {"y0", "5.0"}, {"z0", "5.0"}, {"p_inf", "0.5"}};
}

Prm batchelor_params(){
  return {{"R1", "1.0"}, {"R2", "1.3"}, {"q", "0.8"}, {"rho", "1.2"},
          {"x0", "0.0"}, {"y0", "0.0"}, {"z0", "0.0"}, {"u0", "1.0"},
          {"p_inf", "0.5"}};
}

}  // namespace

TEST_CASE("ABCFlow", "[expr]") {
  Prm p = abc_params();
  Expr_ABCFlow e(p);
  check_all(e, {1.3, 2.7, 4.1}, {2.2, 0.4, 7.9}, 0.);

  // Beltrami: curl u = k u. This is what fixes the sign of the pressure.
  PointValues pv(1.0);
  e.eval({1.3, 2.7, 4.1}, 0., pv);
  const Vector3d curl{pv.gradU(2, 1) - pv.gradU(1, 2),
                      pv.gradU(0, 2) - pv.gradU(2, 0),
                      pv.gradU(1, 0) - pv.gradU(0, 1)};
  const double k = 2*M_PI/10.;
  for (int i = 0; i < 3; ++i)
    REQUIRE(curl[i] == Approx(k*pv.U[i]).margin(1e-12));
}

TEST_CASE("SineFlow", "[expr]") {
  Prm p = {{"flowdir", "1,0"}, {"depdir", "0,1"}, {"chi", "1.2154,3.1199"},
           {"u_inf", "0.7071"}, {"p_inf", "1.0"}, {"tau", "0.5"}, {"rho", "1.2"},
           {"Lx", "1.0"}, {"Ly", "1.0"}, {"Lz", "1.0"}};
  Expr_SineFlow e(p);
  // either side of a phase switch, since the flow depends on t
  check_all(e, {0.31, 0.62, 0.11}, {0.77, 0.05, 0.90}, 0.1);
  check_all(e, {0.31, 0.62, 0.11}, {0.77, 0.05, 0.90}, 0.7);

  PointValues lo(1.0), hi(1.0);
  e.eval({0.31, 0.62, 0.11}, 0.1, lo);
  e.eval({0.31, 0.62, 0.11}, 0.7, hi);
  REQUIRE(lo.U[0] != Approx(hi.U[0]).margin(1e-9));
}

TEST_CASE("BatchelorVortex", "[expr]") {
  Prm p = batchelor_params();
  Expr_BatchelorVortex e(p);
  check_all(e, {0.7, 0.3, 0.2}, {1.9, -0.4, 0.1}, 0.);
  for (double s : {2.0, 0.5, 1e-2, 1e-4})
    check_gradient(e, {s, s/2, 0.3}, 0.);
}

TEST_CASE("BatchelorVortex is regular on the axis", "[expr]") {
  Prm p = batchelor_params();
  Expr_BatchelorVortex e(p);
  const double u0 = 1.0, R1 = 1.0, q = 0.8;

  PointValues pv(1.0);
  e.eval({0., 0., 0.2}, 0., pv);
  for (int i = 0; i < 3; ++i){
    REQUIRE(std::isfinite(pv.U[i]));
    for (int j = 0; j < 3; ++j) REQUIRE(std::isfinite(pv.gradU(i, j)));
  }
  // the core rotates as a solid body at u0/R1
  REQUIRE(pv.U[0] == Approx(0.).margin(1e-12));
  REQUIRE(pv.U[1] == Approx(0.).margin(1e-12));
  REQUIRE(pv.U[2] == Approx(q*u0).margin(1e-12));
  REQUIRE(pv.gradU(0, 1) == Approx(-u0/R1).margin(1e-12));
  REQUIRE(pv.gradU(1, 0) == Approx(u0/R1).margin(1e-12));
}

TEST_CASE("BatchelorVortex has no kink at the series thresholds", "[expr]") {
  Prm p = batchelor_params();
  Expr_BatchelorVortex e(p);
  // phi switches to expm1 at a = 1e-8, dphi at a = 1e-3; R1 = 1 so a = s^2
  for (double a : {1e-3, 1e-8}){
    const double sm = std::sqrt(a*(1. - 1e-9)), sp = std::sqrt(a*(1. + 1e-9));
    PointValues lo(1.0), hi(1.0);
    e.eval({sm/std::sqrt(2.), sm/std::sqrt(2.), 0.2}, 0., lo);
    e.eval({sp/std::sqrt(2.), sp/std::sqrt(2.), 0.2}, 0., hi);
    for (int i = 0; i < 3; ++i){
      for (int j = 0; j < 3; ++j){
        INFO("a = " << a << ", gradU(" << i << ", " << j << ")");
        REQUIRE(hi.gradU(i, j) == Approx(lo.gradU(i, j)).margin(1e-7));
      }
    }
  }
}

TEST_CASE("PlanePoiseuille", "[expr]") {
  Prm p = {{"R", "1.0"}, {"mu", "1.0"}, {"u_inf", "1.0"}, {"p_inf", "0.5"},
           {"rho", "1.2"}, {"x0", "0.0"}, {"y0", "0.0"}, {"z0", "0.0"}};
  Expr_PlanePoiseuille e(p);
  check_all(e, {0.3, 0.2, 0.1}, {-0.4, 0.15, 0.6}, 0.);
}

TEST_CASE("HagenPoiseuille", "[expr]") {
  Prm p = {{"R", "1.0"}, {"mu", "1.0"}, {"u_inf", "1.0"}, {"p_inf", "0.5"},
           {"rho", "1.2"}, {"x0", "0.0"}, {"y0", "0.0"}, {"z0", "0.0"}};
  Expr_HagenPoiseuille e(p);
  check_all(e, {0.3, 0.2, 0.1}, {-0.4, 0.15, 0.6}, 0.);
}

TEST_CASE("InfinitePlate", "[expr]") {
  Prm p = {{"alpha", "0.5"}, {"mu", "1.0"}, {"p_inf", "0.5"}, {"rho", "1.2"},
           {"x0", "0.0"}, {"y0", "0.0"}, {"z0", "0.0"}};
  Expr_InfinitePlane e(p);  // the file is named Plate, the class Plane
  check_all(e, {-0.4, 0.3, 0.1}, {-1.1, -0.2, 0.7}, 0.);
}

TEST_CASE("StokesSphere", "[expr]") {
  Prm p = {{"R", "1.0"}, {"mu", "1.0"}, {"u_inf", "1.0"}, {"p_inf", "0.5"},
           {"rho", "1.2"}, {"x0", "0.0"}, {"y0", "0.0"}, {"z0", "0.0"}};
  Expr_StokesSphere e(p);
  // outside the sphere, where the solution is defined
  check_all(e, {2.1, 1.3, 0.7}, {-1.8, 0.4, 2.2}, 0.);
}

namespace {

Prm brinkman_params(){
  return {{"R", "1.0"}, {"H", "1.0"}, {"mu", "1.0"}, {"u_inf", "1.0"},
          {"p_inf", "1.0"}, {"rho", "1.2"}, {"x0", "0.0"}, {"y0", "0.0"},
          {"z0", "0.0"}, {"rmax_intp", "20.0"}, {"Nintp", "400000"}};
}

}  // namespace

TEST_CASE("BrinkmanCylinder", "[expr]") {
  Prm p = brinkman_params();
  Expr_BrinkmanCylinder e(p);
  // f, f' and f'' come from separate tables, which must resolve a difference
  // outside the cylinder, and away from the interpolation end points
  const Vector3d x{2.1, 1.3, 0.0}, x2{-1.8, 0.4, 0.0};
  check_all(e, x, x2, 0.);
}


namespace {

Prm taylor_couette_params(){
  return {{"R", "2.5"}, {"H", "3.0"}, {"K", "0.95"}, {"a", "0.4"},
          {"c1", "1.8"}, {"c2", "1.2"}, {"c3", "0.17"}, {"rho", "1.2"},
          {"p_inf", "0.5"}, {"x0", "0.0"}, {"y0", "0.0"}, {"z0", "0.0"}};
}

// the flow is given in cylindrical components; recover them from the Cartesian
Vector3d cylindrical(Expr& e, const double s, const double theta, const double z){
  const double c = std::cos(theta), sn = std::sin(theta);
  PointValues pv(1.0);
  e.eval({s*c, s*sn, z}, 0., pv);
  return {pv.U[0]*c + pv.U[1]*sn, -pv.U[0]*sn + pv.U[1]*c, pv.U[2]};
}

}  // namespace

TEST_CASE("TaylorCouette", "[expr]") {
  Prm p = taylor_couette_params();
  Expr_TaylorCouette e(p);
  check_all(e, {1.4, 0.9, 0.3}, {-1.9, 0.7, -0.8}, 0.);
  // across the cell, including near both plates where the corner terms bite
  for (double z : {-1.4, -0.7, 0.0, 0.7, 1.4})
    for (double s : {1.05, 1.5, 2.4}){
      check_gradient(e, {s*0.6, s*0.8, z}, 0.);
      check_incompressible(e, {s*0.6, s*0.8, z}, 0.);
    }
}

TEST_CASE("TaylorCouette holds its boundary conditions", "[expr]") {
  Prm p = taylor_couette_params();
  Expr_TaylorCouette e(p);
  const double R = 2.5, H = 3.0, tol = 1e-12;

  for (double z : {-1.2, -0.5, 0.0, 0.5, 1.2}){
    // the cylinders are impermeable, and the inner one turns at unit speed
    REQUIRE(cylindrical(e, 1.0, 0.9, z)[0] == Approx(0.).margin(tol));
    REQUIRE(cylindrical(e, R, 0.9, z)[0] == Approx(0.).margin(tol));
    REQUIRE(cylindrical(e, 1.0, 0.9, z)[1] == Approx(1.).margin(tol));
    // the fit leaves a small slip on the outer cylinder, from the corner terms
    REQUIRE(std::abs(cylindrical(e, R, 0.9, z)[1]) < 5e-3);
  }
  for (double s : {1.2, 1.8, 2.3}){
    for (double z : {-H/2, H/2}){
      REQUIRE(cylindrical(e, s, 0.4, z)[2] == Approx(0.).margin(tol));  // no flux
      REQUIRE(cylindrical(e, s, 0.4, z)[1] == Approx(0.).margin(tol));  // at rest
    }
  }
}

TEST_CASE("TaylorCouette is finite in the singular corners", "[expr]") {
  // where the turning inner cylinder meets a stationary plate the fit is
  // genuinely discontinuous; it must still not produce a NaN
  Prm p = taylor_couette_params();
  Expr_TaylorCouette e(p);
  for (double z : {-1.5, 1.5}){
    for (double s : {1.0, 1.0 + 1e-14, 2.5}){
      PointValues pv(1.0);
      e.eval({s, 0., z}, 0., pv);
      INFO("s = " << s << ", z = " << z);
      for (int i = 0; i < 3; ++i){
        REQUIRE(std::isfinite(pv.U[i]));
        for (int j = 0; j < 3; ++j) REQUIRE(std::isfinite(pv.gradU(i, j)));
      }
    }
  }
  // and inside the solid inner cylinder, which RK stages can dip into
  for (double s : {0.3, 0.999}){
    PointValues pv(1.0);
    e.eval({s, 0., 0.4}, 0., pv);
    INFO("s = " << s);
    for (int i = 0; i < 3; ++i){
      REQUIRE(std::isfinite(pv.U[i]));
      for (int j = 0; j < 3; ++j) REQUIRE(std::isfinite(pv.gradU(i, j)));
    }
  }
  REQUIRE(e.inside({1.5, 0., 0.}, 0.));
  REQUIRE_FALSE(e.inside({0.5, 0., 0.}, 0.));
  REQUIRE_FALSE(e.inside({1.5, 0., 1.6}, 0.));
}
