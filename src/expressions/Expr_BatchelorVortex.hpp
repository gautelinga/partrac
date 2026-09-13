#include "Expr.hpp"

#ifndef __EXPR_BATCHELORVORTEX_HPP
#define __EXPR_BATCHELORVORTEX_HPP

//using namespace std;

class Expr_BatchelorVortex final : public Expr {
public:
  Expr_BatchelorVortex(std::map<std::string, std::string> &expr_params) : Expr(expr_params) {
    R1 = getd(expr_params, "R1");
    R2 = getd(expr_params, "R2");
    R12 = R1*R1;
    R22 = R2*R2;
    q = getd(expr_params, "q");
    rho_inf = getd(expr_params, "rho");
    x0 = {getd(expr_params, "x0"),
          getd(expr_params, "y0"),
          getd(expr_params, "z0")};
    u0 = getd(expr_params, "u0");
    p_inf = getd(expr_params, "p_inf");
  };
  void eval(const Vector3d &x, const double t __attribute__((unused))) {
    is_inside = true;
    compute(x, U, gradU, P);
  };
  bool inside(const Vector3d &x __attribute__((unused)), const double t __attribute__((unused))) {
    return true;
  };
  void eval(const Vector3d &x, const double t __attribute__((unused)), PointValues& ptvals) {
    compute(x, ptvals.U, ptvals.gradU, ptvals.P);
    ptvals.Rho = rho_inf;
  };
  double ux() { return U[0]; };
  double uy() { return U[1]; };
  double uz() { return U[2]; };
  double rho() { return rho_inf; };
  double p() { return P; };
  double uxx() { return gradU(0, 0); };
  double uxy() { return gradU(0, 1); };
  double uxz() { return gradU(0, 2); };
  double uyx() { return gradU(1, 0); };
  double uyy() { return gradU(1, 1); };
  double uyz() { return gradU(1, 2); };
  double uzx() { return gradU(2, 0); };
  double uzy() { return gradU(2, 1); };
  double uzz() { return gradU(2, 2); };
private:
  // phi(a) = (1 - exp(-a))/a and its derivative, both regular at a = 0
  static double phi(const double a) {
    if (a < 1e-8) return 1. - a/2. + a*a/6.;
    return -std::expm1(-a)/a;
  };
  static double dphi(const double a) {
    if (a < 1e-3) return -0.5 + a/3. - a*a/8.;
    return (a*std::exp(-a) + std::expm1(-a))/(a*a);
  };
  // Writes no members: the PointValues overload runs inside an omp for.
  // phi(s^2/R1^2) absorbs the 1/s^2, so the axis is regular: solid-body at u0/R1.
  void compute(const Vector3d &x, Vector3d& U_, Matrix3d& gradU_, double& P_) const {
    Vector3d r = x-x0;
    double s2 = r[0]*r[0]+r[1]*r[1];
    double a = s2/R12;
    double ph = phi(a);
    double dph = dphi(a);
    double eta22 = exp(-s2/R22);

    U_ = {-u0*r[1]*ph/R1,
           u0*r[0]*ph/R1,
           q*u0*eta22};
    P_ = p_inf;

    double c = 2*u0*dph/(R1*R12);
    // gradU_(i, j) is dU_i/dx_j
    gradU_ <<
      -c*r[0]*r[1],                        -u0*(ph + 2*r[1]*r[1]*dph/R12)/R1, 0.,
       u0*(ph + 2*r[0]*r[0]*dph/R12)/R1,    c*r[0]*r[1],                      0.,
      -2*r[0]*q*u0*eta22/R22,              -2*r[1]*q*u0*eta22/R22,            0.;
  };
  double R1;  // Radius of sphere
  double R2;
  double R12;
  double R22;
  Vector3d x0;  // Center of vortex
  double u0;  // Far-field velocity
  double q; // Amplification
  double p_inf;  // Far-field pressure
  double rho_inf;
  // Useful quantities
  Vector3d U;
  Matrix3d gradU;
  double P;
};

#endif
