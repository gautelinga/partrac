#include "Expr.hpp"

#ifndef __EXPR_ABCFlOW_HPP
#define __EXPR_ABCFLOW_HPP

//using namespace std;

class Expr_ABCFlow : public Expr {
public:
  Expr_ABCFlow(std::map<std::string, std::string> &expr_params) : Expr(expr_params) {
    A = getd(expr_params, "A");
    B = getd(expr_params, "B");
    C = getd(expr_params, "C");
    L = getd(expr_params, "L");
    rho_inf = getd(expr_params, "rho");
    x0 = {getd(expr_params, "x0"),
          getd(expr_params, "y0"),
          getd(expr_params, "z0")};
    p_inf = getd(expr_params, "p_inf");
    // Optional forcing of A; zero leaves the steady flow
    A_amp = getd(expr_params, "A_amp", 0.);
    A_omega = getd(expr_params, "A_omega", 0.);
  };
  void eval(const Vector3d &x, const double t) {
    is_inside = true;
    compute(x, t, U, gradU, P);
  };
  bool inside(const Vector3d &x __attribute__((unused)), const double t __attribute__((unused))) {
    return true;
  };
  void eval(const Vector3d &x, const double t, PointValues& ptvals) {
    compute(x, t, ptvals.U, ptvals.gradU, ptvals.P);
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
  // Writes no members: the PointValues overload runs inside an omp for
  void compute(const Vector3d &x, const double t, Vector3d& U_, Matrix3d& gradU_, double& P_) const {
    Vector3d r = x-x0;
    double k = 2*M_PI/L;
    double At = A + A_amp*sin(A_omega*t);

    double sx = B*sin(k*r[0]);
    double sy = C*sin(k*r[1]);
    double sz = At*sin(k*r[2]);
    double cx = B*cos(k*r[0]);
    double cy = C*cos(k*r[1]);
    double cz = At*cos(k*r[2]);

    U_ = {sz + cy, sx + cz, sy + cx};
    // The Bernoulli sum: the pressure only while A_amp is zero
    P_ = p_inf - rho_inf*(sz*cy + sx*cz + sy*cx);

    // gradU_(i, j) is dU_i/dx_j
    gradU_ <<
       0.,     -k*sy,   k*cz,
       k*cx,    0.,    -k*sz,
      -k*sx,    k*cy,   0.;
  };
  Vector3d x0;  // Center of vortex
  double p_inf;  // Far-field pressure
  double rho_inf;
  double A, B, C;
  double A_amp, A_omega;  // forcing of A
  double L;
  // Useful quantities
  Vector3d U;
  Matrix3d gradU;
  double P;
};

#endif
