#include "Expr.hpp"

#ifndef __EXPR_LINEARFLOW_HPP
#define __EXPR_LINEARFLOW_HPP

// Linear flow u = A (x - x0), entries Axy = dU_x/dy
class Expr_LinearFlow final : public Expr {
public:
  Expr_LinearFlow(std::map<std::string, std::string> &expr_params) : Expr(expr_params) {
    x0 = {getd(expr_params, "x0", 0.),
          getd(expr_params, "y0", 0.),
          getd(expr_params, "z0", 0.)};
    const char* axis = "xyz";
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        A(i, j) = getd(expr_params, std::string("A") + axis[i] + axis[j], 0.);
    p_inf = getd(expr_params, "p_inf", 0.);
    Rho = getd(expr_params, "rho", 1.);
  };
  void eval(const Vector3d &x, const double t __attribute__((unused))) {
    is_inside = true;
    U = A * (x - x0);
  };
  bool inside(const Vector3d &x __attribute__((unused)), const double t __attribute__((unused))) {
    return true;
  };
  void eval(const Vector3d &x, const double t __attribute__((unused)), PointValues& ptvals) {
    ptvals.U = A * (x - x0);
    ptvals.P = p_inf;
    ptvals.Rho = Rho;
    ptvals.gradU = A;
  };
  double ux() { return U[0]; };
  double uy() { return U[1]; };
  double uz() { return U[2]; };
  double rho() { return Rho; };
  double p() { return p_inf; };
  double uxx() { return A(0, 0); };
  double uxy() { return A(0, 1); };
  double uxz() { return A(0, 2); };
  double uyx() { return A(1, 0); };
  double uyy() { return A(1, 1); };
  double uyz() { return A(1, 2); };
  double uzx() { return A(2, 0); };
  double uzy() { return A(2, 1); };
  double uzz() { return A(2, 2); };
private:
  Matrix3d A;
  Vector3d x0;
  Vector3d U;
  double p_inf;
  double Rho;
};

#endif
