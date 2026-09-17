#include "Expr.hpp"

#ifndef __EXPR_LINEARFLOW_HPP
#define __EXPR_LINEARFLOW_HPP

// Linear flow u = A (x - x0), entries Axy = dU_x/dy
class Expr_LinearFlow final : public Expr {
public:
  Expr_LinearFlow(const partrac::Params& expr_params) : Expr(expr_params) {
    x0 = {expr_params.get<double>("x0"),
          expr_params.get<double>("y0"),
          expr_params.get<double>("z0")};
    const char* axis = "xyz";
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        A(i, j) = expr_params.get<double>(std::string("A") + axis[i] + axis[j]);
    p_inf = expr_params.get<double>("p_inf");
    Rho = expr_params.get<double>("rho");
  };
  // Keys this expression reads
  static void add_params(partrac::Schema& s) {
    s.opt<double>("Axx", 0., "velocity gradient dU_x/dx");
    s.opt<double>("Axy", 0., "velocity gradient dU_x/dy");
    s.opt<double>("Axz", 0., "velocity gradient dU_x/dz");
    s.opt<double>("Ayx", 0., "velocity gradient dU_y/dx");
    s.opt<double>("Ayy", 0., "velocity gradient dU_y/dy");
    s.opt<double>("Ayz", 0., "velocity gradient dU_y/dz");
    s.opt<double>("Azx", 0., "velocity gradient dU_z/dx");
    s.opt<double>("Azy", 0., "velocity gradient dU_z/dy");
    s.opt<double>("Azz", 0., "velocity gradient dU_z/dz");
    s.opt<double>("x0", 0., "centre x");
    s.opt<double>("y0", 0., "centre y");
    s.opt<double>("z0", 0., "centre z");
    s.opt<double>("p_inf", 0., "pressure");
    s.opt<double>("rho", 1., "density");
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
