#include "Expr.hpp"

#ifndef __EXPR_INFINITEPLANE_HPP
#define __EXPR_INFINITEPLANE_HPP

//using namespace std;

class Expr_InfinitePlane final : public Expr {
public:
  Expr_InfinitePlane(const partrac::Params& expr_params) : Expr(expr_params) {
    //R = expr_params.get<double>("R");
    x0 = {expr_params.get<double>("x0"),
          expr_params.get<double>("y0"),
          expr_params.get<double>("z0")};
    //u_inf = expr_params.get<double>("u_inf");
    alpha = expr_params.get<double>("alpha");
    mu = expr_params.get<double>("mu");
    p_inf = expr_params.get<double>("p_inf");
    Rho = expr_params.get<double>("rho");
  };
  // Keys this expression reads
  static void add_params(partrac::Schema& s) {
    s.require<double>("alpha", "stagnation strain");
    s.require<double>("mu", "viscosity");
    s.require<double>("p_inf", "pressure");
    s.require<double>("rho", "density");
    s.require<double>("x0", "centre x");
    s.require<double>("y0", "centre y");
    s.require<double>("z0", "centre z");
  };
  void eval(const Vector3d &x, const double t __attribute__((unused))) {
    //cout << "x = " << x << endl;
    //cout << "x0 = " << x0 << endl;
    //cout << "u_inf = " << u_inf << endl;

    Vector3d r = x-x0;
    is_inside = r[0] <= 0.;

    Ux = alpha * r[0] * r[0];
    Uy = - 2 * alpha * r[0] * r[1];
    Uz = 0; // * alpha * r[0] * r[2];
    P = p_inf + 2.0 * mu * alpha * r[0];

    // Hardcoded -- copied from consistency-checked Sympy code
    Uxx = 2 * alpha * r[0];
    Uxy = 0.;
    Uxz = 0.;
    Uyx = - 2 * alpha * r[1];
    Uyy = - 2 * alpha * r[0];
    Uyz = 0.;
    Uzx = 0.;
    Uzy = 0.;
    Uzz = 0.;
  };
  // Distance to the plate, positive in the fluid
  bool has_wall() const { return true; };
  double sdf(const Vector3d &x) const { return x0[0] - x[0]; };
  Vector3d sdf_grad(const Vector3d &x __attribute__((unused))) const { return {-1., 0., 0.}; };
  bool inside(const Vector3d &x, const double t __attribute__((unused))) {
    Vector3d r = x-x0;
    bool _is_inside = r[0] <= 0.;
    return _is_inside;
  };
  void eval(const Vector3d &x, const double t __attribute__((unused)), PointValues& ptvals) {
    Vector3d r = x-x0;

    ptvals.U = {alpha * r[0] * r[0], - 2 * alpha * r[0] * r[1], 0};  // * alpha * r[0] * r[2];
    ptvals.P = p_inf + 2.0 * mu * alpha * r[0];
    ptvals.Rho = Rho;

    // Hardcoded -- copied from consistency-checked Sympy code
    ptvals.gradU << 2 * alpha * r[0], 0., 0.,
                    - 2 * alpha * r[1], - 2 * alpha * r[0], 0.,
                    0., 0., 0.;
  };
  double ux() { return Ux; };
  double uy() { return Uy; };
  double uz() { return Uz; };
  double rho() { return Rho; };
  double p() { return P; };
  double uxx() { return Uxx; };
  double uxy() { return Uxy; };
  double uxz() { return Uxz; };
  double uyx() { return Uyx; };
  double uyy() { return Uyy; };
  double uyz() { return Uyz; };
  double uzx() { return Uzx; };
  double uzy() { return Uzy; };
  double uzz() { return Uzz; };
private:
  double mu;  // Viscosity
  double alpha;  // alpha parameter = u''(x)
  Vector3d x0;  // origin
  //double u_inf;  // Far-field velocity
  double p_inf;  // Far-field pressure
  double Rho;
  // Useful quantitites
  double Ux;
  double Uy;
  double Uz;
  double P;
  double Uxx, Uxy, Uxz;
  double Uyx, Uyy, Uyz;
  double Uzx, Uzy, Uzz;
};

#endif
