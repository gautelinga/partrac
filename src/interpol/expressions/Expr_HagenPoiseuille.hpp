#include "Expr.hpp"

#ifndef __EXPR_HAGENPOISEUILLE_HPP
#define __EXPR_HAGENPOISEUILLE_HPP

//using namespace std;

class Expr_HagenPoiseuille final : public Expr {
public:
  Expr_HagenPoiseuille(const partrac::Params& expr_params) : Expr(expr_params) {
    //R = expr_params.get<double>("R");
    x0 = {expr_params.get<double>("x0"),
          expr_params.get<double>("y0"),
          expr_params.get<double>("z0")};
    u_inf = expr_params.get<double>("u_inf");
    R = expr_params.get<double>("R");
    //alpha = expr_params.get<double>("alpha");
    mu = expr_params.get<double>("mu");
    p_inf = expr_params.get<double>("p_inf");
    Rho = expr_params.get<double>("rho");
  };
  // Keys this expression reads
  static void add_params(partrac::Schema& s) {
    s.require<double>("R", "pipe radius");
    s.require<double>("u_inf", "centreline velocity");
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
    double chi = pow(r[0]/R, 2) + pow(r[1]/R, 2);
    is_inside = chi <= 1.;

    Ux = 0.;
    Uy = 0.;
    Uz = 2*u_inf*(1.0 - chi); // * alpha * r[0] * r[2];
    P = p_inf;

    Uxx = 0.;
    Uxy = 0.;
    Uxz = 0.;
    Uyx = 0.;
    Uyy = 0.;
    Uyz = 0.;
    Uzx = - 4 * u_inf * r[0]/pow(R, 2);
    Uzy = - 4 * u_inf * r[1]/pow(R, 2);
    Uzz = 0.;
  };
  // Distance to the pipe wall, positive in the fluid
  bool has_wall() const { return true; };
  double sdf(const Vector3d &x) const {
    const Vector3d r = x-x0;
    return R - std::hypot(r[0], r[1]);
  };
  Vector3d sdf_grad(const Vector3d &x) const {
    const Vector3d r = x-x0;
    const double s = std::hypot(r[0], r[1]);
    return s > 0. ? Vector3d(-r[0]/s, -r[1]/s, 0.) : Vector3d::Zero();
  };
  bool inside(const Vector3d &x, const double t __attribute__((unused))) {
    Vector3d r = x-x0;
    double chi = pow(r[0]/R, 2) + pow(r[1]/R, 2);
    bool _is_inside = chi <= 1.;
    return _is_inside;
  };
  void eval(const Vector3d &x, const double t __attribute__((unused)), PointValues& ptvals) {
    Vector3d r = x-x0;
    double chi = pow(r[0]/R, 2) + pow(r[1]/R, 2);
    ptvals.U = {0., 0., 2*u_inf*(1.0 - chi)}; // * alpha * r[0] * r[2];
    ptvals.P = p_inf;
    ptvals.Rho = Rho;
    ptvals.gradU << 0., 0., 0.,
                    0., 0., 0.,
                    - 4 * u_inf * r[0]/pow(R, 2),
                    - 4 * u_inf * r[1]/pow(R, 2),
                    0.;
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
  //double alpha;  // alpha parameter = u''(x)
  Vector3d x0;  // origin
  double u_inf;  // Far-field velocity
  double R; //
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
