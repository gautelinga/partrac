#include "Expr.hpp"

#ifndef __EXPR_STOKESSPHERE_HPP
#define __EXPR_STOKESSPHERE_HPP

//using namespace std;

class Expr_StokesSphere : public Expr {
public:
  Expr_StokesSphere(std::map<std::string, std::string> &expr_params) : Expr(expr_params) {
    R = getd(expr_params, "R");
    mu = getd(expr_params, "mu");
    x0 = {getd(expr_params, "x0"),
          getd(expr_params, "y0"),
          getd(expr_params, "z0")};
    u_inf = getd(expr_params, "u_inf");
    p_inf = getd(expr_params, "p_inf");
  };
  void eval(const Vector3d &x, const double t __attribute__((unused))) {
    //cout << "x = " << x << endl;
    //cout << "x0 = " << x0 << endl;
    //cout << "u_inf = " << u_inf << endl;

    Vector3d r = x-x0;
    double r2 = r.squaredNorm();
    double R2 = R*R;
    is_inside = r2 >= R2;

    double eta = R/sqrt(r2);
    double eta3 = pow(eta, 3);

    double f_r = 1.0 + 1./2. * eta3 - 3./2. * eta;
    double f_theta = -1.0 + 1./4. * eta3 + 3./4. * eta;

    double r7 = pow(r2, 7.0/2.0);

    Ux = r[0]*r[2]/r2 * (f_r + f_theta) * u_inf;
    Uy = r[1]*r[2]/r2 * (f_r + f_theta) * u_inf;
    Uz = r[2]*r[2]/r2 * (f_r + f_theta) * u_inf - f_theta * u_inf;
    P = p_inf - 3.0/2.0 * mu * u_inf * r[2] * eta / r2;

    // Hardcoded -- copied from consistency-checked Sympy code
    Uxx = (3.0/4.0)*R*r[2]*u_inf*(-2*pow(r[0], 2)*(R2 - r2) - pow(r[0], 2)*(3*R2 - r2) + r2*(R2 - r2))/r7;
    Uxy = (3.0/4.0)*R*r[0]*r[1]*r[2]*u_inf*(-5*R2 + 3*r2)/r7;
    Uxz = (3.0/4.0)*R*r[0]*u_inf*(-2*pow(r[2], 2)*(R2 - r2) - pow(r[2], 2)*(3*R2 - r2) + r2*(R2 - r2))/r7;
    Uyx = (3.0/4.0)*R*r[0]*r[1]*r[2]*u_inf*(-5*R2 + 3*r2)/r7;
    Uyy = (3.0/4.0)*R*r[2]*u_inf*(-2*pow(r[1], 2)*(R2 - r2) - pow(r[1], 2)*(3*R2 - r2) + r2*(R2 - r2))/r7;
    Uyz = (3.0/4.0)*R*r[1]*u_inf*(-2*pow(r[2], 2)*(R2 - r2) - pow(r[2], 2)*(3*R2 - r2) + r2*(R2 - r2))/r7;
    Uzx = (3.0/4.0)*R*r[0]*u_inf*(-2*pow(r[2], 2)*(R2 - r2) - pow(r[2], 2)*(3*R2 - r2) + r2*(R2 + r2))/r7;
    Uzy = (3.0/4.0)*R*r[1]*u_inf*(-2*pow(r[2], 2)*(R2 - r2) - pow(r[2], 2)*(3*R2 - r2) + r2*(R2 + r2))/r7;
    Uzz = (3.0/4.0)*R*r[2]*u_inf*(-2*pow(r[2], 2)*(R2 - r2) - pow(r[2], 2)*(3*R2 - r2) + 2*r2*(R2 - r2) + r2*(R2 + r2))/r7;
  };
  bool inside(const Vector3d &x, const double t __attribute__((unused))) {
    Vector3d r = x-x0;
    double r2 = r.squaredNorm();
    double R2 = R*R;
    bool _is_inside = r2 >= R2;
    return _is_inside;
  };
  void eval(const Vector3d &x, const double t __attribute__((unused)), PointValues& ptvals) {

    Vector3d r = x-x0;
    double r2 = r.squaredNorm();
    double R2 = R*R;
    // is_inside = r2 >= R2;

    double eta = R/sqrt(r2);
    double eta3 = pow(eta, 3);

    double f_r = 1.0 + 1./2. * eta3 - 3./2. * eta;
    double f_theta = -1.0 + 1./4. * eta3 + 3./4. * eta;

    double r7 = pow(r2, 7.0/2.0);

    ptvals.U = {r[0]*r[2]/r2 * (f_r + f_theta) * u_inf,
                r[1]*r[2]/r2 * (f_r + f_theta) * u_inf,
                r[2]*r[2]/r2 * (f_r + f_theta) * u_inf - f_theta * u_inf};
    ptvals.P = p_inf - 3.0/2.0 * mu * u_inf * r[2] * eta / r2;

    // Hardcoded -- copied from consistency-checked Sympy code
    ptvals.gradU << (3.0/4.0)*R*r[2]*u_inf*(-2*pow(r[0], 2)*(R2 - r2) - pow(r[0], 2)*(3*R2 - r2) + r2*(R2 - r2))/r7,
                    (3.0/4.0)*R*r[0]*r[1]*r[2]*u_inf*(-5*R2 + 3*r2)/r7,
                    (3.0/4.0)*R*r[0]*u_inf*(-2*pow(r[2], 2)*(R2 - r2) - pow(r[2], 2)*(3*R2 - r2) + r2*(R2 - r2))/r7,
                    (3.0/4.0)*R*r[0]*r[1]*r[2]*u_inf*(-5*R2 + 3*r2)/r7,
                    (3.0/4.0)*R*r[2]*u_inf*(-2*pow(r[1], 2)*(R2 - r2) - pow(r[1], 2)*(3*R2 - r2) + r2*(R2 - r2))/r7,
                    (3.0/4.0)*R*r[1]*u_inf*(-2*pow(r[2], 2)*(R2 - r2) - pow(r[2], 2)*(3*R2 - r2) + r2*(R2 - r2))/r7,
                    (3.0/4.0)*R*r[0]*u_inf*(-2*pow(r[2], 2)*(R2 - r2) - pow(r[2], 2)*(3*R2 - r2) + r2*(R2 + r2))/r7,
                    (3.0/4.0)*R*r[1]*u_inf*(-2*pow(r[2], 2)*(R2 - r2) - pow(r[2], 2)*(3*R2 - r2) + r2*(R2 + r2))/r7,
                    (3.0/4.0)*R*r[2]*u_inf*(-2*pow(r[2], 2)*(R2 - r2) - pow(r[2], 2)*(3*R2 - r2) + 2*r2*(R2 - r2) + r2*(R2 + r2))/r7;
  };
  double ux() { return Ux; };
  double uy() { return Uy; };
  double uz() { return Uz; };
  double rho() { return getd(expr_params, "rho"); };
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
  double R;  // Radius of sphere
  double mu;  // Viscosity
  Vector3d x0;  // Center of sphere
  double u_inf;  // Far-field velocity
  double p_inf;  // Far-field pressure
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
