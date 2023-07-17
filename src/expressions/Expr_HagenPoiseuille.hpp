#include "Expr.hpp"

#ifndef __EXPR_HAGENPOISEUILLE_HPP
#define __EXPR_HAGENPOISEUILLE_HPP

//using namespace std;

class Expr_HagenPoiseuille : public Expr {
public:
  Expr_HagenPoiseuille(std::map<std::string, std::string> &expr_params) : Expr(expr_params) {
    //R = getd(expr_params, "R");
    x0 = {getd(expr_params, "x0"),
          getd(expr_params, "y0"),
          getd(expr_params, "z0")};
    u_inf = getd(expr_params, "u_inf");
    R = getd(expr_params, "R");
    //alpha = getd(expr_params, "alpha");
    mu = getd(expr_params, "mu");
    p_inf = getd(expr_params, "p_inf");
    Rho = getd(expr_params, "rho");
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
    Uzx = - 2 * u_inf * r[0]/pow(R, 2);
    Uzy = - 2 * u_inf * r[1]/pow(R, 2);
    Uzz = 0.;
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
