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
  };
  void eval(const Vector3d &x, const double t __attribute__((unused))) {
    Vector3d r = x-x0;

    double sx = B*sin(2*M_PI*r[0]/L);
    double sy = C*sin(2*M_PI*r[1]/L);
    double sz = A*sin(2*M_PI*r[2]/L);
    double cx = B*cos(2*M_PI*r[0]/L);
    double cy = C*cos(2*M_PI*r[1]/L);
    double cz = A*cos(2*M_PI*r[2]/L);

    Ux = sz + cy;
    Uy = sx + cz;
    Uz = sy + cx;
    P = p_inf + rho_inf*(sz*cy + sx*cz + sy*cx);

    Uxx = 0.;
    Uxy = -sy*2*M_PI/L;
    Uxz = cz*2*M_PI/L;
    Uyx = -sx*2*M_PI/L;
    Uyy = 0.;
    Uyz = cz*2*M_PI/L;
    Uzx = -sx*2*M_PI/L;
    Uzy = cy*2*M_PI/L;
    Uzz = 0.;
  };
  double ux() { return Ux; };
  double uy() { return Uy; };
  double uz() { return Uz; };
  double rho() { return rho_inf; };
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
  Vector3d x0;  // Center of vortex
  double p_inf;  // Far-field pressure
  double rho_inf;
  double A, B, C;
  double L;
  // Useful quantities
  double Ux;
  double Uy;
  double Uz;
  double P;
  double Uxx, Uxy, Uxz;
  double Uyx, Uyy, Uyz;
  double Uzx, Uzy, Uzz;
};

#endif
