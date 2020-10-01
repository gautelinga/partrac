#include "Expr.hpp"

#ifndef __EXPR_BATCHELORVORTEX_HPP
#define __EXPR_BATCHELORVORTEX_HPP

//using namespace std;

class Expr_BatchelorVortex : public Expr {
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
    Vector3d r = x-x0;
    double s2 = r[0]*r[0]+r[1]*r[1];
    double s4 = s2*s2;

    double eta12 = exp(-s2/R12);
    double eta22 = exp(-s2/R22);

    Ux = -R1*r[1]*u0*(1.0 - eta12)/s2;
    Uy =  R1*r[0]*u0*(1.0 - eta12)/s2;
    Uz = q*u0*eta22;
    P = p_inf;

    Uxx = 2*r[0]*r[1]*u0*(R12*(1.0 - eta12) - s2*eta12)/(R1*s4);
    Uxy = u0*(2*R12*pow(r[1], 2)*(1.0 - eta12) - R12*(1.0 - eta12)*s2 - 2*pow(r[1], 2)*s2*eta12)/(R1*s4);
    Uxz = 0;
    Uyx = u0*(-2*R12*pow(r[0], 2)*(1.0 - eta12) + R12*(1.0 - eta12)*s2 + 2*pow(r[0], 2)*s2*eta12)/(R1*s4);
    Uyy = -Uxx;
    Uyz = 0;
    Uzx = -2*r[0]*q*u0*eta22/R22;
    Uzy = -2*r[1]*q*u0*eta22/R22;
    Uzz = 0;
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
  double Ux;
  double Uy;
  double Uz;
  double P;
  double Uxx, Uxy, Uxz;
  double Uyx, Uyy, Uyz;
  double Uzx, Uzy, Uzz;
};

#endif
