#include "Expr.hpp"
#include <cmath>
#include <functional>

#ifndef __EXPR_BRINKMANCYLINDER_HPP
#define __EXPR_BRINKMANCYLINDER_HPP

//using namespace std;

double beta(double zeta, double r){
  return 2*std::cyl_bessel_k(1, zeta * r)/(zeta * std::cyl_bessel_k(0, zeta));
}

double betar(double zeta, double r){
  return - 2 * std::cyl_bessel_k(0, zeta * r)/std::cyl_bessel_k(0, zeta) - beta(zeta, r)/r;
}

double betarr(double zeta, double r){
  return (pow(zeta, 2) + pow(r, -2))*beta(zeta, r) - betar(zeta, r)/r;
}

double f(double zeta, double r){
  return r - (1+beta(zeta, 1.0))/r + beta(zeta, r);
}

double fr(double zeta, double r){
  return 1 + (1+beta(zeta, 1.0))/pow(r, 2) + betar(zeta, r);
}

double frr(double zeta, double r){
  return -2.*(1 + beta(zeta, 1.0))/pow(r, 3) + betarr(zeta, r);
}

double prf(double zeta, double r){
  return pow(zeta, 2)*(r + (1+beta(zeta, 1.0))/r);
}


class LinIntp {
public:
  LinIntp() {};
  ~LinIntp() {};
  void load(std::function<double(double, double)> func, const double x0, const double x1, const Uint N, const double arg0) {
    y.clear();
    this->N = N;
    this->x0 = x0;
    this->x1 = x1;
    dx = (x1-x0)/N;
    for (Uint i=0; i <= N; ++i){
      double x = x0 + i * dx;
      y.push_back(func(arg0, x));
    }
  };
  double eval(double x){
    if (x <= x0){
      return y[0];
    }
    if (x >= x1){
      return y[N];
    }
    double iest = (x - x0)/dx;
    int i0 = floor(iest);
    int i1 = ceil(iest);
    double alpha = iest-i0;
    return alpha*y[i0] + (1-alpha)*y[i1];
  };
private:
  std::vector<double> y;
  double x0, x1, dx;
  Uint N;
};

class Expr_BrinkmanCylinder : public Expr {
public:
  Expr_BrinkmanCylinder(std::map<std::string, std::string> &expr_params) : Expr(expr_params) {
    R = getd(expr_params, "R");
    H = getd(expr_params, "H");
    mu = getd(expr_params, "mu");
    x0 = {getd(expr_params, "x0"),
          getd(expr_params, "y0"),
          getd(expr_params, "z0")};
    u_inf = getd(expr_params, "u_inf");
    p_inf = getd(expr_params, "p_inf");
    zeta = sqrt(3)*2*R/H;

    double rmax = getd(expr_params, "rmax_intp");
    int N = geti(expr_params, "Nintp");

    f_intp.load(&f, 1., rmax, N, zeta);
    fr_intp.load(&fr, 1., rmax, N, zeta);
    frr_intp.load(&frr, 1., rmax, N, zeta);
    prf_intp.load(&prf, 1., rmax, N, zeta);
  };
  void eval(const Vector3d &x, const double t __attribute__((unused))) {
    //cout << "x = " << x << endl;
    //cout << "x0 = " << x0 << endl;
    //cout << "u_inf = " << u_inf << endl;

    Vector3d s = x-x0;
    s[2] = 0.;
    Vector3d r = s/R;
    double r2 = r.squaredNorm();
    //double R2 = R*R;
    is_inside = r2 >= 1.;

    double rabs = sqrt(r2);

    /*
    Ux = pow(sin(theta), 2)*fr(r) + f(r)*pow(cos(theta), 2)/r;
    Uy = -sin(theta)*cos(theta)*fr(r) + f(r)*sin(theta)*cos(theta)/r;
    // Hardcoded -- copied from consistency-checked Sympy code
    Uxx = (pow(r, 2)*pow(sin(theta), 2)*frr(r) - 3*r*pow(sin(theta), 2)*fr(r) + r*fr(r) + 3*f(r)*pow(sin(theta), 2) - f(r))*cos(theta)/pow(r, 2);
    Uxy = (pow(r, 2)*pow(sin(theta), 2)*frr(r) + 3*r*pow(cos(theta), 2)*fr(r) - 3*f(r)*pow(cos(theta), 2))*sin(theta)/pow(r, 2);
    Uyx = (-pow(r, 2)*pow(cos(theta), 2)*frr(r) + 3*r*pow(cos(theta), 2)*fr(r) - r*fr(r) - 3*f(r)*pow(cos(theta), 2) + f(r))*sin(theta)/pow(r, 2);
    Uyy = (-pow(r, 2)*pow(sin(theta), 2)*frr(r) + 3*r*pow(sin(theta), 2)*fr(r) - r*fr(r) - 3*f(r)*pow(sin(theta), 2) + f(r))*cos(theta)/pow(r, 2);
    P = p_inf - pow(zeta, 2)*cos(theta)*(r + (1+beta(zeta, 1.0))/r);
    */

    double st = r[1]/rabs; // sin(theta);
    double ct = r[0]/rabs; // cos(theta);
    double s2t = st*st;
    double c2t = ct*ct;

    // These are costly!
    double f_r = f_intp.eval(rabs); // f(zeta, rabs);
    double fr_r = fr_intp.eval(rabs); // fr(zeta, rabs);
    double frr_r = frr_intp.eval(rabs); // frr(zeta, rabs);
    double prf_r = prf_intp.eval(rabs); // prf(zeta, rabs);

    Ux = u_inf * s2t*fr_r + f_r*c2t/rabs;
    Uy = u_inf * st*ct*(-fr_r + f_r/rabs);

    // Hardcoded -- copied from consistency-checked Sympy code
    Uxx = u_inf/R * (r2*s2t*frr_r - 3*rabs*s2t*fr_r + rabs*fr_r + 3*f_r*s2t - f_r)*ct/r2;
    Uxy = u_inf/R * (r2*s2t*frr_r + 3*rabs*c2t*fr_r - 3*f_r*c2t)*st/r2;
    Uyx = u_inf/R * (-r2*c2t*frr_r + 3*rabs*c2t*fr_r - rabs*fr_r - 3*f_r*c2t + f_r)*st/r2;
    Uyy = u_inf/R * (-r2*s2t*frr_r + 3*rabs*s2t*fr_r - rabs*fr_r - 3*f_r*s2t + f_r)*ct/r2;

    P = p_inf - mu * u_inf / R * prf_r * ct;
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
  double R;  // Radius of cylinder
  double H;  // Height of cylinder
  double mu;  // Viscosity
  Vector3d x0;  // Center of cylinder
  double u_inf;  // Far-field velocity
  double p_inf;  // Far-field pressure
  double zeta;
  // Useful quantitites
  double Ux;
  double Uy;
  double Uz;
  double P;
  double Uxx, Uxy, Uyx, Uyy;
  double Uxz = 0.;
  double Uyz = 0.;
  double Uzx = 0.;
  double Uzy = 0.;
  double Uzz = 0.;
  LinIntp f_intp, fr_intp, frr_intp, prf_intp;
};

#endif
