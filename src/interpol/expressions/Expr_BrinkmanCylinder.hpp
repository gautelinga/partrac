#include "Expr.hpp"
#include <cmath>
#include <functional>

#ifndef __EXPR_BRINKMANCYLINDER_HPP
#define __EXPR_BRINKMANCYLINDER_HPP

//using namespace std;

inline double beta(double zeta, double r){
  return 2*std::cyl_bessel_k(1.0, zeta * r)/(zeta * std::cyl_bessel_k(0.0, zeta));
}

inline double betar(double zeta, double r){
  return - 2 * std::cyl_bessel_k(0.0, zeta * r)/std::cyl_bessel_k(0.0, zeta) - beta(zeta, r)/r;
}

inline double betarr(double zeta, double r){
  return (pow(zeta, 2) + pow(r, -2))*beta(zeta, r) - betar(zeta, r)/r;
}

inline double f(double zeta, double r){
  return r - (1+beta(zeta, 1.0))/r + beta(zeta, r);
}

inline double fr(double zeta, double r){
  return 1 + (1+beta(zeta, 1.0))/pow(r, 2) + betar(zeta, r);
}

inline double frr(double zeta, double r){
  return -2.*(1 + beta(zeta, 1.0))/pow(r, 3) + betarr(zeta, r);
}

inline double prf(double zeta, double r){
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
    return (1-alpha)*y[i0] + alpha*y[i1];
  };
private:
  std::vector<double> y;
  double x0, x1, dx;
  Uint N;
};

class Expr_BrinkmanCylinder final : public Expr {
public:
  Expr_BrinkmanCylinder(const partrac::Params& expr_params) : Expr(expr_params) {
    R = expr_params.get<double>("R");
    H = expr_params.get<double>("H");
    mu = expr_params.get<double>("mu");
    Rho = expr_params.get<double>("rho");
    x0 = {expr_params.get<double>("x0"),
          expr_params.get<double>("y0"),
          expr_params.get<double>("z0")};
    u_inf = expr_params.get<double>("u_inf");
    p_inf = expr_params.get<double>("p_inf");
    zeta = sqrt(3)*2*R/H;

    double rmax = expr_params.get<double>("rmax_intp");
    int N = expr_params.get<int>("Nintp");

    f_intp.load(&f, 1., rmax, N, zeta);
    fr_intp.load(&fr, 1., rmax, N, zeta);
    frr_intp.load(&frr, 1., rmax, N, zeta);
    prf_intp.load(&prf, 1., rmax, N, zeta);
  };
  // Keys this expression reads
  static void add_params(partrac::Schema& s) {
    s.require<double>("R", "radius");
    s.require<double>("H", "height of the cylinder");
    s.require<double>("mu", "viscosity");
    s.require<double>("u_inf", "far-field velocity");
    s.require<double>("rmax_intp", "outer radius of the tabulated profile");
    s.require<int>("Nintp", "points in the tabulated profile");
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

    Vector3d s = x-x0;
    s[2] = 0.;
    Vector3d r = s/R;
    double r2 = r.squaredNorm();
    //double R2 = R*R;
    is_inside = r2 >= 1.;

    double rabs = sqrt(r2);


    double st = r[1]/rabs; // sin(theta);
    double ct = r[0]/rabs; // cos(theta);
    double s2t = st*st;
    double c2t = ct*ct;

    // These are costly!
    double f_r = f_intp.eval(rabs); // f(zeta, rabs);
    double fr_r = fr_intp.eval(rabs); // fr(zeta, rabs);
    double frr_r = frr_intp.eval(rabs); // frr(zeta, rabs);
    double prf_r = prf_intp.eval(rabs); // prf(zeta, rabs);

    Ux = u_inf * (s2t*fr_r + f_r*c2t/rabs);
    Uy = u_inf * st*ct*(-fr_r + f_r/rabs);

    // Hardcoded -- copied from consistency-checked Sympy code
    Uxx = u_inf/R * (r2*s2t*frr_r - 3*rabs*s2t*fr_r + rabs*fr_r + 3*f_r*s2t - f_r)*ct/r2;
    Uxy = u_inf/R * (r2*s2t*frr_r + 3*rabs*c2t*fr_r - 3*f_r*c2t)*st/r2;
    Uyx = u_inf/R * (-r2*c2t*frr_r + 3*rabs*c2t*fr_r - rabs*fr_r - 3*f_r*c2t + f_r)*st/r2;
    Uyy = u_inf/R * (-r2*s2t*frr_r + 3*rabs*s2t*fr_r - rabs*fr_r - 3*f_r*s2t + f_r)*ct/r2;

    P = p_inf - mu * u_inf / R * prf_r * ct;
  };
  // Distance to the cylinder, positive in the fluid
  bool has_wall() const { return true; };
  double sdf(const Vector3d &x) const { return std::hypot(x[0]-x0[0], x[1]-x0[1]) - R; };
  Vector3d sdf_grad(const Vector3d &x) const {
    const Vector3d r(x[0]-x0[0], x[1]-x0[1], 0.);
    return r.normalized();
  };
  bool inside(const Vector3d &x, const double t __attribute__((unused))) {
    Vector3d s = x-x0;
    s[2] = 0.;
    Vector3d r = s/R;
    double r2 = r.squaredNorm();
    bool _is_inside = r2 >= 1.;
    return _is_inside;
  };
  void eval(const Vector3d &x, const double t __attribute__((unused)), PointValues& ptvals) {
    Vector3d s = x-x0;
    s[2] = 0.;
    Vector3d r = s/R;
    double r2 = r.squaredNorm();
    
    double rabs = sqrt(r2);

    double st = r[1]/rabs; // sin(theta);
    double ct = r[0]/rabs; // cos(theta);
    double s2t = st*st;
    double c2t = ct*ct;

    // These are costly!
    double f_r = f_intp.eval(rabs); // f(zeta, rabs);
    double fr_r = fr_intp.eval(rabs); // fr(zeta, rabs);
    double frr_r = frr_intp.eval(rabs); // frr(zeta, rabs);
    double prf_r = prf_intp.eval(rabs); // prf(zeta, rabs);

    ptvals.U = {u_inf * (s2t*fr_r + f_r*c2t/rabs),
                u_inf * st*ct*(-fr_r + f_r/rabs),
                0.};

    // Hardcoded -- copied from consistency-checked Sympy code
    ptvals.gradU << 
      u_inf/R * (r2*s2t*frr_r - 3*rabs*s2t*fr_r + rabs*fr_r + 3*f_r*s2t - f_r)*ct/r2,
      u_inf/R * (r2*s2t*frr_r + 3*rabs*c2t*fr_r - 3*f_r*c2t)*st/r2,
      0.0,
      u_inf/R * (-r2*c2t*frr_r + 3*rabs*c2t*fr_r - rabs*fr_r - 3*f_r*c2t + f_r)*st/r2,
      u_inf/R * (-r2*s2t*frr_r + 3*rabs*s2t*fr_r - rabs*fr_r - 3*f_r*s2t + f_r)*ct/r2,
      0.0,
      0.0,
      0.0,
      0.0;
    ptvals.P = p_inf - mu * u_inf / R * prf_r * ct;
    ptvals.Rho = Rho;
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
  double R;  // Radius of cylinder
  double H;  // Height of cylinder
  double mu;  // Viscosity
  double Rho;
  Vector3d x0;  // Center of cylinder
  double u_inf;  // Far-field velocity
  double p_inf;  // Far-field pressure
  double zeta;
  // Useful quantitites
  double Ux;
  double Uy;
  double Uz = 0.;  // the flow is in the plane
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
