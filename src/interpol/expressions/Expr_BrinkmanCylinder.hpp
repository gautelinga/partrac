#include "Expr.hpp"
#include <cmath>
#include <vector>

#ifndef __EXPR_BRINKMANCYLINDER_HPP
#define __EXPR_BRINKMANCYLINDER_HPP

// Radial profile from K0, K1 at zeta r
struct BrinkmanProfile {
  double zeta, k0z, b1;
  double beta(double k1) const {
    return 2*k1/(zeta * k0z);
  }
  double betar(double r, double k0, double k1) const {
    return - 2 * k0/k0z - beta(k1)/r;
  }
  double betarr(double r, double k0, double k1) const {
    return (pow(zeta, 2) + pow(r, -2))*beta(k1) - betar(r, k0, k1)/r;
  }
  double f(double r, double k1) const {
    return r - (1+b1)/r + beta(k1);
  }
  double fr(double r, double k0, double k1) const {
    return 1 + (1+b1)/pow(r, 2) + betar(r, k0, k1);
  }
  double frr(double r, double k0, double k1) const {
    return -2.*(1 + b1)/pow(r, 3) + betarr(r, k0, k1);
  }
  double prf(double r) const {
    return pow(zeta, 2)*(r + (1+b1)/r);
  }
};


class LinIntp {
public:
  LinIntp() {};
  ~LinIntp() {};
  // N + 1 equispaced values on [x0, x1]
  void load(std::vector<double> values, const double x0, const double x1) {
    y = std::move(values);
    N = Uint(y.size() - 1);
    this->x0 = x0;
    this->x1 = x1;
    dx = (x1-x0)/N;
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

    // Tables from K0 and K1, once a point
    const double k0z = std::cyl_bessel_k(0.0, zeta);
    const BrinkmanProfile pr{zeta, k0z, 2*std::cyl_bessel_k(1.0, zeta)/(zeta * k0z)};
    const double dr = (rmax - 1.)/Uint(N);
    std::vector<double> f(N + 1), fr(N + 1), frr(N + 1), prf(N + 1);
#pragma omp parallel for schedule(static)
    for (int i = 0; i <= N; ++i){
      const double r = 1. + Uint(i) * dr;
      const double k0 = std::cyl_bessel_k(0.0, zeta * r), k1 = std::cyl_bessel_k(1.0, zeta * r);
      f[i] = pr.f(r, k1);
      fr[i] = pr.fr(r, k0, k1);
      frr[i] = pr.frr(r, k0, k1);
      prf[i] = pr.prf(r);
    }
    f_intp.load(std::move(f), 1., rmax);
    fr_intp.load(std::move(fr), 1., rmax);
    frr_intp.load(std::move(frr), 1., rmax);
    prf_intp.load(std::move(prf), 1., rmax);
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
    double f_r = f_intp.eval(rabs);
    double fr_r = fr_intp.eval(rabs);
    double frr_r = frr_intp.eval(rabs);
    double prf_r = prf_intp.eval(rabs);

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
private:
  double R;  // Radius of cylinder
  double H;  // Height of cylinder
  double mu;  // Viscosity
  double Rho;
  Vector3d x0;  // Center of cylinder
  double u_inf;  // Far-field velocity
  double p_inf;  // Far-field pressure
  double zeta;
  LinIntp f_intp, fr_intp, frr_intp, prf_intp;
};

#endif
