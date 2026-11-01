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
  double betarrr(double r, double k0, double k1) const {
    return pow(zeta, 2)*betar(r, k0, k1) - betarr(r, k0, k1)/r + 2*betar(r, k0, k1)/pow(r, 2)
         - 2*beta(k1)/pow(r, 3);
  }
  double frrr(double r, double k0, double k1) const {
    return 6.*(1 + b1)/pow(r, 4) + betarrr(r, k0, k1);
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
  __attribute__((always_inline))
  double eval(double x) const {
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
    std::vector<double> f(N + 1), fr(N + 1), frr(N + 1), frrr(N + 1), prf(N + 1);
#pragma omp parallel for schedule(static)
    for (int i = 0; i <= N; ++i){
      const double r = 1. + Uint(i) * dr;
      const double k0 = std::cyl_bessel_k(0.0, zeta * r), k1 = std::cyl_bessel_k(1.0, zeta * r);
      f[i] = pr.f(r, k1);
      fr[i] = pr.fr(r, k0, k1);
      frr[i] = pr.frr(r, k0, k1);
      frrr[i] = pr.frrr(r, k0, k1);
      prf[i] = pr.prf(r);
    }
    f_intp.load(std::move(f), 1., rmax);
    fr_intp.load(std::move(fr), 1., rmax);
    frr_intp.load(std::move(frr), 1., rmax);
    frrr_intp.load(std::move(frrr), 1., rmax);
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
    if (convective)
      convective_terms(s[0], s[1], rabs, f_r, fr_r, frr_r, ptvals.convGradU);
    ptvals.P = p_inf - mu * u_inf / R * prf_r * ct;
    ptvals.Rho = Rho;
  };
  // (U.grad) gradU, out of the hot eval
  __attribute__((noinline))
  void convective_terms(const double r_0, const double r_1, const double rabs, const double f_r, const double fr_r, const double frr_r, Matrix3d& out) const {
    const double frrr_r = frrr_intp.eval(rabs);
    // (U.grad) gradU, generated by Sympy (python/codegen/conv_gradu.py)
    const double w0 = ((u_inf)*(u_inf));
    const double w1 = ((r_0)*(r_0));
    const double w2 = ((r_1)*(r_1));
    const double w3 = w1 + w2;
    const double w4 = (1.0/((w3)*(w3)));
    const double w5 = fr_r*w4;
    const double w6 = 2*w5;
    const double w7 = (1.0/((w3)*(w3)*sqrt(w3)));
    const double w8 = R*f_r;
    const double w9 = w7*w8;
    const double w10 = (1.0/((w3)*sqrt(w3)));
    const double w11 = (1.0/(R));
    const double w12 = frr_r*w11;
    const double w13 = w10*w12;
    const double w14 = 2*w13;
    const double w15 = fr_r/((w3)*(w3)*(w3));
    const double w16 = 7*w15;
    const double w17 = 8*w15;
    const double w18 = 15*w1;
    const double w19 = w8/((w3)*(w3)*(w3)*sqrt(w3));
    const double w20 = (1.0/((R)*(R)));
    const double w21 = frrr_r*w20;
    const double w22 = w21*w4;
    const double w23 = w12*w7;
    const double w24 = 5*w23;
    const double w25 = -w1*w16 + w1*w23 + w14 + w17*w2 + w18*w19 + w2*w22 - w2*w24 - w6 - 6*w9;
    const double w26 = (1.0/(w3));
    const double w27 = (1.0/sqrt(w3));
    const double w28 = -fr_r + w27*w8;
    const double w29 = w26*w28;
    const double w30 = w1*w2;
    const double w31 = w29*w30;
    const double w32 = fr_r*w26;
    const double w33 = w10*w8;
    const double w34 = w1*w33 + w2*w32;
    const double w35 = w1*w5;
    const double w36 = ((r_0)*(r_0)*(r_0)*(r_0));
    const double w37 = 15*w19;
    const double w38 = w13*w2;
    const double w39 = ((r_1)*(r_1)*(r_1)*(r_1));
    const double w40 = 3*w9;
    const double w41 = w1*w40;
    const double w42 = r_0*r_1;
    const double w43 = w12*w27 - w32 + w33;
    const double w44 = w1*w43;
    const double w45 = 4*w26;
    const double w46 = 8*w28*w4;
    const double w47 = -3*R*f_r*w10 - 3*frr_r*w11*w27 - 6*w26*w28 + 3*w32;
    const double w48 = 2*w29;
    const double w49 = w2*w43;
    const double w50 = R*f_r*w27 - fr_r + 4*w1*w2*w26*w43 + 8*w1*w2*w28*w4 + w1*w2*(w14 - w21*w26 + w40 - 3*w5) - w1*w48 - w2*w48 - w44 - w49;
    const double w51 = w0*w26;
    out <<
      w0*(w25*w31 + w34*(-w16*w36 + w17*w30 - w18*w9 - w2*w6 + w22*w30 + w23*w36 - w24*w30 + 2*w33 + 5*w35 + w36*w37 + w38)),
      w0*w42*(w25*w34 + w29*(-w16*w30 + w17*w39 - 10*w2*w5 + w22*w39 + w23*w30 - w24*w39 + w30*w37 + 2*w32 + w35 + 5*w38 - w41)),
      0,
      w42*w51*(w26*w28*w50 - w34*(3*fr_r*w1*w4 + frrr_r*w1*w20*w26 - w1*w14 - w1*w46 - w41 - w44*w45 - w47)),
      w51*(-w31*(3*fr_r*w2*w4 + frrr_r*w2*w20*w26 - w2*w40 - w2*w46 - 2*w38 - w45*w49 - w47) + w34*w50),
      0,
      0,
      0,
      0;
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
  LinIntp f_intp, fr_intp, frr_intp, frrr_intp, prf_intp;
};

#endif
