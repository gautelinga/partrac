#include "Expr.hpp"

#ifndef __EXPR_BATCHELORVORTEX_HPP
#define __EXPR_BATCHELORVORTEX_HPP

//using namespace std;

class Expr_BatchelorVortex final : public Expr {
public:
  Expr_BatchelorVortex(const partrac::Params& expr_params) : Expr(expr_params) {
    R1 = expr_params.get<double>("R1");
    R2 = expr_params.get<double>("R2");
    R12 = R1*R1;
    R22 = R2*R2;
    q = expr_params.get<double>("q");
    rho_inf = expr_params.get<double>("rho");
    x0 = {expr_params.get<double>("x0"),
          expr_params.get<double>("y0"),
          expr_params.get<double>("z0")};
    u0 = expr_params.get<double>("u0");
    p_inf = expr_params.get<double>("p_inf");
  };
  // Keys this expression reads
  static void add_params(partrac::Schema& s) {
    s.require<double>("u0", "velocity scale");
    s.require<double>("q", "swirl parameter");
    s.require<double>("R1", "radius of the first gaussian");
    s.require<double>("R2", "radius of the second gaussian");
    s.require<double>("p_inf", "pressure");
    s.require<double>("rho", "density");
    s.require<double>("x0", "centre x");
    s.require<double>("y0", "centre y");
    s.require<double>("z0", "centre z");
  };
  bool inside(const Vector3d &x __attribute__((unused)), const double t __attribute__((unused))) {
    return true;
  };
  void eval(const Vector3d &x, const double t __attribute__((unused)), PointValues& ptvals) {
    compute(x, ptvals.U, ptvals.gradU, ptvals.P);
    ptvals.Rho = rho_inf;
  };
private:
  // phi(a) = (1 - exp(-a))/a and its derivative, both regular at a = 0
  static double phi(const double a) {
    if (a < 1e-8) return 1. - a/2. + a*a/6.;
    return -std::expm1(-a)/a;
  };
  static double dphi(const double a) {
    if (a < 1e-3) return -0.5 + a/3. - a*a/8.;
    return (a*std::exp(-a) + std::expm1(-a))/(a*a);
  };
  // No member writes: called inside omp for. Regular on the axis
  void compute(const Vector3d &x, Vector3d& U_, Matrix3d& gradU_, double& P_) const {
    Vector3d r = x-x0;
    double s2 = r[0]*r[0]+r[1]*r[1];
    double a = s2/R12;
    double ph = phi(a);
    double dph = dphi(a);
    double eta22 = exp(-s2/R22);

    U_ = {-u0*r[1]*ph/R1,
           u0*r[0]*ph/R1,
           q*u0*eta22};
    P_ = p_inf;

    double c = 2*u0*dph/(R1*R12);
    // gradU_(i, j) is dU_i/dx_j
    gradU_ <<
      -c*r[0]*r[1],                        -u0*(ph + 2*r[1]*r[1]*dph/R12)/R1, 0.,
       u0*(ph + 2*r[0]*r[0]*dph/R12)/R1,    c*r[0]*r[1],                      0.,
      -2*r[0]*q*u0*eta22/R22,              -2*r[1]*q*u0*eta22/R22,            0.;
  };
  double R1;  // Radius of sphere
  double R2;
  double R12;
  double R22;
  Vector3d x0;  // Center of vortex
  double u0;  // Far-field velocity
  double q; // Amplification
  double p_inf;  // Far-field pressure
  double rho_inf;
};

#endif
