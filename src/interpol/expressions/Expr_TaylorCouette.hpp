#include "Expr.hpp"

#ifndef __EXPR_TAYLORCOUETTE_HPP
#define __EXPR_TAYLORCOUETTE_HPP

// Martinez-Ruiz et al., J. Fluid Mech. 837 (2018), eq. (3.2); lengths in inner radii, time in 1/Omega
class Expr_TaylorCouette final : public Expr {
public:
  Expr_TaylorCouette(const partrac::Params& expr_params) : Expr(expr_params) {
    R = expr_params.get<double>("R");
    H = expr_params.get<double>("H");
    K = expr_params.get<double>("K");
    a = expr_params.get<double>("a");
    c1 = expr_params.get<double>("c1");
    c2 = expr_params.get<double>("c2");
    c3 = expr_params.get<double>("c3");
    rho_inf = expr_params.get<double>("rho");
    p_inf = expr_params.get<double>("p_inf");
    x0 = {expr_params.get<double>("x0"),
          expr_params.get<double>("y0"),
          expr_params.get<double>("z0")};
    kr = M_PI/(R-1);
    kz = 2*M_PI/H;
  };
  // Keys this expression reads
  static void add_params(partrac::Schema& s) {
    s.require<double>("R", "outer radius, in inner radii");
    s.require<double>("H", "height");
    s.require<double>("K", "roll amplitude");
    s.require<double>("a", "second harmonic of the rolls");
    s.require<double>("c1", "decay of the azimuthal velocity from the inner cylinder");
    s.require<double>("c2", "end-plate decay");
    s.require<double>("c3", "azimuthal modulation by the rolls");
    s.require<double>("p_inf", "pressure");
    s.require<double>("rho", "density");
    s.require<double>("x0", "centre x");
    s.require<double>("y0", "centre y");
    s.require<double>("z0", "centre z");
  };
  void eval(const Vector3d &x, const double t __attribute__((unused))) {
    is_inside = inside(x, 0.);
    compute(x, U, gradU);
  };
  // Distance to the nearest of cylinders and plates, positive in the fluid
  bool has_wall() const { return true; };
  double sdf(const Vector3d &x) const {
    const Vector3d r = x-x0;
    const double s = std::hypot(r[0], r[1]);
    return std::min({s - 1., R - s, H/2 - std::abs(r[2])});
  };
  Vector3d sdf_grad(const Vector3d &x) const {
    const Vector3d r = x-x0;
    const double s = std::hypot(r[0], r[1]);
    const double inner = s - 1., outer = R - s, plate = H/2 - std::abs(r[2]);
    if (plate <= inner && plate <= outer) return {0., 0., r[2] > 0. ? -1. : 1.};
    const Vector3d e = s > 0. ? Vector3d(r[0]/s, r[1]/s, 0.) : Vector3d::Zero();
    return inner <= outer ? e : Vector3d(-e);
  };
  bool inside(const Vector3d &x, const double t __attribute__((unused))) {
    Vector3d r = x-x0;
    double s = sqrt(r[0]*r[0]+r[1]*r[1]);
    return s >= 1. && s <= R && 2*std::abs(r[2]) <= H;
  };
  void eval(const Vector3d &x, const double t __attribute__((unused)), PointValues& ptvals) {
    compute(x, ptvals.U, ptvals.gradU);
    ptvals.P = p_inf;
    ptvals.Rho = rho_inf;
  };
  double ux() { return U(0); };
  double uy() { return U(1); };
  double uz() { return U(2); };
  double rho() { return rho_inf; };
  double p() { return p_inf; };
  double uxx() { return gradU(0, 0); };
  double uxy() { return gradU(0, 1); };
  double uxz() { return gradU(0, 2); };
  double uyx() { return gradU(1, 0); };
  double uyy() { return gradU(1, 1); };
  double uyz() { return gradU(1, 2); };
  double uzx() { return gradU(2, 0); };
  double uzy() { return gradU(2, 1); };
  double uzz() { return gradU(2, 2); };
private:
  // No member writes: called inside omp for
  void compute(const Vector3d &x, Vector3d& U_, Matrix3d& gradU_) const {
    Vector3d d = x-x0;
    double s = sqrt(d[0]*d[0]+d[1]*d[1]);
    double z = d[2];

    // Clamp to the walls
    double dr = std::max(s-1., 0.);
    double t1 = tanh(c2*std::max(H/2-z, 0.));
    double t2 = tanh(c2*std::max(H/2+z, 0.));
    double co = 1./t1 + 1./t2;                       // infinite at either plate
    double eE = dr > 0. ? exp(-c1*dr*co) : 1.;

    double q = kr*dr, p = kz*z;
    double cq = cos(q), sq = sin(q);
    double F = cos(p) + a*cos(2*p);                  // dG/dz / kz
    double G = sin(p) + a*sin(2*p)/2;
    double Fz = -kz*(sin(p) + 2*a*sin(2*p));
    double W = 1. - 4*z*z/(H*H);
    double A = K*kz, B = K*kr;

    double u_r = A*sq*F/s;
    double u_z = -B*cq*G/s;
    double u_t = eE + c3*sq*W;

    double dur_dr = A*F*(kr*cq/s - sq/(s*s));
    double dur_dz = A*sq*Fz/s;
    double duz_dr = B*G*(kr*sq/s + cq/(s*s));
    double duz_dz = -B*kz*cq*F/s;
    double dut_dr = c3*kr*cq*W;
    double dut_dz = -8*c3*sq*z/(H*H);
    if (eE > 0. && std::isfinite(co)){
      dut_dr += -c1*co*eE;
      dut_dz += -c1*c2*dr*((1-t1*t1)/(t1*t1) - (1-t2*t2)/(t2*t2))*eE;
    }

    double c = d[0]/s, sn = d[1]/s;
    U_ = {u_r*c - u_t*sn, u_r*sn + u_t*c, u_z};
    // gradU_(i, j) is dU_i/dx_j
    gradU_ <<
      c*c*dur_dr - c*sn*dut_dr + (sn*sn*u_r + sn*c*u_t)/s,
      sn*c*dur_dr - sn*sn*dut_dr - (sn*c*u_r + c*c*u_t)/s,
      c*dur_dz - sn*dut_dz,
      c*sn*dur_dr + c*c*dut_dr - (sn*c*u_r - sn*sn*u_t)/s,
      sn*sn*dur_dr + sn*c*dut_dr + (c*c*u_r - sn*c*u_t)/s,
      sn*dur_dz + c*dut_dz,
      c*duz_dr,
      sn*duz_dr,
      duz_dz;
  };
  double R, H, K, a, c1, c2, c3;
  double kr, kz;
  double p_inf;
  double rho_inf;
  Vector3d x0;
  // Useful quantitites
  Vector3d U;
  Matrix3d gradU;
};

#endif
