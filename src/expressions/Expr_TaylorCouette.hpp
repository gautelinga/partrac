#include "Expr.hpp"

#ifndef __EXPR_TAYLORCOUETTE_HPP
#define __EXPR_TAYLORCOUETTE_HPP

// Martinez-Ruiz et al., J. Fluid Mech. 837 (2018) 230-257, equation (3.2):
// the analytic fit to their measured Taylor-Couette flow, with the inner
// cylinder at r = 1 turning at unit speed and the outer cylinder, the top and
// the bottom at rest. Lengths are in units of the inner radius and time in
// units of 1/Omega, so the cell is 1 <= r <= R, |z| <= H/2 about x0.
//
// The meridional part comes from a stream function, u_r = -(1/r) dPsi/dz and
// u_z = (1/r) dPsi/dr, so it is divergence free by construction rather than
// to within a fit. Ekman pumping off the end plates is the a-mode; the two
// exponentials in u_theta carry the corner singularities where the turning
// inner cylinder slides past the stationary plates, c3 the bulk rotation.
class Expr_TaylorCouette final : public Expr {
public:
  Expr_TaylorCouette(std::map<std::string, std::string> &expr_params) : Expr(expr_params) {
    R = getd(expr_params, "R");
    H = getd(expr_params, "H");
    K = getd(expr_params, "K");
    a = getd(expr_params, "a");
    c1 = getd(expr_params, "c1");
    c2 = getd(expr_params, "c2");
    c3 = getd(expr_params, "c3");
    rho_inf = getd(expr_params, "rho");
    p_inf = getd(expr_params, "p_inf");
    x0 = {getd(expr_params, "x0"),
          getd(expr_params, "y0"),
          getd(expr_params, "z0")};
    kr = M_PI/(R-1);
    kz = 2*M_PI/H;
  };
  void eval(const Vector3d &x, const double t __attribute__((unused))) {
    is_inside = inside(x, 0.);
    compute(x, U, gradU);
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
  // Writes no members: the PointValues overload runs inside an omp for
  void compute(const Vector3d &x, Vector3d& U_, Matrix3d& gradU_) const {
    Vector3d d = x-x0;
    double s = sqrt(d[0]*d[0]+d[1]*d[1]);
    double z = d[2];

    // the walls are streamlines, so a particle only leaves them by integration
    // error; clamping keeps the corner terms finite for those excursions
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
    // gradU_(i, j) is dU_i/dx_j; the basis vectors turn with theta, which is
    // what the u_r/s and u_t/s terms are
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
