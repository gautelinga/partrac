#include "Expr.hpp"

#ifndef __EXPR_STOKESSPHERE_HPP
#define __EXPR_STOKESSPHERE_HPP

//using namespace std;

class Expr_StokesSphere final : public Expr {
public:
  Expr_StokesSphere(const partrac::Params& expr_params) : Expr(expr_params) {
    R = expr_params.get<double>("R");
    mu = expr_params.get<double>("mu");
    Rho = expr_params.get<double>("rho");
    x0 = {expr_params.get<double>("x0"),
          expr_params.get<double>("y0"),
          expr_params.get<double>("z0")};
    u_inf = expr_params.get<double>("u_inf");
    p_inf = expr_params.get<double>("p_inf");
  };
  // Keys this expression reads
  static void add_params(partrac::Schema& s) {
    s.require<double>("R", "sphere radius");
    s.require<double>("u_inf", "far-field velocity");
    s.require<double>("mu", "viscosity");
    s.require<double>("p_inf", "pressure");
    s.require<double>("rho", "density");
    s.require<double>("x0", "centre x");
    s.require<double>("y0", "centre y");
    s.require<double>("z0", "centre z");
  };
  // Distance to the sphere, positive in the fluid
  bool has_wall() const { return true; };
  double sdf(const Vector3d &x) const { return (x-x0).norm() - R; };
  Vector3d sdf_grad(const Vector3d &x) const { return (x-x0).normalized(); };
  bool inside(const Vector3d &x, const double t __attribute__((unused))) {
    Vector3d r = x-x0;
    double r2 = r.squaredNorm();
    double R2 = R*R;
    bool _is_inside = r2 >= R2;
    return _is_inside;
  };
  void eval(const Vector3d &x, const double t __attribute__((unused)), PointValues& ptvals) {

    Vector3d r = x-x0;
    double r2 = r.squaredNorm();
    double R2 = R*R;

    double eta = R/sqrt(r2);
    double eta3 = pow(eta, 3);

    double f_r = 1.0 + 1./2. * eta3 - 3./2. * eta;
    double f_theta = -1.0 + 1./4. * eta3 + 3./4. * eta;

    double r7 = pow(r2, 7.0/2.0);

    ptvals.U = {r[0]*r[2]/r2 * (f_r + f_theta) * u_inf,
                r[1]*r[2]/r2 * (f_r + f_theta) * u_inf,
                r[2]*r[2]/r2 * (f_r + f_theta) * u_inf - f_theta * u_inf};
    ptvals.P = p_inf - 3.0/2.0 * mu * u_inf * r[2] * eta / r2;
    ptvals.Rho = Rho;

    // Hardcoded -- copied from consistency-checked Sympy code
    ptvals.gradU << (3.0/4.0)*R*r[2]*u_inf*(-2*pow(r[0], 2)*(R2 - r2) - pow(r[0], 2)*(3*R2 - r2) + r2*(R2 - r2))/r7,
                    (3.0/4.0)*R*r[0]*r[1]*r[2]*u_inf*(-5*R2 + 3*r2)/r7,
                    (3.0/4.0)*R*r[0]*u_inf*(-2*pow(r[2], 2)*(R2 - r2) - pow(r[2], 2)*(3*R2 - r2) + r2*(R2 - r2))/r7,
                    (3.0/4.0)*R*r[0]*r[1]*r[2]*u_inf*(-5*R2 + 3*r2)/r7,
                    (3.0/4.0)*R*r[2]*u_inf*(-2*pow(r[1], 2)*(R2 - r2) - pow(r[1], 2)*(3*R2 - r2) + r2*(R2 - r2))/r7,
                    (3.0/4.0)*R*r[1]*u_inf*(-2*pow(r[2], 2)*(R2 - r2) - pow(r[2], 2)*(3*R2 - r2) + r2*(R2 - r2))/r7,
                    (3.0/4.0)*R*r[0]*u_inf*(-2*pow(r[2], 2)*(R2 - r2) - pow(r[2], 2)*(3*R2 - r2) + r2*(R2 + r2))/r7,
                    (3.0/4.0)*R*r[1]*u_inf*(-2*pow(r[2], 2)*(R2 - r2) - pow(r[2], 2)*(3*R2 - r2) + r2*(R2 + r2))/r7,
                    (3.0/4.0)*R*r[2]*u_inf*(-2*pow(r[2], 2)*(R2 - r2) - pow(r[2], 2)*(3*R2 - r2) + 2*r2*(R2 - r2) + r2*(R2 + r2))/r7;
  };
private:
  double R;  // Radius of sphere
  double mu;  // Viscosity
  Vector3d x0;  // Center of sphere
  double u_inf;  // Far-field velocity
  double p_inf;  // Far-field pressure
  double Rho;
};

#endif
