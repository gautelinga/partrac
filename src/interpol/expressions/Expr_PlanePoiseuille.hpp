#include "Expr.hpp"

#ifndef __EXPR_PLANEPOISEUILLE_HPP
#define __EXPR_PLANEPOISEUILLE_HPP

//using namespace std;

class Expr_PlanePoiseuille final : public Expr {
public:
  Expr_PlanePoiseuille(const partrac::Params& expr_params) : Expr(expr_params) {
    //R = expr_params.get<double>("R");
    x0 = {expr_params.get<double>("x0"),
          expr_params.get<double>("y0"),
          expr_params.get<double>("z0")};
    u_inf = expr_params.get<double>("u_inf");
    R = expr_params.get<double>("R");
    //alpha = expr_params.get<double>("alpha");
    mu = expr_params.get<double>("mu");
    p_inf = expr_params.get<double>("p_inf");
    Rho = expr_params.get<double>("rho");
  };
  // Keys this expression reads
  static void add_params(partrac::Schema& s) {
    s.require<double>("R", "channel half-width");
    s.require<double>("u_inf", "centreline velocity");
    s.require<double>("mu", "viscosity");
    s.require<double>("p_inf", "pressure");
    s.require<double>("rho", "density");
    s.require<double>("x0", "centre x");
    s.require<double>("y0", "centre y");
    s.require<double>("z0", "centre z");
  };
  // Distance to the nearer plate, positive in the fluid
  bool has_wall() const { return true; };
  double sdf(const Vector3d &x) const { return R - std::abs(x[0]-x0[0]); };
  Vector3d sdf_grad(const Vector3d &x) const { return {x[0] > x0[0] ? -1. : 1., 0., 0.}; };
  bool inside(const Vector3d &x, const double t __attribute__((unused))) {
    Vector3d r = x-x0;
    double chi = pow(r[0]/R, 2);
    bool _is_inside = chi <= 1.;
    return _is_inside;
  };
  void eval(const Vector3d &x, const double t __attribute__((unused)), PointValues& ptvals) {
    Vector3d r = x-x0;
    double chi = pow(r[0]/R, 2);
    ptvals.U = {0., 0., 3./2*u_inf*(1.0 - chi)}; // * alpha * r[0] * r[2];
    ptvals.P = p_inf;
    ptvals.Rho = Rho;
    ptvals.gradU << 0., 0., 0.,
                    0., 0., 0.,
                    - 3 * u_inf * r[0]/pow(R, 2),
                    0.,
                    0.;
  };
private:
  double mu;  // Viscosity
  //double alpha;  // alpha parameter = u''(x)
  Vector3d x0;  // origin
  double u_inf;  // Far-field velocity
  double R; //
  double p_inf;  // Far-field pressure
  double Rho;
};

#endif
