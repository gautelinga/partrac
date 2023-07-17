
#include "Expr.hpp"

#ifndef __EXPR_SINEFLOW_HPP
#define __EXPR_SINEFLOW_HPP

//using namespace std;

class Expr_SineFlow : public Expr {
public:
  Expr_SineFlow(std::map<std::string, std::string> &expr_params) : Expr(expr_params) {
    flowdir_ = getivec(expr_params, "flowdir");
    depdir_ = getivec(expr_params, "depdir");
    chi_ = getdvec(expr_params, "chi");
    u_inf = getd(expr_params, "u_inf");
    p_inf = getd(expr_params, "p_inf");
    tau = getd(expr_params, "tau");
    rho_inf = getd(expr_params, "rho");
    L = {getd(expr_params, "Lx"),
         getd(expr_params, "Ly"),
         getd(expr_params, "Lz")};
  };
  void eval(const Vector3d &x, const double t) {
    is_inside = true;

    int i = floor(t/tau);
    double chi = chi_[i % chi_.size()];
    int j = flowdir_[i % flowdir_.size()];
    int k = depdir_[i % flowdir_.size()];
    assert(j < 2 && j >= 0);
    assert(k < 2 && k >= 0);
    // cout << ":: " << chi << " " << j << endl;

    U = {0., 0., 0.};
    Ujk <<
      0., 0., 0.,
      0., 0., 0.,
      0., 0., 0.;
    U(j) = u_inf * sin(2*M_PI*x[k]/L[k] + chi);
    Ujk(j, k) = u_inf * 2*M_PI/L[k] * cos(2*M_PI*x[k]/L[k] + chi);
  };
  bool inside(const Vector3d &x, const double t __attribute__((unused))) {
    std::cout << "Not implemented yet!" << std::endl;
    exit(0);
    return false;
  };
  void eval(const Vector3d &x, const double t __attribute__((unused)), PointValues& ptvals) {
    std::cout << "Not implemented yet!" << std::endl;
    exit(0);
  };
  double ux() { return U(0); };
  double uy() { return U(1); };
  double uz() { return U(2); };
  double rho() { return rho_inf; };
  double p() { return p_inf; };
  double uxx() { return Ujk(0, 0); };
  double uxy() { return Ujk(0, 1); };
  double uxz() { return Ujk(0, 2); };
  double uyx() { return Ujk(1, 0); };
  double uyy() { return Ujk(1, 1); };
  double uyz() { return Ujk(1, 2); };
  double uzx() { return Ujk(2, 0); };
  double uzy() { return Ujk(2, 1); };
  double uzz() { return Ujk(2, 2); };
private:
  std::vector<double> chi_;
  std::vector<int> flowdir_;
  std::vector<int> depdir_;
  double tau;
  double u_inf;  // Far-field velocity
  double p_inf;  // Far-field pressure
  double rho_inf;
  // Useful quantitites
  Vector3d U;
  Matrix3d Ujk;
  Vector3d L;
};

#endif
