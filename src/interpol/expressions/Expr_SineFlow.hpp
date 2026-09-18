
#include "Expr.hpp"
#include "strings.hpp"

#ifndef __EXPR_SINEFLOW_HPP
#define __EXPR_SINEFLOW_HPP

//using namespace std;

// Comma separated lists
inline std::vector<int> int_list(const std::string& s){
  std::vector<int> v;
  for (const auto& tok : split_string(s, ","))
    v.push_back(stoi(tok));
  return v;
}
inline std::vector<double> real_list(const std::string& s){
  std::vector<double> v;
  for (const auto& tok : split_string(s, ","))
    v.push_back(stod(tok));
  return v;
}

class Expr_SineFlow final : public Expr {
public:
  Expr_SineFlow(const partrac::Params& expr_params) : Expr(expr_params) {
    flowdir_ = int_list(expr_params.get<std::string>("flowdir"));
    depdir_ = int_list(expr_params.get<std::string>("depdir"));
    chi_ = real_list(expr_params.get<std::string>("chi"));
    u_inf = expr_params.get<double>("u_inf");
    p_inf = expr_params.get<double>("p_inf");
    tau = expr_params.get<double>("tau");
    rho_inf = expr_params.get<double>("rho");
    L = {expr_params.get<double>("Lx"),
         expr_params.get<double>("Ly"),
         expr_params.get<double>("Lz")};
  };
  // Keys this expression reads
  static void add_params(partrac::Schema& s) {
    s.require<std::string>("flowdir", "flow direction of each phase, comma separated");
    s.require<std::string>("depdir", "direction each phase varies along, comma separated");
    s.require<std::string>("chi", "phase shifts, comma separated");
    s.require<double>("u_inf", "velocity amplitude");
    s.require<double>("tau", "duration of each phase");
    s.require<double>("Lx", "domain length x");
    s.require<double>("Ly", "domain length y");
    s.require<double>("Lz", "domain length z");
    s.require<double>("p_inf", "pressure");
    s.require<double>("rho", "density");
  };
  bool inside(const Vector3d &x __attribute__((unused)), const double t __attribute__((unused))) {
    return true;
  };
  void eval(const Vector3d &x, const double t, PointValues& ptvals) {
    compute(x, t, ptvals.U, ptvals.gradU);
    ptvals.P = p_inf;
    ptvals.Rho = rho_inf;
  };
private:
  // No member writes: called inside omp for
  void compute(const Vector3d &x, const double t, Vector3d& U_, Matrix3d& gradU_) const {
    int i = floor(t/tau);
    double chi = chi_[i % chi_.size()];
    int j = flowdir_[i % flowdir_.size()];
    int k = depdir_[i % depdir_.size()];
    assert(j < 3 && j >= 0);
    assert(k < 3 && k >= 0);

    U_ = {0., 0., 0.};
    gradU_ <<
      0., 0., 0.,
      0., 0., 0.,
      0., 0., 0.;
    U_(j) = u_inf * sin(2*M_PI*x[k]/L[k] + chi);
    gradU_(j, k) = u_inf * 2*M_PI/L[k] * cos(2*M_PI*x[k]/L[k] + chi);
  };
  std::vector<double> chi_;
  std::vector<int> flowdir_;
  std::vector<int> depdir_;
  double tau;
  double u_inf;  // Far-field velocity
  double p_inf;  // Far-field pressure
  double rho_inf;
  Vector3d L;
};

#endif
