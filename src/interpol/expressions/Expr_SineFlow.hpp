
#include "Expr.hpp"
#include <stdexcept>
#include <vector>
#include "strings.hpp"

#ifndef __EXPR_SINEFLOW_HPP
#define __EXPR_SINEFLOW_HPP

// Comma separated list; empty if not all numbers
inline std::vector<int> int_list(const std::string& s){
  std::vector<int> v;
  for (const auto& tok : split_string(s, ",")){
    std::size_t end = 0;
    try { v.push_back(stoi(tok, &end)); } catch (const std::logic_error&) { return {}; }
    if (end != tok.size()) return {};
  }
  return v;
}
inline std::vector<double> real_list(const std::string& s){
  std::vector<double> v;
  for (const auto& tok : split_string(s, ",")){
    std::size_t end = 0;
    try { v.push_back(stod(tok, &end)); } catch (const std::logic_error&) { return {}; }
    if (end != tok.size()) return {};
  }
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
    s.check([](const partrac::Params& p){
              return !int_list(p.get<std::string>("flowdir")).empty()
                  && !int_list(p.get<std::string>("depdir")).empty()
                  && !real_list(p.get<std::string>("chi")).empty();
            },
            "flowdir, depdir and chi: lists of at least one entry, comma separated numbers");
    s.check([](const partrac::Params& p){
              for (const char* key : {"flowdir", "depdir"})
                for (const int d : int_list(p.get<std::string>(key)))
                  if (d < 0 || d > 2) return false;
              return true;
            },
            "flowdir and depdir: each entry 0, 1 or 2 (x, y or z)");
    // flowdir and depdir per phase, repeating
    s.check([](const partrac::Params& p){
              const std::vector<int> f = int_list(p.get<std::string>("flowdir"));
              const std::vector<int> d = int_list(p.get<std::string>("depdir"));
              for (std::size_t i = 0; i < f.size()*d.size(); ++i)
                if (f[i % f.size()] == d[i % d.size()]) return false;
              return true;
            },
            "sine_flow: a phase flows along the direction it varies along (flowdir equals depdir)");
    s.check([](const partrac::Params& p){ return p.get<double>("tau") > 0.; },
            "sine_flow: tau must be positive");
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
    const long i = long(std::floor(t/tau));   // phase; negative before t = 0
    double chi = chi_[phase(i, chi_.size())];
    int j = flowdir_[phase(i, flowdir_.size())];
    int k = depdir_[phase(i, depdir_.size())];
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
  // Index of phase i, repeating
  static std::size_t phase(const long i, const std::size_t n){
    const long m = i % long(n);
    return std::size_t(m < 0 ? m + long(n) : m);
  }
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
