#include "Interpol.hpp"
#include "Params.hpp"
#include "files.hpp"
#include "PointValues.hpp"
#include "expressions/Expr_StokesSphere.hpp"
#include "expressions/Expr_SineFlow.hpp"
#include "expressions/Expr_BatchelorVortex.hpp"
#include "expressions/Expr_ABCFlow.hpp"
#include "expressions/Expr_InfinitePlate.hpp"
#include "expressions/Expr_HagenPoiseuille.hpp"
#include "expressions/Expr_PlanePoiseuille.hpp"
#include "expressions/Expr_BrinkmanCylinder.hpp"
#include "expressions/Expr_TaylorCouette.hpp"
#include "expressions/Expr_LinearFlow.hpp"

#ifndef __ANALYTICINTERPOL_HPP
#define __ANALYTICINTERPOL_HPP

//using namespace std;

class AnalyticInterpol final : public Interpol {
public:
  AnalyticInterpol(const std::string infilename);
  void update(const double t) { this->t_update=t; };
  bool locate(const Vector3d &x, const double t, CellPos& pos) {
    return expr->inside(x, t);
  };
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& ptvals) {
    expr->eval(x, t, ptvals);
  };
  // Walk a move off the walls
  __attribute__((noinline))
  bool reflect(const Vector3d& x, Vector3d& dx, CellPos& pos) {
    constexpr int max_bounces = 8;
    Vector3d p = x;
    Vector3d d = dx;
    Vector3d walked = Vector3d::Zero();
    for (int bounce = 0; bounce <= max_bounces; ++bounce){
      if (expr->inside(x + (walked + d), 0.)){
        dx = walked + d;
        return true;
      }
      // Crossing by bisection, to 1e-9 of the step: sdf(p) >= 0 > sdf(p + d)
      double lo = 0., hi = 1.;
      for (int k = 0; k < 30; ++k){
        const double mid = 0.5*(lo + hi);
        (expr->sdf(p + mid*d) >= 0. ? lo : hi) = mid;
      }
      const Vector3d part = lo*d;
      Vector3d rest = d - part;
      p += part;
      walked += part;
      const Vector3d n = expr->sdf_grad(p).normalized();
      const double rn = n.dot(rest);
      if (rn < 0.)
        rest -= 2*rn*n;
      d = rest;
    }
    return false;
  };
  void enable_reflection() { can_reflect = expr->has_wall(); };
  double get_t_min() { return expr_params.get<double>("t_min"); };
  double get_t_max() { return expr_params.get<double>("t_max"); };
  using Interpol::locate;
  using Interpol::evaluate;
protected:
  partrac::Params expr_params;
  std::shared_ptr<Expr> expr;
};

// The expressions by name: their keys and their constructor
struct ExprKind {
  const char* name;
  const char* alias;
  void (*add_params)(partrac::Schema&);
  std::shared_ptr<Expr> (*make)(const partrac::Params&);
};

template<typename E>
std::shared_ptr<Expr> make_expr(const partrac::Params& prm){ return std::make_shared<E>(prm); }

inline const std::vector<ExprKind>& expr_kinds(){
  static const std::vector<ExprKind> kinds = {
    {"stokes_sphere", "StokesSphere", Expr_StokesSphere::add_params, make_expr<Expr_StokesSphere>},
    {"sine_flow", "SineFlow", Expr_SineFlow::add_params, make_expr<Expr_SineFlow>},
    {"batchelor_vortex", "BatchelorVortex", Expr_BatchelorVortex::add_params, make_expr<Expr_BatchelorVortex>},
    {"abc_flow", "ABCFlow", Expr_ABCFlow::add_params, make_expr<Expr_ABCFlow>},
    {"infinite_plate", "InfinitePlate", Expr_InfinitePlane::add_params, make_expr<Expr_InfinitePlane>},
    {"hagen_poiseuille", "HagenPoiseuille", Expr_HagenPoiseuille::add_params, make_expr<Expr_HagenPoiseuille>},
    {"plane_poiseuille", "PlanePoiseuille", Expr_PlanePoiseuille::add_params, make_expr<Expr_PlanePoiseuille>},
    {"brinkman_cylinder", "BrinkmanCylinder", Expr_BrinkmanCylinder::add_params, make_expr<Expr_BrinkmanCylinder>},
    {"taylor_couette", "TaylorCouette", Expr_TaylorCouette::add_params, make_expr<Expr_TaylorCouette>},
    {"linear_flow", "LinearFlow", Expr_LinearFlow::add_params, make_expr<Expr_LinearFlow>},
  };
  return kinds;
}

inline const ExprKind* find_expr_kind(const std::string& name){
  for (const auto& k : expr_kinds())
    if (name == k.name || name == k.alias) return &k;
  return nullptr;
}

// Keys of expr_params.dat for one expression
inline partrac::Schema expression_schema(const ExprKind& kind){
  partrac::Schema s(std::string("expr_params.dat, expression=") + kind.name, "");
  s.require<std::string>("expression", "the flow");
  s.require<double>("t_min", "start of the flow's time interval");
  s.require<double>("t_max", "end of the flow's time interval");
  s.require<double>("x_min", "domain box");
  s.require<double>("y_min", "domain box");
  s.require<double>("z_min", "domain box");
  s.require<double>("x_max", "domain box");
  s.require<double>("y_max", "domain box");
  s.require<double>("z_max", "domain box");
  kind.add_params(s);
  // Most example files carry the box lengths; only some flows read them
  for (const char* l : {"Lx", "Ly", "Lz"})
    if (!s.declares(l)) s.optional<double>(l, "box length, not read by this flow");
  return s;
}

inline AnalyticInterpol::AnalyticInterpol(const std::string infilename) : Interpol(infilename) {
  verify_file_exists(infilename);
  const std::string name = partrac::peek_file(infilename, "expression");
  if (name.empty()){
    std::cout << "No expression= in " << infilename << std::endl;
    exit(1);
  }
  const ExprKind* kind = find_expr_kind(name);
  if (!kind){
    std::cout << "Unknown expression " << name << " in " << infilename << std::endl;
    exit(1);
  }
  expr_params = partrac::parse_file_or_exit(expression_schema(*kind), infilename);

  std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  this->x_min << expr_params.get<double>("x_min"), expr_params.get<double>("y_min"), expr_params.get<double>("z_min");
  this->x_max << expr_params.get<double>("x_max"), expr_params.get<double>("y_max"), expr_params.get<double>("z_max");

  std::cout << kind->alias << " selected" << std::endl;
  expr = kind->make(expr_params);
}

#endif
