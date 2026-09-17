#include <boost/algorithm/string.hpp>
#include "Interpol.hpp"
#include "utils.hpp"
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
#include <fstream>

#ifndef __ANALYTICINTERPOL_HPP
#define __ANALYTICINTERPOL_HPP

//using namespace std;

class AnalyticInterpol : public Interpol {
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
  double get_t_min() { return getd(expr_params, "t_min"); };
  double get_t_max() { return getd(expr_params, "t_max"); };
  using Interpol::locate;
  using Interpol::evaluate;
protected:
  std::map<std::string, std::string> expr_params;
  std::shared_ptr<Expr> expr;
};

inline AnalyticInterpol::AnalyticInterpol(const std::string infilename) : Interpol(infilename) {
  std::ifstream input(infilename);
  if (!input){
    std::cout << "File " << infilename <<" doesn't exist." << std::endl;
    exit(1);
  }
  size_t found;
  std::string key, val;
  for (std::string line; getline(input, line); ){
    found = line.find('=');
    if (found != std::string::npos){
      key = line.substr(0, found);
      val = line.substr(found+1);
      boost::trim(key);
      boost::trim(val);
      expr_params[key] = val;
    }
  }

  std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  this->x_min << getd(expr_params, "x_min"), getd(expr_params, "y_min"), getd(expr_params, "z_min");
  this->x_max << getd(expr_params, "x_max"), getd(expr_params, "y_max"), getd(expr_params, "z_max");
  
  if (expr_params["expression"] == "stokes_sphere" ||
      expr_params["expression"] == "StokesSphere"){
    std::cout << "StokesSphere selected" << std::endl;
    expr = std::make_shared<Expr_StokesSphere>(expr_params);
  }
  else if (expr_params["expression"] == "sine_flow" ||
           expr_params["expression"] == "SineFlow"){
    std::cout << "SineFlow selected" << std::endl;
    expr = std::make_shared<Expr_SineFlow>(expr_params);
  }
  else if (expr_params["expression"] == "batchelor_vortex" ||
           expr_params["expression"] == "BatchelorVortex"){
    std::cout << "BatchelorVortex selected" << std::endl;
    expr = std::make_shared<Expr_BatchelorVortex>(expr_params);
  }
  else if (expr_params["expression"] == "abc_flow" ||
           expr_params["expression"] == "ABCFlow"){
    std::cout << "ABCFlow selected" << std::endl;
    expr = std::make_shared<Expr_ABCFlow>(expr_params);
  }
  else if (expr_params["expression"] == "infinite_plate" ||
           expr_params["expression"] == "InfinitePlate"){
    std::cout << "InfinitePlate selected" << std::endl;
    expr = std::make_shared<Expr_InfinitePlane>(expr_params);
  }
  else if (expr_params["expression"] == "hagen_poiseuille" ||
           expr_params["expression"] == "HagenPoiseuille"){
    std::cout << "HagenPoiseuille selected" << std::endl;
    expr = std::make_shared<Expr_HagenPoiseuille>(expr_params);
  }
  else if (expr_params["expression"] == "plane_poiseuille" ||
           expr_params["expression"] == "PlanePoiseuille"){
    std::cout << "PlanePoiseuille selected" << std::endl;
    expr = std::make_shared<Expr_PlanePoiseuille>(expr_params);
  }
  else if (expr_params["expression"] == "brinkman_cylinder" ||
           expr_params["expression"] == "BrinkmanCylinder"){
    std::cout << "BrinkmanCylinder selected" << std::endl;
    expr = std::make_shared<Expr_BrinkmanCylinder>(expr_params);
  }
  else if (expr_params["expression"] == "taylor_couette" ||
           expr_params["expression"] == "TaylorCouette"){
    std::cout << "TaylorCouette selected" << std::endl;
    expr = std::make_shared<Expr_TaylorCouette>(expr_params);
  }
  else if (expr_params["expression"] == "linear_flow" ||
           expr_params["expression"] == "LinearFlow"){
    std::cout << "LinearFlow selected" << std::endl;
    expr = std::make_shared<Expr_LinearFlow>(expr_params);
  }
  else {
    std::cout << "Could not find expression: " << expr_params["expression"] << std::endl;
    exit(1);
  }
}

#endif
