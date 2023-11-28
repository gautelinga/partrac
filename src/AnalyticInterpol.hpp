#include "Interpol.hpp"
#include "utils.hpp"
#include "expressions/Expr_StokesSphere.hpp"
#include "expressions/Expr_SineFlow.hpp"
#include "expressions/Expr_BatchelorVortex.hpp"
#include "expressions/Expr_ABCFlow.hpp"
#include "expressions/Expr_InfinitePlate.hpp"
#include "expressions/Expr_HagenPoiseuille.hpp"
#include "expressions/Expr_BrinkmanCylinder.hpp"
#include <fstream>

#ifndef __ANALYTICINTERPOL_HPP
#define __ANALYTICINTERPOL_HPP

//using namespace std;

class AnalyticInterpol : public Interpol {
public:
  AnalyticInterpol(const std::string infilename);
  void update(const double t) { this->t_update=t; };
  void probe(const Vector3d &x, const double t) { expr->eval(x, t); };
  void probe(const Vector3d &x, const double t, int& cell_id) { expr->eval(x, t); };
  bool probe_light(const Vector3d &x, const double t, int& cell_id) {
    return expr->inside(x, t);
  };
  void probe_heavy(const Vector3d &x, const double t, const int cell_id, PointValues& ptvals) {
    expr->eval(x, t, ptvals);
  };
  bool inside_domain() const { return expr->inside(); };
  double get_ux() { return expr->ux(); };
  double get_uy() { return expr->uy(); };
  double get_uz() { return expr->uz(); };
  double get_ax() { return expr->ax(); };
  double get_ay() { return expr->ay(); };
  double get_az() { return expr->az(); };
  double get_t_min() { return getd(expr_params, "t_min"); };
  double get_t_max() { return getd(expr_params, "t_max"); };
  double get_rho() { return expr->rho(); };
  double get_p() { return expr->p(); };
  double get_uxx() { return expr->uxx(); };
  double get_uxy() { return expr->uxy(); };
  double get_uxz() { return expr->uxz(); };
  double get_uyx() { return expr->uyx(); };
  double get_uyy() { return expr->uyy(); };
  double get_uyz() { return expr->uyz(); };
  double get_uzx() { return expr->uzx(); };
  double get_uzy() { return expr->uzy(); };
  double get_uzz() { return expr->uzz(); };
  Matrix3d get_grada() { return expr->grada(); };
  void probe(const Vector3d &x){ probe(x, this->t_update); };
protected:
  std::map<std::string, std::string> expr_params;
  std::shared_ptr<Expr> expr;
};

AnalyticInterpol::AnalyticInterpol(const std::string infilename) : Interpol(infilename) {
  std::ifstream input(infilename);
  if (!input){
    std::cout << "File " << infilename <<" doesn't exist." << std::endl;
    exit(0);
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
  else
  if (expr_params["expression"] == "brinkman_cylinder" ||
           expr_params["expression"] == "BrinkmanCylinder"){
    std::cout << "BrinkmanCylinder selected" << std::endl;
    expr = std::make_shared<Expr_BrinkmanCylinder>(expr_params);
  }
  else {
    std::cout << "Could not find expression: " << expr_params["expression"] << std::endl;
    exit(0);
  }
}

#endif
