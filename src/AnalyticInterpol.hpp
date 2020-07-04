#include "Interpol.hpp"
#include "utils.hpp"
#include "expressions/Expr_StokesSphere.hpp"
#include <fstream>

#ifndef __ANALYTICINTERPOL_HPP
#define __ANALYTICINTERPOL_HPP

using namespace std;

class AnalyticInterpol : public Interpol {
public:
  AnalyticInterpol(const string infilename);
  void update(const double t) { this->t=t; };
  void probe(const Vector3d &x) { expr->eval(x, t); };
  bool inside_domain() { return expr->inside(); };
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
protected:
  double t;
  map<string, string> expr_params;
  Expr* expr;
};

AnalyticInterpol::AnalyticInterpol(const string infilename) : Interpol(infilename) {
  ifstream input(infilename);
  if (!input){
    cout << "File " << infilename <<" doesn't exist." << endl;
    exit(0);
  }
  size_t found;
  string key, val;
  for (string line; getline(input, line); ){
    found = line.find('=');
    if (found != string::npos){
      key = line.substr(0, found);
      val = line.substr(found+1);
      boost::trim(key);
      boost::trim(val);
      expr_params[key] = val;
    }
  }

  std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  this->t = get_t_min();
  this->nx = 0;
  this->ny = 0;
  this->nz = 0;
  Lx = getd(expr_params, "Lx");
  Ly = getd(expr_params, "Ly");
  Lz = getd(expr_params, "Lz");

  if (expr_params["expression"] == "stokes_sphere" || expr_params["expression"] == "StokesSphere"){
    cout << "StokesSphere selected" << endl;
    expr = new Expr_StokesSphere(expr_params);
  }
  else {
    cout << "Could not find expression: " << expr_params[""] << endl;
    exit(0);
  }
}

#endif
