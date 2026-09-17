#ifndef __EXPR_HPP
#define __EXPR_HPP

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include "typedefs.hpp"
#include "Params.hpp"
#include "PointValues.hpp"
//using namespace std;

class Expr {
public:
  Expr(const partrac::Params& expr_params) {
    this->expr_params = expr_params;
  };
  ~Expr() {};
  virtual void eval(const Vector3d &x, const double t) = 0;
  virtual void eval(const Vector3d &x, const double t, PointValues& ptvals) = 0;
  virtual bool inside() { return is_inside; };
  virtual bool inside(const Vector3d &x, const double t) = 0;
  // Walls as a signed distance, positive in the fluid
  virtual bool has_wall() const { return false; };
  virtual double sdf(const Vector3d &x __attribute__((unused))) const { return 1.; };
  virtual Vector3d sdf_grad(const Vector3d &x __attribute__((unused))) const { return Vector3d::Zero(); };
  virtual double ux() { return 0.; };
  virtual double uy() { return 0.; };
  virtual double uz() { return 0.; };
  virtual double ax() { return 0.; };
  virtual double ay() { return 0.; };
  virtual double az() { return 0.; };
  virtual double rho() { return 1.; };
  virtual double p() { return 1.; };
  virtual double uxx() { return 0.; };
  virtual double uxy() { return 0.; };
  virtual double uxz() { return 0.; };
  virtual double uyx() { return 0.; };
  virtual double uyy() { return 0.; };
  virtual double uyz() { return 0.; };
  virtual double uzx() { return 0.; };
  virtual double uzy() { return 0.; };
  virtual double uzz() { return 0.; };
  virtual Matrix3d grada() {
    Matrix3d Da;
    Da <<
      0., 0., 0.,
      0., 0., 0.,
      0., 0., 0.;
    return Da; };
protected:
  partrac::Params expr_params;
  bool is_inside = true;
};

#endif
