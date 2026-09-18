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
  virtual void eval(const Vector3d &x, const double t, PointValues& ptvals) = 0;
  virtual bool inside(const Vector3d &x, const double t) = 0;
  // Walls as a signed distance, positive in the fluid
  virtual bool has_wall() const { return false; };
  virtual double sdf(const Vector3d &x __attribute__((unused))) const { return 1.; };
  virtual Vector3d sdf_grad(const Vector3d &x __attribute__((unused))) const { return Vector3d::Zero(); };
protected:
  partrac::Params expr_params;
};

#endif
