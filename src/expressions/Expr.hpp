#ifndef __EXPR_HPP
#define __EXPR_HPP

//using namespace std;

class Expr {
public:
  Expr(std::map<std::string, std::string> &expr_params) {
    this->expr_params = expr_params;
  };
  ~Expr() {};
  virtual void eval(const Vector3d &x, const double t) = 0;
  virtual bool inside() { return is_inside; };
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
  std::map<std::string, std::string> expr_params;
  bool is_inside = true;
};

#endif
