#ifndef __INTERPOL_HPP
#define __INTERPOL_HPP

#include "io.hpp"
#include "typedefs.hpp"

//using namespace std;

class Interpol {  // Abstract base class
public:
  Interpol(const std::string infilename) { this->infilename=infilename; };
  virtual ~Interpol() = default;
  void set_folder(std::string folder){ this->folder=folder; };
  std::string get_folder(){ return folder; };
  void set_U0(const double U0) { this->U0 = U0; this->U02 = U0*U0; };
  void set_int_order(const int int_order) { this->int_order = int_order; };
  //
  Vector3d get_u() { return { U0*get_ux(), U0 * get_uy(), U0 * get_uz()}; };
  Vector3d get_a() { return { U0 * get_ax(), U0 * get_ay(), U0 * get_az()}; };
  double get_Lx() { return x_max[0]-x_min[0]; };
  double get_Ly() { return x_max[1]-x_min[1]; };
  double get_Lz() { return x_max[2]-x_min[2]; };
  Vector3d get_x_min() const { return x_min; };
  Vector3d get_x_max() const { return x_max; };
  double get_divu() { return U0 * (get_uxx()+get_uyy()+get_uzz()); };
  double get_vortx() { return U0 * (get_uzy()-get_uyz()); };
  double get_vorty() { return U0 * (get_uxz()-get_uzx()); };
  double get_vortz() { return U0 * (get_uyx()-get_uxy()); };
  double get_Jux() { return U02 * (get_uxx()*get_ux()
                             + get_uxy()*get_uy()
                             + get_uxz()*get_uz()); };
  double get_Juy() { return U02 * (get_uyx()*get_ux()
                             + get_uyy()*get_uy()
                             + get_uyz()*get_uz()); };
  double get_Juz() { return U02 * (get_uzx()*get_ux()
                             + get_uzy()*get_uy()
                             + get_uzz()*get_uz()); };
  Vector3d get_Ju() { return {get_Jux(), get_Juy(), get_Juz()}; };
  Matrix3d get_J() {
    Matrix3d J;
    J <<
      U0 * get_uxx(), U0 * get_uxy(), U0 * get_uxz(),
      U0 * get_uyx(), U0 * get_uyy(), U0 * get_uyz(),
      U0 * get_uzx(), U0 * get_uzy(), U0 * get_uzz();
    return J; };
  void probe(const Vector3d &x){ probe(x, t_update); };
  //
  virtual double get_t_min() = 0;
  virtual double get_t_max() = 0;
  //
  virtual void update(const double t) = 0;
  virtual void probe(const Vector3d &x, const double t) = 0;
  //
  virtual bool inside_domain() = 0;
  virtual double get_ux() = 0;
  virtual double get_uy() = 0;
  virtual double get_uz() = 0;
  //
  virtual double get_ax() = 0;
  virtual double get_ay() = 0;
  virtual double get_az() = 0;
  //
  virtual double get_rho() = 0;
  virtual double get_p() = 0;
  virtual double get_uxx() = 0;
  virtual double get_uxy() = 0;
  virtual double get_uxz() = 0;
  virtual double get_uyx() = 0;
  virtual double get_uyy() = 0;
  virtual double get_uyz() = 0;
  virtual double get_uzx() = 0;
  virtual double get_uzy() = 0;
  virtual double get_uzz() = 0;
  //
  virtual Matrix3d get_grada() = 0;
  //
protected:
  std::string infilename;
  std::string folder;
  bool is_initialized = false;
  bool verbose = true;
  int int_order = 1;
  //double Lx = 0;
  //double Ly = 0;
  //double Lz = 0;
  Vector3d x_min;
  Vector3d x_max;
  double U0 = 1.0;
  double U02 = 1.0;
  double t_update;
};

#endif
