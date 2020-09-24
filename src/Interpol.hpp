#include "io.hpp"

#ifndef __INTERPOL_HPP
#define __INTERPOL_HPP

//using namespace std;

class Interpol {  // Abstract base class
public:
  Interpol(const std::string infilename) { this->infilename=infilename; };
  ~Interpol() {};
  void set_folder(std::string folder){ this->folder=folder; };
  std::string get_folder(){ return folder; };
  Vector3d get_u() { return {get_ux(), get_uy(), get_uz()}; };
  Vector3d get_a() { return {get_ax(), get_ay(), get_az()}; };
  double get_Lx() { return Lx; };
  double get_Ly() { return Ly; };
  double get_Lz() { return Lz; };
  Uint get_nx() { return nx; };
  Uint get_ny() { return ny; };
  Uint get_nz() { return nz; };
  double get_divu() { return get_uxx()+get_uyy()+get_uzz(); };
  double get_vortx() { return get_uzy()-get_uyz(); };
  double get_vorty() { return get_uxz()-get_uzx(); };
  double get_vortz() { return get_uyx()-get_uxy(); };
  double get_Jux() { return (get_uxx()*get_ux()
                             + get_uxy()*get_uy()
                             + get_uxz()*get_uz()); };
  double get_Juy() { return (get_uyx()*get_ux()
                             + get_uyy()*get_uy()
                             + get_uyz()*get_uz()); };
  double get_Juz() { return (get_uzx()*get_ux()
                             + get_uzy()*get_uy()
                             + get_uzz()*get_uz()); };
  Vector3d get_Ju() { return {get_Jux(), get_Juy(), get_Juz()}; };
  Matrix3d get_J() {
    Matrix3d J;
    J <<
      get_uxx(), get_uxy(), get_uxz(),
      get_uyx(), get_uyy(), get_uyz(),
      get_uzx(), get_uzy(), get_uzz();
    return J; };
  //
  virtual double get_t_min() = 0;
  virtual double get_t_max() = 0;
  //
  virtual void update(const double t) = 0;
  virtual void probe(const Vector3d &x) = 0;
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
protected:
  std::string infilename;
  std::string folder;
  bool is_initialized = false;
  bool verbose = true;
  double Lx = 0;
  double Ly = 0;
  double Lz = 0;
  Uint nx = 0;
  Uint ny = 0;
  Uint nz = 0;
};

#endif
