#ifndef __INTERPOL_HPP
#define __INTERPOL_HPP

#include <iostream>
#include <string>
#include "typedefs.hpp"
#include "PointValues.hpp"

//using namespace std;

// Time weight and rate between stamps; zero rate on a single stamp
inline double stamp_weight(const double t, const double t_prev, const double t_next){
  return (t_next > t_prev) ? (t - t_prev)/(t_next - t_prev) : 0.;
}
inline double stamp_rate(const double next, const double prev, const double t_prev, const double t_next){
  return (t_next > t_prev) ? (next - prev)/(t_next - t_prev) : 0.;
}
inline Vector3d stamp_rate(const Vector3d& next, const Vector3d& prev, const double t_prev, const double t_next){
  if (t_next > t_prev) return (next - prev)/(t_next - t_prev);
  return Vector3d::Zero();
}
inline Matrix3d stamp_rate(const Matrix3d& next, const Matrix3d& prev, const double t_prev, const double t_next){
  if (t_next > t_prev) return (next - prev)/(t_next - t_prev);
  return Matrix3d::Zero();
}

class Interpol {  // Abstract base class
public:
  Interpol(const std::string& infilename) { this->infilename=infilename; };
  virtual ~Interpol() = default;
  void set_folder(const std::string& folder){ this->folder=folder; };
  std::string get_folder() const { return folder; };
  void set_U0(const double U0) { this->U0 = U0; };
  double get_U0() { return this->U0; };
  void set_int_order(const int int_order) { this->int_order = int_order; check_gradient(); };
  // Gradient needed for int_order > 1 or carried elements
  void set_needs_gradient(const bool b) { needs_gradient_ = b; check_gradient(); };
  bool wants_gradient() const { return int_order > 1 || needs_gradient_; };
  // A field without a gradient refuses a run that wants one
  virtual void check_gradient() const {};
  //
  double get_Lx() { return x_max[0]-x_min[0]; };
  double get_Ly() { return x_max[1]-x_min[1]; };
  double get_Lz() { return x_max[2]-x_min[2]; };
  Vector3d get_x_min() const { return x_min; };
  Vector3d get_x_max() const { return x_max; };
  // Serial versions: no cell cache, time from update()
  bool locate(const Vector3d &x){ CellPos pos; return locate(x, t_update, pos); };
  bool evaluate(const Vector3d &x, PointValues& ptvals){
    CellPos pos;
    const bool inside = locate(x, t_update, pos);
    if (inside)
      evaluate(x, t_update, pos, ptvals);
    return inside;
  };
  // Locate only
  bool locate(const Vector3d &x, const double t, int& cell_id){
    CellPos pos;
    pos.id = cell_id;
    const bool inside = locate(x, t, pos);
    cell_id = pos.id;
    return inside;
  };
  //
  virtual double get_t_min() = 0;
  virtual double get_t_max() = 0;
  //
  virtual void update(const double t) = 0;
  // After locate, pos describes x in pos.id, inside or not; on failure id is unchanged
  virtual bool locate(const Vector3d &x, const double t, CellPos& pos) = 0;
  // Only after a successful locate: outside the fluid the velocity is zero
  virtual void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& ptvals) = 0;
  //
  virtual Vector3d get_boundary_normal(const Vector3d &x, int& cell_id) { return {0., 0., 0.}; }; // should be overloaded
  // Walk dx from x, located in pos, off the walls; dx and pos give the end
  virtual bool reflect(const Vector3d& x, Vector3d& dx, CellPos& pos) { return false; };
  // Build what reflect needs
  virtual void enable_reflection() {};
  bool can_reflect = false;
  virtual double hmin() const { return 0.; };   // 0: no mesh scale
protected:
  std::string infilename;
  std::string folder;
  bool is_initialized = false;
  bool verbose = true;
  int int_order = 1;
  bool needs_gradient_ = false;
  //double Lx = 0;
  //double Ly = 0;
  //double Lz = 0;
  Vector3d x_min;
  Vector3d x_max;
  double U0 = 1.0;
  double t_update;
};


#endif
