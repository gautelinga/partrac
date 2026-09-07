#ifndef __INTERPOL_HPP
#define __INTERPOL_HPP

#include "io.hpp"
#include "typedefs.hpp"
#include "utils.hpp"

//using namespace std;

class Interpol {  // Abstract base class
public:
  Interpol(const std::string& infilename) { this->infilename=infilename; };
  //virtual ~Interpol() = default;
  virtual ~Interpol(){ std::cout << "Destructing Interpol." << std::endl; };
  void set_folder(const std::string& folder){ this->folder=folder; };
  std::string get_folder() const { return folder; };
  void set_U0(const double U0) { this->U0 = U0; };
  double get_U0() { return this->U0; };
  void set_int_order(const int int_order) { this->int_order = int_order; };
  //
  double get_Lx() { return x_max[0]-x_min[0]; };
  double get_Ly() { return x_max[1]-x_min[1]; };
  double get_Lz() { return x_max[2]-x_min[2]; };
  Vector3d get_x_min() const { return x_min; };
  Vector3d get_x_max() const { return x_max; };
  // For the serial initialization path: no cell cache, and the time last set
  // by update(). Both return whether the point is inside the domain.
  bool locate(const Vector3d &x){ int cell_id = -1; return locate(x, t_update, cell_id); };
  bool evaluate(const Vector3d &x, PointValues& ptvals){
    int cell_id = -1;
    const bool inside = locate(x, t_update, cell_id);
    evaluate(x, t_update, cell_id, ptvals);
    return inside;
  };
  //
  virtual double get_t_min() = 0;
  virtual double get_t_max() = 0;
  //
  virtual void update(const double t) = 0;
  virtual bool locate(const Vector3d &x, const double t, int& cell_id) = 0;
  virtual void evaluate(const Vector3d &x, const double t, const int cell_id, PointValues& ptvals) = 0;
  //
  virtual Vector3d get_boundary_normal(const Vector3d &x, int& cell_id) { return {0., 0., 0.}; }; // should be overloaded
  //
  template<typename T>
  void assign_fields(T&, const std::map<std::string, bool>& output_fields);
  virtual void reflect(Vector3d &x, Vector3d &dx_new, const double t, const double dt, int& cell_id) { };
  bool can_reflect = false;
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
  double t_update;
};

template<typename T>
void Interpol::assign_fields(T& ps, const std::map<std::string, bool>& output_fields){
  double t = t_update;
  #pragma omp parallel for
  for ( auto & particle : ps.particles() ){
    PointValues ptvals(get_U0());
    int cell_id = particle.cell_id();
    bool inside = true;
    if (cell_id == -1 )
      inside = locate(particle.get_x(), t, cell_id);
    if (inside) {
      evaluate(particle.get_x(), t, cell_id, ptvals);
      if (output_fields.find("u")->second)
        particle.u() = ptvals.get_u();
      if (output_fields.find("rho")->second)
        particle.rho() = ptvals.get_rho();
      if (output_fields.find("p")->second)
        particle.p() = ptvals.get_p();
      if (output_fields.find("cell_type")->second)
        particle.cell_type() = ptvals.get_cell_type();
      if (output_fields.find("J")->second)
        particle.J() = ptvals.get_J();
    }
  }
}

#endif
