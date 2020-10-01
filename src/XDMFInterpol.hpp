#include <memory>
#include "Interpol.hpp"
#include "Timestamps.hpp"

#include "H5Cpp.h"
#define hid_t aa_hid_t
#undef hid_t
#include <dolfin.h>
#define hid_t bb_hid_t
//#include <dolfin/io/HDF5File.h>
//#include <dolfin/io/XDMFFile.h>
#undef hid_t
#define hid_t ambiguous use aa_hid_t or bb_hid_t

#include "dolfin_elements/vP1_2.h"
#include "dolfin_elements/vP2_2.h"
#include "dolfin_elements/vP3_2.h"
#include "dolfin_elements/vP1_3.h"
#include "dolfin_elements/vP2_3.h"
#include "dolfin_elements/vP3_3.h"
#include "dolfin_elements/P1_2.h"
#include "dolfin_elements/P2_2.h"
#include "dolfin_elements/P3_2.h"
#include "dolfin_elements/P1_3.h"
#include "dolfin_elements/P2_3.h"
#include "dolfin_elements/P3_3.h"

#ifndef __XDMFINTERPOL_HPP
#define __XDMFINTERPOL_HPP

using namespace H5;
//using namespace std;
//using namespace dolfin;


class PBC : public dolfin::SubDomain {
public:
  PBC(const std::vector<bool> &periodic,
      const Vector3d &x_min, const Vector3d &x_max, const Uint dim) {
    periodic_x = periodic[0] && dim > 0;
    periodic_y = periodic[1] && dim > 1;
    periodic_z = periodic[2] && dim > 2;
    this->x_min = x_min;
    this->x_max = x_max;
  };
  bool inside(const dolfin::Array<double>& x, bool on_boundary) const {
    return (on_boundary &&
            ((periodic_x && x[0] < x_min[0] + DOLFIN_EPS_LARGE) ||
             (periodic_y && x[1] < x_min[1] + DOLFIN_EPS_LARGE) ||
             (periodic_z && x[2] < x_min[2] + DOLFIN_EPS_LARGE)) &&
            !(
              (periodic_x && periodic_y &&
               x[0] < x_min[0] + DOLFIN_EPS_LARGE &&
               x[1] > x_max[1] - DOLFIN_EPS_LARGE) ||
              (periodic_x && periodic_y &&
               x[0] > x_max[0] - DOLFIN_EPS_LARGE &&
               x[1] < x_min[1] + DOLFIN_EPS_LARGE) ||
              (periodic_x && periodic_z &&
               x[0] < x_min[0] + DOLFIN_EPS_LARGE &&
               x[2] > x_max[2] - DOLFIN_EPS_LARGE) ||
              (periodic_x && periodic_z &&
               x[0] > x_max[0] - DOLFIN_EPS_LARGE &&
               x[2] < x_min[2] + DOLFIN_EPS_LARGE) ||
              (periodic_y && periodic_z &&
               x[1] < x_min[1] + DOLFIN_EPS_LARGE &&
               x[2] > x_max[2] - DOLFIN_EPS_LARGE) ||
              (periodic_y && periodic_z &&
               x[1] > x_max[1] - DOLFIN_EPS_LARGE &&
               x[2] < x_min[2] + DOLFIN_EPS_LARGE)
              )
            );
  }
  void map(const dolfin::Array<double>& x, dolfin::Array<double>& y) const {
    if (periodic_x && periodic_y && periodic_z &&
        x[0] > x_max[0] - DOLFIN_EPS_LARGE &&
        x[1] > x_max[1] - DOLFIN_EPS_LARGE &&
        x[2] > x_max[2] - DOLFIN_EPS_LARGE){
      y[0] = x[0] - (x_max[0]-x_min[0]);
      y[1] = x[1] - (x_max[1]-x_min[1]);
      y[2] = x[2] - (x_max[2]-x_min[2]);
    }
    else if (periodic_x && periodic_y &&
             x[0] > x_max[0] - DOLFIN_EPS_LARGE &&
             x[1] > x_max[1] - DOLFIN_EPS_LARGE){
      y[0] = x[0] - (x_max[0]-x_min[0]);
      y[1] = x[1] - (x_max[1]-x_min[1]);
      y[2] = x[2];
    }
    else if (periodic_x && periodic_z &&
             x[0] > x_max[0] - DOLFIN_EPS_LARGE &&
             x[2] > x_max[2] - DOLFIN_EPS_LARGE){
      y[0] = x[0] - (x_max[0]-x_min[0]);
      y[1] = x[1];
      y[2] = x[2] - (x_max[2]-x_min[2]);
    }
    else if (periodic_y && periodic_z && 
             x[1] > x_max[1] - DOLFIN_EPS_LARGE &&
             x[2] > x_max[2] - DOLFIN_EPS_LARGE){
      y[0] = x[0];
      y[1] = x[1] - (x_max[1]-x_min[1]);
      y[2] = x[2] - (x_max[2]-x_min[2]);
    }
    else if (periodic_x && x[0] > x_max[0] - DOLFIN_EPS_LARGE){
      y[0] = x[0] - (x_max[0]-x_min[0]);
      y[1] = x[1];
      y[2] = x[2];
    }
    else if (periodic_y && x[1] > x_max[1] - DOLFIN_EPS_LARGE){
      y[0] = x[0];
      y[1] = x[1] - (x_max[1]-x_min[1]);
      y[2] = x[2];
    }
    else if (periodic_z && x[2] > x_max[2] - DOLFIN_EPS_LARGE){
      y[0] = x[0];
      y[1] = x[1];
      y[2] = x[2] - (x_max[2]-x_min[2]);
    }
    else {
      y[0] = x[0]-1000;
      y[1] = x[1]-1000;
      y[2] = x[2]-1000;
    }
  }
private:
  double periodic_x;
  double periodic_y;
  double periodic_z;
  Vector3d x_min;
  Vector3d x_max;
};



class XDMFInterpol : public Interpol {
public:
  XDMFInterpol(const std::string infilename);
  void update(const double t);
  void probe(const Vector3d &x);
  bool inside_domain() { return inside; };
  double get_ux(){ return U[0]; };
  double get_uy(){ return U[1]; };
  double get_uz(){ return U[2]; };
  double get_ax() { return A[0]; };
  double get_ay() { return A[1];};
  double get_az() { return A[2];};
  double get_t_min() { return ts.get_t_min(); };
  double get_t_max() { return ts.get_t_max(); };
  double get_rho() {
    if (contains(dolfin_params, std::string("rho")))
      return stod(dolfin_params["rho"]);
    else {
      std::cout << "dolfin_params does not contain \"rho\"" << std::endl;
      exit(0);
    }
  };
  double get_p() { return P;};
  double get_uxx() { return gradU(0, 0); };
  double get_uxy() { return gradU(0, 1); };
  double get_uxz() { return gradU(0, 2); };
  double get_uyx() { return gradU(1, 0); };
  double get_uyy() { return gradU(1, 1); };
  double get_uyz() { return gradU(1, 2); };
  double get_uzx() { return gradU(2, 0); };
  double get_uzy() { return gradU(2, 1); };
  double get_uzz() { return gradU(2, 2); };
  Matrix3d get_grada() {
    return gradA;
  };
protected:
  Timestamps ts;
  double t_prev = 0.;
  double t_next = 0.;
  double alpha_t;

  std::vector<bool> periodic = {false, false, false};
  //Vector3d x_min = {0., 0., 0.};
  //Vector3d x_max = {0., 0., 0.};

  Vector3d U = {0., 0., 0.};
  //double Uy = 0.;
  //double Uz = 0.;
  Vector3d A = {0., 0., 0.};
  double P = 0.;
  Matrix3d gradU, gradA;

  bool inside;

  std::map<std::string, std::string> dolfin_params;

  std::shared_ptr<dolfin::Mesh> mesh;
  Uint dim;

  std::shared_ptr<dolfin::FunctionSpace> u_space;
  std::shared_ptr<dolfin::FunctionSpace> p_space;

  std::shared_ptr<dolfin::Function> u_prev_;
  std::shared_ptr<dolfin::Function> u_next_;
  std::shared_ptr<dolfin::Function> p_prev_;
  std::shared_ptr<dolfin::Function> p_next_;
};

XDMFInterpol::XDMFInterpol(const std::string infilename) : Interpol(infilename) {
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
      dolfin_params[key] = val;
    }
  }

  std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  ts.initialize(get_folder() + "/" + dolfin_params["timestamps"]);

  if (dolfin_params["periodic_x"] == "true"){
    periodic[0] = true;
  }
  if (dolfin_params["periodic_y"] == "true"){
    periodic[1] = true;
  }
  if (dolfin_params["periodic_z"] == "true"){
    periodic[2] = true;
  }

  //std::cout << dolfin_params["velocity"] << std::endl;

  //std::string xdmffname = get_folder() + "/" + dolfin_params["velocity"];
  //dolfin::XDMFFile xdmff(MPI_COMM_WORLD, xdmffname);
  //std::string h5filename = get_folder() + "/" + ts.get(0).prev.filename;
  //std::cout << h5filename << std::endl;
  //MPI_Comm mpi_comm;

  std::string meshfilename = get_folder() + "/" + dolfin_params["mesh"];
  dolfin::HDF5File meshfile(MPI_COMM_WORLD, meshfilename, "r");

  dolfin::Mesh mesh_in;
  meshfile.read(mesh_in, "mesh", false);

  mesh = std::make_shared<dolfin::Mesh>(mesh_in);
  dim = mesh->geometry().dim();

  std::vector<double> xx = mesh->coordinates();

  for (Uint i=0; i<dim; ++i){
    x_min[i] = xx[i];
    x_max[i] = xx[i];
  }

  for (Uint i=0; i<xx.size(); ++i){
    Uint i_loc = i % dim;
    x_min[i_loc] = std::min(x_min[i_loc], xx[i]);
    x_max[i_loc] = std::max(x_max[i_loc], xx[i]);
  }
  //std::cout << x_min << std::endl;
  //std::cout << x_max << std::endl;
  //this->Lx = x_max[0]-x_min[0];
  //this->Ly = x_max[1]-x_min[1];
  //this->Lz = x_max[2]-x_min[2];

  auto constrained_domain = std::make_shared<PBC>(periodic, x_min, x_max, dim);

  std::string u_el = dolfin_params["velocity_space"];
  std::string p_el = dolfin_params["pressure_space"];

  switch (dim){
  case 2:
    // Velocity
    if (u_el == "P1"){
      u_space = std::make_shared<vP1_2::FunctionSpace>(mesh, constrained_domain);
    }
    else if (u_el == "P2"){
      u_space = std::make_shared<vP2_2::FunctionSpace>(mesh, constrained_domain);
    }
    else if (u_el == "P3"){
      u_space = std::make_shared<vP3_2::FunctionSpace>(mesh, constrained_domain);
    }
    else {
      std::cout << "Unrecognized velocity element: " << u_el << std::endl;
      exit(0);
    }
    // Pressure
    if (p_el == "P1"){
      p_space = std::make_shared<P1_2::FunctionSpace>(mesh, constrained_domain);
    }
    else if (u_el == "P2"){
      p_space = std::make_shared<P2_2::FunctionSpace>(mesh, constrained_domain);
    }
    else if (u_el == "P3"){
      p_space = std::make_shared<P3_2::FunctionSpace>(mesh, constrained_domain);
    }
    else {
      std::cout << "Unrecognized pressure element: " << u_el << std::endl;
      exit(0);
    }
    break;
  case 3:
    // Velocity
    if (u_el == "P1"){
      u_space = std::make_shared<vP1_3::FunctionSpace>(mesh, constrained_domain);
    }
    else if (u_el == "P2"){
      u_space = std::make_shared<vP2_3::FunctionSpace>(mesh, constrained_domain);
    }
    else if (u_el == "P3"){
      u_space = std::make_shared<vP3_3::FunctionSpace>(mesh, constrained_domain);
    }
    else {
      std::cout << "Unrecognized velocity element: " << u_el << std::endl;
      exit(0);
    }
    // Pressure
    if (p_el == "P1"){
      p_space = std::make_shared<P1_3::FunctionSpace>(mesh, constrained_domain);
    }
    else if (u_el == "P2"){
      p_space = std::make_shared<P2_3::FunctionSpace>(mesh, constrained_domain);
    }
    else if (u_el == "P3"){
      p_space = std::make_shared<P3_3::FunctionSpace>(mesh, constrained_domain);
    }
    else {
      std::cout << "Unrecognized pressure element: " << u_el << std::endl;
      exit(0);
    }
    break;
  default:
    std::cout << "Unsupported dimensionality" << std::endl;
    exit(0);
  }

  u_prev_ = std::make_shared<dolfin::Function>(u_space);
  u_next_ = std::make_shared<dolfin::Function>(u_space);
  p_prev_ = std::make_shared<dolfin::Function>(p_space);
  p_next_ = std::make_shared<dolfin::Function>(p_space);

  //std::cout << "GOT THIS FAR" << std::endl;
}

void XDMFInterpol::update(const double t){
  StampPair sp = ts.get(t);
  // std::cout << sp.prev.filename << " " << sp.next.filename << std::endl;

  if (!is_initialized || t_prev != sp.prev.t || t_next != sp.next.t){
    std::cout << "Prev: Timestep = " << sp.prev.t << ", filename = " << sp.prev.filename << std::endl;
    dolfin::HDF5File prevfile(MPI_COMM_WORLD, get_folder() + "/" + sp.prev.filename, "r");
    prevfile.read(*u_prev_, "u");
    prevfile.read(*p_prev_, "p");

    std::cout << "Next: Timestep = " << sp.next.t << ", filename = " << sp.next.filename << std::endl;
    dolfin::HDF5File nextfile(MPI_COMM_WORLD, get_folder() + "/" + sp.prev.filename, "r");
    nextfile.read(*u_next_, "u");
    nextfile.read(*p_next_, "p");

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  alpha_t = sp.weight_next(t);
}

void XDMFInterpol::probe(const Vector3d &x){
  //const double* _x = x.data();
  dolfin::Array<double> x_loc(dim);
  for (Uint i=0; i<dim; ++i){
    if (periodic[i]){
      x_loc[i] = x_min[i] + modulox(x[i]-x_min[i], x_max[i]-x_min[i]);
    }
    else {
      x_loc[i] = x[i];
    }
  }

  //std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
  //std::cout << x_loc[0] << " " << x_loc[1] << " " << x_loc[2] << std::endl;

  double* _x = x_loc.data();
  const dolfin::Point point(dim, _x);
  // index of cell containing point
  unsigned int id
    = mesh->bounding_box_tree()->compute_first_entity_collision(point);

  inside = (id != std::numeric_limits<unsigned int>::max());
  if (inside){
    const dolfin::Cell dolfin_cell(*mesh, id);
    ufc::cell ufc_cell;
    dolfin_cell.get_cell_data(ufc_cell);

    // dolfin::Array<double> vals(3);
    // u_prev_->eval(vals, x_loc, dolfin_cell, ufc_cell);
    // std::cout << "U_prev_a = " << vals[0] << " " << vals[1] << " " << vals[2] << std::endl;

    const dolfin::FiniteElement u_element = *u_space->element();
    const dolfin::FiniteElement p_element = *p_space->element();

    Uint u_dim = u_element.space_dimension();
    Uint p_dim = p_element.space_dimension();

    //std::cout << "u_element.space_dimension() = " << u_dim << std::endl;
    //std::cout << "p_element.space_dimension() = " << p_dim << std::endl;

    // Work vectors for expansion coefficients
    std::vector<double> u_prev_coefficients(u_dim);
    std::vector<double> u_next_coefficients(u_dim);
    std::vector<double> p_prev_coefficients(p_dim);
    std::vector<double> p_next_coefficients(p_dim);

    // Cell coordinates
    std::vector<double> coordinate_dofs;
    dolfin_cell.get_coordinate_dofs(coordinate_dofs);
    // std::cout << "coordinate_dofs = ";
    // for (int i=0; i < coordinate_dofs.size(); ++i){
    //   std::cout << coordinate_dofs[i] << " ";
    // }
    // std::cout << std::endl;

    u_prev_->restrict(u_prev_coefficients.data(), u_element, dolfin_cell,
                      coordinate_dofs.data(), ufc_cell);
    u_next_->restrict(u_next_coefficients.data(), u_element, dolfin_cell,
                      coordinate_dofs.data(), ufc_cell);

    p_prev_->restrict(p_prev_coefficients.data(), p_element, dolfin_cell,
                      coordinate_dofs.data(), ufc_cell);
    p_next_->restrict(p_next_coefficients.data(), p_element, dolfin_cell,
                      coordinate_dofs.data(), ufc_cell);

    // std::cout << "u_prev_coeff: " << std::endl;
    // for (int i=0; i<dim; ++i)
    //   std::cout << u_prev_coefficients[i] << " ";
    // std::cout << std::endl;

    Vector3d U_prev = {0., 0., 0.};
    Vector3d U_next = {0., 0., 0.};
    double P_prev = 0.;
    double P_next = 0.;

    Matrix3d gradU_prev;
    gradU_prev << 0., 0., 0., 0., 0., 0., 0., 0., 0.;
    Matrix3d gradU_next;
    gradU_next << 0., 0., 0., 0., 0., 0., 0., 0., 0.;
    Vector3d gradP_prev = {0., 0., 0.};
    Vector3d gradP_next = {0., 0., 0.};

    // Work vector for basis
    std::vector<double> u_basis(dim);
    std::vector<double> p_basis(1);
    std::vector<double> gradu_basis(dim*dim);
    std::vector<double> gradp_basis(dim);

    for (Uint i=0; i<u_dim; ++i){
      u_element.evaluate_basis(i, u_basis.data(), _x,
                               coordinate_dofs.data(),
                               ufc_cell.orientation);
      for (Uint j=0; j<dim; ++j){
        U_prev[j] += u_prev_coefficients[i]*u_basis[j];
        U_next[j] += u_next_coefficients[i]*u_basis[j];
      }
    }
    for (Uint i=0; i<p_dim; ++i){
      p_element.evaluate_basis(i, p_basis.data(), _x,
                               coordinate_dofs.data(),
                               ufc_cell.orientation);
      P_prev += p_prev_coefficients[i]*p_basis[0];
      P_next += p_next_coefficients[i]*p_basis[0];
    }

    if (this->int_order > 1){
      for (Uint i=0; i<u_dim; ++i){
        u_element.evaluate_basis_derivatives(i, 1,
                                             gradu_basis.data(), _x,
                                             coordinate_dofs.data(),
                                             ufc_cell.orientation);
        for (Uint j=0; j<dim; ++j){
          for (Uint k=0; k<dim; ++k){
            gradU_prev(j, k) += u_prev_coefficients[i]*gradu_basis[dim*j+k];
            gradU_next(j, k) += u_next_coefficients[i]*gradu_basis[dim*j+k];
          }
        }
      }
      /* // GradP not used for anything...
      for (Uint i=0; i<p_dim; ++i){
        p_element.evaluate_basis_derivatives(i, 1,
                                             gradp_basis.data(), _x,
                                             coordinate_dofs.data(),
                                             ufc_cell.orientation);
        for (Uint j=0; j<dim; ++j){
          gradP_prev[j] += p_prev_coefficients[i]*gradp_basis[j];
          gradP_next[j] += p_next_coefficients[i]*gradp_basis[j];
        }
      }
      */
    }

    U = alpha_t * U_next + (1-alpha_t) * U_prev;
    A = (U_next-U_prev)/(t_next-t_prev);
    P = alpha_t * P_next + (1-alpha_t) * P_prev;
    if (this->int_order > 1){
      gradU = alpha_t * gradU_next + (1-alpha_t) * gradU_prev;
      gradA = (gradU_next-gradU_prev)/(t_next-t_prev);
    }
    //std::cout << x << std::endl;
    //std::cout << U << std::endl;
    //std::cout << A << std::endl;
    //std::cout << "U_prev_b = " << U_prev[0] << " " << U_prev[1] << " " << U_prev[2] << std::endl;
  }
}

#endif
