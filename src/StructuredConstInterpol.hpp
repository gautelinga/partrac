#ifndef __STRUCTUREDCONSTINTERPOL_HPP
#define __STRUCTUREDCONSTINTERPOL_HPP

#include "Interpol.hpp"
#include "StructuredInterpol.hpp"
#include "Timestamps.hpp"

#include "H5Cpp.h"

class StructuredConstInterpol : public Interpol {
public:
  StructuredConstInterpol(const std::string& infilename);
  void update(const double t);
  void probe(const Vector3d &x, const double t);
  void probe(const Vector3d &x, const double t, int &cell_id) { probe(x, t); };
  bool inside_domain() const;
  Uint get_nx() { return n[0]; };
  Uint get_ny() { return n[1]; };
  Uint get_nz() { return n[2]; };
  double get_ux();
  double get_uy();
  double get_uz();
  double get_ax();
  double get_ay();
  double get_az();
  double get_t_min() { return ts.get_t_min(); };
  double get_t_max() { return ts.get_t_max(); };
  double get_rho();
  double get_p();
  double get_uxx() { return 0.; };
  double get_uxy() { return 0.; };
  double get_uxz() { return 0.; };
  double get_uyx() { return 0.; };
  double get_uyy() { return 0.; };
  double get_uyz() { return 0.; };
  double get_uzx() { return 0.; };
  double get_uzy() { return 0.; };
  double get_uzz() { return 0.; };
  Matrix3d get_grada()  { Matrix3d gradA; gradA << 0., 0., 0., 0., 0., 0., 0., 0., 0.; return gradA; };
  void probe(const Vector3d &x){ probe(x, this->t_update); };
protected:
  void probe_space(const Vector3d &x);
  void probe_grad();
  Timestamps ts;
  double t_prev = 0.;
  double t_next = 0.;
  double alpha_t;

  Uint n[3] = {0, 0, 0};
  Vector3d dx;

  int*** isSolid;
  double*** levelZ;
  double*** ux_prev;
  double*** uy_prev;
  double*** uz_prev;
  double*** ux_next;
  double*** uy_next;
  double*** uz_next;
  double*** rho_prev;
  double*** rho_next;
  double*** p_prev;
  double*** p_next;

  Uint ind_pc[3] = {0, 0, 0};  // piecewise constant intp

  std::map<std::string, std::string> felbm_params;

  double Ux = 0.;
  double Uy = 0.;
  double Uz = 0.;
  double Ax = 0.;
  double Ay = 0.;
  double Az = 0.;
  bool inside = false;

  bool ignore_density = false;
  bool ignore_pressure = false;
  bool ignore_uz = false;
};

StructuredConstInterpol::StructuredConstInterpol(const std::string& infilename) : Interpol(infilename) {
  std::ifstream input(infilename);
  if (!input){
    std::cout << "File " << infilename <<" doesn't exist." << std::endl;
    exit(0);
  }
  // Default params
  felbm_params["timestamps"] = "timestamps.dat";
  felbm_params["is_solid_file"] = "output_is_solid.h5";
  felbm_params["boundary_mode"] = "rounded";

  felbm_params["ignore_pressure"] = "false";
  felbm_params["ignore_density"] = "false";
  felbm_params["ignore_uz"] = "false";

  // Overload from file
  size_t found;
  std::string key, val;
  for (std::string line; getline(input, line); ){
    found = line.find('=');
    if (found != std::string::npos){
      key = line.substr(0, found);
      val = line.substr(found+1);
      boost::trim(key);
      boost::trim(val);
      felbm_params[key] = val;
      //std::cout << "--> " << key << ": " << val << std::endl;
    }
  }
  std::cout << "Chosen parameters:" << std::endl;
  for (auto & k : felbm_params) {
    std::cout << k.first << ": " << k.second << std::endl;
  }

  if (felbm_params["ignore_pressure"] == "true" || felbm_params["ignore_pressure"] == "True"){
    ignore_pressure = true;
  }
  if (felbm_params["ignore_density"] == "true" || felbm_params["ignore_density"] == "True"){
    ignore_density = true;
  }
  if (felbm_params["ignore_uz"] == "true" || felbm_params["ignore_uz"] == "True"){
    ignore_uz = true;
  }

  std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));
  ts.initialize(get_folder() + "/" + felbm_params["timestamps"]);

  std::string solid_filename = get_folder() + "/" + felbm_params["is_solid_file"];
  verify_file_exists(solid_filename);

  H5::H5File solid_file(solid_filename, H5F_ACC_RDONLY);
  H5::DataSet dset_solid = solid_file.openDataSet("is_solid");
  H5::DataSpace dspace_solid = dset_solid.getSpace();

  hsize_t dims[3];
  dspace_solid.getSimpleExtentDims(dims, NULL);
  for (Uint i=0; i<3; ++i)
    n[i] = dims[i];
  x_min << 0., 0., 0.;
  x_max << n[0], n[1], n[2];

  dx << this->get_Lx()/n[0], this->get_Ly()/n[1], this->get_Lz()/n[2];

  // Create arrays
  isSolid = new int**[n[0]];
  levelZ = new double**[n[0]];
  ux_prev = new double**[n[0]];
  uy_prev = new double**[n[0]];
  uz_prev = new double**[n[0]];
  ux_next = new double**[n[0]];
  uy_next = new double**[n[0]];
  uz_next = new double**[n[0]];
  rho_prev = new double**[n[0]];
  rho_next = new double**[n[0]];
  p_prev = new double**[n[0]];
  p_next = new double**[n[0]];
  for (Uint ix=0; ix<n[0]; ++ix){
    isSolid[ix] = new int*[n[1]];
    levelZ[ix] = new double*[n[1]];
    ux_prev[ix] = new double*[n[1]];
    uy_prev[ix] = new double*[n[1]];
    uz_prev[ix] = new double*[n[1]];
    ux_next[ix] = new double*[n[1]];
    uy_next[ix] = new double*[n[1]];
    uz_next[ix] = new double*[n[1]];
    rho_prev[ix] = new double*[n[1]];
    rho_next[ix] = new double*[n[1]];
    p_prev[ix] = new double*[n[1]];
    p_next[ix] = new double*[n[1]];
    for (Uint iy=0; iy<n[1]; ++iy){
      isSolid[ix][iy] = new int[n[2]];
      levelZ[ix][iy] = new double[n[2]];
      ux_prev[ix][iy] = new double[n[2]];
      uy_prev[ix][iy] = new double[n[2]];
      uz_prev[ix][iy] = new double[n[2]];
      ux_next[ix][iy] = new double[n[2]];
      uy_next[ix][iy] = new double[n[2]];
      uz_next[ix][iy] = new double[n[2]];
      rho_prev[ix][iy] = new double[n[2]];
      rho_next[ix][iy] = new double[n[2]];
      p_prev[ix][iy] = new double[n[2]];
      p_next[ix][iy] = new double[n[2]];
    }
  }
  load_int_field(solid_file, isSolid, "is_solid", n[0], n[1], n[2]);
}

void StructuredConstInterpol::update(const double t){
  StampPair sp = ts.get(t);

  if (!is_initialized || t_prev != sp.prev.t || t_next != sp.next.t){
    if (is_initialized && t_next == sp.prev.t){
      std::swap(ux_prev, ux_next);
      std::swap(uy_prev, uy_next);
      if (!ignore_uz)
        std::swap(uz_prev, uz_prev);
      if (!ignore_density)
        std::swap(rho_prev, rho_next);
      if (!ignore_pressure)
        std::swap(p_prev, p_next);
    }
    else {
      std::cout << "Previous: Timestep = " << sp.prev.t << ", filename = " << sp.prev.filename << std::endl;
      load_h5(folder + "/" + sp.prev.filename,
              ux_prev, uy_prev, uz_prev, rho_prev, p_prev,
              n[0], n[1], n[2],
              verbose, ignore_density, ignore_pressure, ignore_uz);
    }

    std::cout << "Next: Timestep = " << sp.next.t << ", filename = " << sp.next.filename << std::endl;
    load_h5(folder + "/" + sp.next.filename,
            ux_next, uy_next, uz_next, rho_next, p_next,
            n[0], n[1], n[2], 
            verbose, ignore_density, ignore_pressure, ignore_uz);

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  // alpha_t = sp.weight_next(t);
  t_update = t;
}

void StructuredConstInterpol::probe(const Vector3d &x, const double t){
  assert(t <= t_next && t >= t_prev);
  alpha_t = (t-t_prev)/(t_next-t_prev);

  for (Uint i=0; i<3; ++i){
    ind_pc[i] = imodulo(round(x[i]/dx[i]), n[i]);
  }

  // Precompute inside factor
  inside = !isSolid[ind_pc[0]][ind_pc[1]][ind_pc[2]];

  // Precompute velocities
  double Ux_prev = ux_prev[ind_pc[0]][ind_pc[1]][ind_pc[2]];
  double Ux_next = ux_next[ind_pc[0]][ind_pc[1]][ind_pc[2]];
  Ux = alpha_t * Ux_next + (1-alpha_t) * Ux_prev;
  Ax = (Ux_next-Ux_prev)/(t_next-t_prev);

  double Uy_prev = uy_prev[ind_pc[0]][ind_pc[1]][ind_pc[2]];
  double Uy_next = uy_next[ind_pc[0]][ind_pc[1]][ind_pc[2]];
  Uy = alpha_t * Uy_next + (1-alpha_t) * Uy_prev;
  Ay = (Uy_next-Uy_prev)/(t_next-t_prev);

  if (!ignore_uz){
    double Uz_prev = uz_prev[ind_pc[0]][ind_pc[1]][ind_pc[2]];
    double Uz_next = uz_next[ind_pc[0]][ind_pc[1]][ind_pc[2]];
    Uz = alpha_t * Uz_next + (1-alpha_t) * Uz_prev;
    Az = (Uz_next-Uz_prev)/(t_next-t_prev);
  }
}

bool StructuredConstInterpol::inside_domain() const {
  return inside;
}

// Interpolate in space and time and enforce BCs
double StructuredConstInterpol::get_ux(){
  return Ux;
}

double StructuredConstInterpol::get_uy(){
  return Uy;
}

double StructuredConstInterpol::get_uz(){
  return Uz;
}

double StructuredConstInterpol::get_ax(){
  return Ax;
}

double StructuredConstInterpol::get_ay(){
  return Ay;
}

double StructuredConstInterpol::get_az(){
  return Az;
}

double StructuredConstInterpol::get_rho(){
  double Rho_prev = rho_prev[ind_pc[0]][ind_pc[1]][ind_pc[2]];
  double Rho_next = rho_next[ind_pc[0]][ind_pc[1]][ind_pc[2]];
  return alpha_t * Rho_next + (1-alpha_t) * Rho_prev;
}

double StructuredConstInterpol::get_p(){
  double P_prev = p_prev[ind_pc[0]][ind_pc[1]][ind_pc[2]];
  double P_next = p_next[ind_pc[0]][ind_pc[1]][ind_pc[2]];
  return alpha_t * P_next + (1-alpha_t) * P_prev;
}

#endif
