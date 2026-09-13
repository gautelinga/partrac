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
  Uint get_nx() { return n[0]; };
  Uint get_ny() { return n[1]; };
  Uint get_nz() { return n[2]; };
  double get_t_min() { return ts.get_t_min(); };
  double get_t_max() { return ts.get_t_max(); };
  using Interpol::locate;
  using Interpol::evaluate;
protected:
  void probe_space(const Vector3d &x);
  void probe_grad();
  Timestamps ts;
  double t_prev = 0.;
  double t_next = 0.;
  double alpha_t;

  Uint n[3] = {0, 0, 0};
  Vector3d dx;

  GridBlock<int> solid_;
  GridBlock<double> fields_;
  Grid3<int> isSolid;
  Grid3<double> ux_prev, uy_prev, uz_prev;
  Grid3<double> ux_next, uy_next, uz_next;
  Grid3<double> rho_prev, rho_next;
  Grid3<double> p_prev, p_next;

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
    exit(1);
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
  solid_.resize(n[0], n[1], n[2], 1);
  isSolid = solid_.field(0);
  fields_.resize(n[0], n[1], n[2], 10);
  Uint f = 0;
  for (Grid3<double>* g : {&ux_prev, &ux_next, &uy_prev, &uy_next, &uz_prev, &uz_next,
                           &rho_prev, &rho_next, &p_prev, &p_next})
    *g = fields_.field(f++);
  load_int_field(solid_file, isSolid, "is_solid", n[0], n[1], n[2]);
}

void StructuredConstInterpol::update(const double t){
  StampPair sp = ts.get(t);

  if (!is_initialized || t_prev != sp.prev.t || t_next != sp.next.t){
    if (is_initialized && t_next == sp.prev.t){
      std::swap(ux_prev, ux_next);
      std::swap(uy_prev, uy_next);
      if (!ignore_uz)
        std::swap(uz_prev, uz_next);
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



// Interpolate in space and time and enforce BCs








#endif
