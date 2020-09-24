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

#include "dolfin_elements/vP2_3.h"
#include "dolfin_elements/P1_3.h"

#ifndef __XDMFINTERPOL_HPP
#define __XDMFINTERPOL_HPP

using namespace H5;
//using namespace std;
//using namespace dolfin;

class XDMFInterpol : public Interpol {
public:
  XDMFInterpol(const std::string infilename);
  void update(const double t);
  void probe(const Vector3d &x);
  bool inside_domain();
  double get_ux();
  double get_uy();
  double get_uz();
  double get_ax() { return 0.; };
  double get_ay() { return 0.;};
  double get_az() { return 0.;};
  double get_t_min() { return ts.get_t_min(); };
  double get_t_max() { return ts.get_t_max(); };
  double get_rho() { return 0.;};
  double get_p() {return 0.;};
  double get_uxx() {return 0.;};
  double get_uxy() {return 0.;};
  double get_uxz() {return 0.;};
  double get_uyx() {return 0.;};
  double get_uyy() {return 0.;};
  double get_uyz() {return 0.;};
  double get_uzx() {return 0.;};
  double get_uzy() {return 0.;};
  double get_uzz() {return 0.;};
  Matrix3d get_grada() {
    Matrix3d aa;
    aa <<
      0., 0., 0.,
      0., 0., 0.,
      0., 0., 0.;
    return aa;
  };
protected:
  Timestamps ts;
  double t_prev = 0.;
  double t_next = 0.;
  double alpha_t;

  double Ux = 0.;
  double Uy = 0.;
  double Uz = 0.;
  double Ax = 0.;
  double Ay = 0.;
  double Az = 0.;

  std::map<std::string, std::string> dolfin_params;
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

  std::cout << dolfin_params["velocity"] << std::endl;

  //std::string xdmffname = get_folder() + "/" + dolfin_params["velocity"];
  //dolfin::XDMFFile xdmff(MPI_COMM_WORLD, xdmffname);
  //std::string h5filename = get_folder() + "/" + ts.get(0).prev.filename;
  //std::cout << h5filename << std::endl;
  //MPI_Comm mpi_comm;

  std::string meshfilename = get_folder() + "/" + dolfin_params["mesh"];
  dolfin::HDF5File meshfile(MPI_COMM_WORLD, meshfilename, "r");

  dolfin::Mesh mesh;
  meshfile.read(mesh, "mesh", false);

  auto mesh2 = std::make_shared<dolfin::Mesh>(mesh);

  std::vector<double> xx = mesh.coordinates();
  dolfin::MeshGeometry mesh_geom = mesh.geometry();
  dolfin::MeshTopology mesh_topo = mesh.topology();

  //std::cout << "Degree: " << mesh_geom.degree() << std::endl;
  std::cout << "Dim:    " << mesh_geom.dim() << std::endl;

  std::cout << xx[0] << " " << xx.size() << std::endl;

  auto V = std::make_shared<vP2_3::FunctionSpace>(mesh2);
  auto P = std::make_shared<P1_3::FunctionSpace>(mesh2);
  auto u_ = std::make_shared<dolfin::Function>(V);
  auto p_ = std::make_shared<dolfin::Function>(P);

  //dolfin::MeshFunction<double> ux(mesh2, 1);
  //xdmff.read(ux);

  std::cout << "GOT THIS FAR" << std::endl;
}

void XDMFInterpol::update(const double t){
  StampPair sp = ts.get(t);
  // std::cout << sp.prev.filename << " " << sp.next.filename << std::endl;

  std::cout << "SGDSG" << std::endl;

  if (!is_initialized || t_prev != sp.prev.t || t_next != sp.next.t){
    std::cout << "Prev: Timestep = " << sp.prev.t << ", filename = " << sp.prev.filename << std::endl;
    dolfin::HDF5File prevfile(MPI_COMM_WORLD, get_folder() + "/" + sp.prev.filename, "r");
    exit(0);

    std::cout << "Next: Timestep = " << sp.next.t << ", filename = " << sp.next.filename << std::endl;

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  alpha_t = sp.weight_next(t);
}

void XDMFInterpol::probe(const Vector3d &x){
}

bool XDMFInterpol::inside_domain(){
  return true;
}

// Interpolate in space and time and enforce BCs
double XDMFInterpol::get_ux(){
  return 0.;
}

double XDMFInterpol::get_uy(){
  return 0.;
}

double XDMFInterpol::get_uz(){
  return 0.;
}

#endif
 
