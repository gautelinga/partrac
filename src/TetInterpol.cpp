#include "TetInterpol.hpp"
#include "Timestamps.hpp"
#include "H5Cpp.h"
#include <boost/algorithm/string.hpp>
#include <cassert>
#include "dolfin_elements/P1_3.h"
#include "dolfin_elements/vP2_3.h"
#include "PBC.hpp"

using namespace H5;
//using namespace std;
//using namespace dolfin;


TetInterpol::TetInterpol(const std::string infilename)
  : Interpol(infilename)
{
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
      boost::algorithm::trim(key);
      boost::algorithm::trim(val);
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
  //dolfin::BasicFile xdmff(MPI_COMM_WORLD, xdmffname);
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

  // Precompute all tets Taylor-Hood P2-P1
  // FIXME compute on the fly and save
  tets_.resize(mesh->num_cells());

  for (std::size_t i = 0; i < mesh->num_cells(); ++i)
  {
    tets_[i] = Tet(dolfin::Cell(*mesh, i));
  }

  auto constrained_domain = std::make_shared<PBC>(periodic, x_min, x_max, dim);
  u_space = std::make_shared<vP2_3::FunctionSpace>(mesh, constrained_domain);
  p_space = std::make_shared<P1_3::FunctionSpace>(mesh, constrained_domain);

  u_prev_ = std::make_shared<dolfin::Function>(u_space);
  u_next_ = std::make_shared<dolfin::Function>(u_space);
  p_prev_ = std::make_shared<dolfin::Function>(p_space);
  p_next_ = std::make_shared<dolfin::Function>(p_space);

}

void TetInterpol::update(const double t)
{

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

void TetInterpol::probe(const Vector3d &x)
{
  dolfin::Array<double> x_loc(dim);
  for (std::size_t i=0; i<dim; ++i){
    if (periodic[i]){
      x_loc[i] = x_min[i] + modulox(x[i]-x_min[i], x_max[i]-x_min[i]);
    }
    else {
      x_loc[i] = x[i];
    }
  }

  // Index of cell containing point
  const dolfin::Point point(dim, x_loc.data());
  unsigned int id
    = mesh->bounding_box_tree()->compute_first_entity_collision(point);

  inside = (id != std::numeric_limits<unsigned int>::max());
  if (inside)
  {
    // Compute P2-P1 basis at x
    double r,s,t,u;
    tets_[id].xyz2bary(x_loc[0], x_loc[1], x_loc[2], r,s,t,u);
    std::array<double, 10> N10;
    std::array<double, 4> N4;
    tets_[id].quadbasis(r,s,t,u, N10);
    tets_[id].linearbasis(r,s,t,u, N4);

    // Restrict solution to cell
    const dolfin::Cell dolfin_cell(*mesh, id);
    ufc::cell ufc_cell;
    dolfin_cell.get_cell_data(ufc_cell);
    std::vector<double> coordinate_dofs;
    dolfin_cell.get_coordinate_dofs(coordinate_dofs);

    const dolfin::FiniteElement u_element = *u_space->element();
    const dolfin::FiniteElement p_element = *p_space->element();
    std::size_t u_dim = u_element.space_dimension();
    std::size_t p_dim = p_element.space_dimension();
    std::vector<double> u_prev_coefficients(u_dim);
    std::vector<double> u_next_coefficients(u_dim);
    std::vector<double> p_prev_coefficients(p_dim);
    std::vector<double> p_next_coefficients(p_dim);

    u_prev_->restrict(u_prev_coefficients.data(), u_element, dolfin_cell,
                      coordinate_dofs.data(), ufc_cell);
    u_next_->restrict(u_next_coefficients.data(), u_element, dolfin_cell,
                      coordinate_dofs.data(), ufc_cell);
    p_prev_->restrict(p_prev_coefficients.data(), p_element, dolfin_cell,
                      coordinate_dofs.data(), ufc_cell);
    p_next_->restrict(p_next_coefficients.data(), p_element, dolfin_cell,
                      coordinate_dofs.data(), ufc_cell);

    // Evaluate
    double P_prev = std::inner_product(N4.begin(), N4.end(), p_prev_coefficients.begin(), 0.0);
    double P_next = std::inner_product(N4.begin(), N4.end(), p_next_coefficients.begin(), 0.0);

    Vector3d U_prev = {0., 0., 0.};
    Vector3d U_next = {0., 0., 0.};
    // for (std::size_t i = 0, idx = 0; i < 10; ++i)
    //   for (std::size_t j = 0; j < 3; ++j, ++idx)
    //   {
    //     U_prev[j] += u_prev_coefficients[idx]*N10[i];
    //     U_next[j] += u_next_coefficients[idx]*N10[i];
    //   }
    for (std::size_t j = 0, idx = 0; j < 3; ++j)
      for (std::size_t i = 0; i < 10; ++i, ++idx)
      {
        U_prev[j] += u_prev_coefficients[idx]*N10[i];
        U_next[j] += u_next_coefficients[idx]*N10[i];
      }
    assert(true);
    assert(false);
    // Update
    U = alpha_t * U_next + (1-alpha_t) * U_prev;
    A = (U_next-U_prev)/(t_next-t_prev);
    P = alpha_t * P_next + (1-alpha_t) * P_prev;
    assert(this->int_order == 1);
    // if (this->int_order > 1){
    //   gradU = alpha_t * gradU_next + (1-alpha_t) * gradU_prev;
    //   gradA = (gradU_next-gradU_prev)/(t_next-t_prev);
    // }
    //std::cout<<x[0]<<' '<<x[1]<<' '<<x[2]<<"   "<<x_loc[0]<<' '<<x_loc[1]<<' '<<x_loc[2]<<"    "<<id<<' '<<inside<<"   "<<P_prev<<' '<<P_next<<"    "<<U_prev[0]<<' '<<U_prev[1]<<' '<<U_prev[2]<<"   "<<U_next[0]<<' '<<U_next[1]<<' '<<U_next[2]<<"   "<<alpha_t<<std::endl;
  }

}
