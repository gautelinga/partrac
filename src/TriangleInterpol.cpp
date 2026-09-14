#ifdef USE_DOLFIN
#include "TriangleInterpol.hpp"
#include <array>
#include <numeric>
#include "Timestamps.hpp"
//#include "H5Cpp.h"
#include <boost/algorithm/string.hpp>
#include <cassert>
#include "dolfin_elements/P1_2.h"
#include "dolfin_elements/P2_2.h"
#include "dolfin_elements/vP1_2.h"
#include "dolfin_elements/vP2_2.h"
#include "PeriodicBC.hpp"
#include "dolfin_helpers.hpp"

//using namespace H5;

TriangleInterpol::TriangleInterpol(const std::string& infilename)
  : Interpol(infilename)
{
  std::ifstream input(infilename);
  if (!input){
    std::cout << "File " << infilename <<" doesn't exist." << std::endl;
    exit(1);
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
  if (dolfin_params["ignore_pressure"] == "true"){
    include_pressure = false;
  }
  std::string meshfilename = get_folder() + "/" + dolfin_params["mesh"];
  dolfin::HDF5File meshfile(MPI_COMM_WORLD, meshfilename, "r");

  dolfin::Mesh mesh_in;
  meshfile.read(mesh_in, "mesh", false);

  mesh = std::make_shared<dolfin::Mesh>(mesh_in);
  dim = mesh->geometry().dim();
  mesh->init();
  mesh->bounding_box_tree();

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

  // Precompute all triangles Taylor-Hood P2-P1
  // FIXME compute on the fly and save
  triangles_.resize(mesh->num_cells());
  dolfin_cells_.resize(mesh->num_cells());
  cell2cells_.resize(mesh->num_cells());

  for (std::size_t i = 0; i < mesh->num_cells(); ++i)
  {
    dolfin::Cell dolfin_cell(*mesh, i);
    triangles_[i] = Triangle(dolfin_cell);
    dolfin_cells_[i] = dolfin_cell;
  }
  // Build cell neighbour list for lookup speed
  build_neighbor_list(cell2cells_, mesh, dolfin_cells_);

  std::cout << "Built neighbour list" << std::endl;
  
  auto constrained_domain = std::make_shared<PeriodicBC>(periodic, x_min, x_max, dim);
  std::cout << "Made periodic domain." << std::endl;

  std::string u_el = dolfin_params["velocity_space"];
  std::string p_el = dolfin_params["pressure_space"];
  
  // Velocity
  if (u_el == "P1"){
    u_space_ = std::make_shared<vP1_2::FunctionSpace>(mesh, constrained_domain);
    ncoeffs_u = 3;
  }
  else if (u_el == "P2"){
    u_space_ = std::make_shared<vP2_2::FunctionSpace>(mesh, constrained_domain);
    ncoeffs_u = 6;
  }
  else {
    std::cout << "Unrecognized velocity element: " << u_el << std::endl;
    exit(1);
  }

  // Pressure
  if (include_pressure){
    if (p_el == "P1"){
      p_space_ = std::make_shared<P1_2::FunctionSpace>(mesh, constrained_domain);
      ncoeffs_p = 3;
    }
    else if (p_el == "P2"){
      p_space_ = std::make_shared<P2_2::FunctionSpace>(mesh, constrained_domain);
      ncoeffs_p = 6;
    }
    else {
      std::cout << "Unrecognized pressure element: " << p_el << std::endl;
      exit(1);
    }
  }
  else {
    std::cout << "Note: Ignoring pressure." << std::endl;
  }

  u_prev_ = std::make_shared<dolfin::Function>(u_space_);
  u_next_ = std::make_shared<dolfin::Function>(u_space_);

  u_prev_coefficients_.resize(dim*ncoeffs_u);
  u_next_coefficients_.resize(dim*ncoeffs_u);

  Nu_.resize(ncoeffs_u);

  Nux_.resize(ncoeffs_u);
  Nuy_.resize(ncoeffs_u);

  if (include_pressure){
    p_prev_ = std::make_shared<dolfin::Function>(p_space_);
    p_next_ = std::make_shared<dolfin::Function>(p_space_);

    p_prev_coefficients_.resize(ncoeffs_p);
    p_next_coefficients_.resize(ncoeffs_p);
    
    Np_.resize(ncoeffs_p);
  }

  // Precompute dofs of all cells
  u_dofs_.build(*u_space_->dofmap(), dolfin_cells_, "TriangleInterpol");
  if (include_pressure)
    p_dofs_.build(*p_space_->dofmap(), dolfin_cells_, "TriangleInterpol");

  check_dofs_fit(ncoeffs_u, ncoeffs_p, Triangle::n_dofs_max, "TriangleInterpol");

  std::cout << "Setting max threads: " << omp_get_max_threads() << std::endl;

  found_.resize(omp_get_max_threads());
}

void TriangleInterpol::update(const double t)
{
  // TODO: swapping!

  StampPair sp = ts.get(t);
  // std::cout << sp.prev.filename << " " << sp.next.filename << std::endl;

  if (!is_initialized || t_prev != sp.prev.t || t_next != sp.next.t){
    // Swap if possible
    if (is_initialized && t_next == sp.prev.t){
      std::cout << "Prev: Timestep = " << sp.prev.t << ", swapping... " << std::endl;
      u_prev_data_.swap(u_next_data_);
      if (include_pressure)
        p_prev_data_.swap(p_next_data_);
    }
    else {
      std::cout << "Prev: Timestep = " << sp.prev.t << ", filename = " << sp.prev.filename << std::endl;
      dolfin::HDF5File prevfile(MPI_COMM_WORLD, get_folder() + "/" + sp.prev.filename, "r");
      prevfile.read(*u_prev_, "u");
      u_prev_->vector()->get_local(u_prev_data_);
      if (include_pressure){
        prevfile.read(*p_prev_, "p");
        p_prev_->vector()->get_local(p_prev_data_);
      }
    }

    std::cout << "Next: Timestep = " << sp.next.t << ", filename = " << sp.next.filename << std::endl;
    // Single stamp: copy prev
    if (sp.next.filename == sp.prev.filename){
      u_next_data_ = u_prev_data_;
      if (include_pressure)
        p_next_data_ = p_prev_data_;
    }
    else {
      dolfin::HDF5File nextfile(MPI_COMM_WORLD, get_folder() + "/" + sp.next.filename, "r");
      nextfile.read(*u_next_, "u");
      u_next_->vector()->get_local(u_next_data_);
      if (include_pressure){
        nextfile.read(*p_next_, "p");
        p_next_->vector()->get_local(p_next_data_);
      }
    }

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  //alpha_t = sp.weight_next(t);
  t_update = t;
}




Vector3d TriangleInterpol::_modx(const Vector3d &x){
  Vector3d x_loc = x;
  for (std::size_t i=0; i<dim; ++i)
    if (periodic[i])
      x_loc[i] = x_min[i] + modulox(x[i]-x_min[i], x_max[i]-x_min[i]);
  return x_loc;
}

bool TriangleInterpol::locate(const Vector3d &x, const double t, int& id_prev)
{
  assert(t <= t_next && t >= t_prev);
  const Vector3d xx = _modx(x);
  return locate_in_cells(triangles_, cell2cells_, *mesh, dim, xx, id_prev,
                         found_);
}

void TriangleInterpol::evaluate(const Vector3d &x, const double tin, const int id, PointValues& fields)
{
  const Vector3d x_loc = _modx(x);

  const double _alpha_t = stamp_weight(tin, t_prev, t_next);

  // Compute Pk-Pl basis at x
  double r, s, u;
  triangles_[id].xy2bary(x_loc[0], x_loc[1], r, s, u);

  std::array<double, Triangle::n_dofs_max> _Nu_{};
  std::array<double, Triangle::n_dofs_max> _Np_{};
  std::array<double, Triangle::n_dofs_max> _Nux_{};
  std::array<double, Triangle::n_dofs_max> _Nuy_{};

  if (ncoeffs_u == 3){
    triangles_[id].linearbasis(r, s, u, _Nu_.data());
  }
  else if (ncoeffs_u == 6){
    triangles_[id].quadbasis(r, s, u, _Nu_.data());
  }
  else {
    std::cout << "Unrecognized ncoeffs_u = " << ncoeffs_u << std::endl;
    exit(1);
  }
  if (include_pressure){
    if (ncoeffs_p == 3){
      triangles_[id].linearbasis(r, s, u, _Np_.data());
    }
    else if (ncoeffs_p == 6){
      triangles_[id].quadbasis(r, s, u, _Np_.data());
    }
    else {
      std::cout << "Unrecognized ncoeffs_p = " << ncoeffs_p << std::endl;
      exit(1);
    }
  }

  std::array<double, Triangle::n_dofs_max*3> u_prev_block{};
  std::array<double, Triangle::n_dofs_max*3> u_next_block{};
  std::array<double, Triangle::n_dofs_max> p_prev_block{};
  std::array<double, Triangle::n_dofs_max> p_next_block{};

  // Restrict solution to cell
  const std::uint32_t* u_dofs = u_dofs_[id];
  for (std::size_t i = 0; i < u_dofs_.stride(); ++i){
      u_prev_block[i] = u_prev_data_[u_dofs[i]];
      u_next_block[i] = u_next_data_[u_dofs[i]];
  }
  if (include_pressure){
      const std::uint32_t* p_dofs = p_dofs_[id];
      for (std::size_t i = 0; i < p_dofs_.stride(); ++i){
          p_prev_block[i] = p_prev_data_[p_dofs[i]];
          p_next_block[i] = p_next_data_[p_dofs[i]];
      }
  }

  // Evaluate
  Vector3d U_prev = {std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, u_prev_block.begin(), 0.0),
                     std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, &u_prev_block[ncoeffs_u], 0.0),
                     0.0};
  Vector3d U_next = {std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, u_next_block.begin(), 0.0),
                     std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, &u_next_block[ncoeffs_u], 0.0),
                     0.0 };


  // Update
  fields.U = _alpha_t * U_next + (1-_alpha_t) * U_prev;
  fields.A = stamp_rate(U_next, U_prev, t_prev, t_next);

  if (include_pressure){
    // Evaluate
    double P_prev = std::inner_product(_Np_.data(), _Np_.data()+ncoeffs_p, p_prev_block.begin(), 0.0);
    double P_next = std::inner_product(_Np_.data(), _Np_.data()+ncoeffs_p, p_next_block.begin(), 0.0);
    fields.P = _alpha_t * P_next + (1-_alpha_t) * P_prev;
  }

  if (wants_gradient()){
    if (ncoeffs_u == 3){
      triangles_[id].linearderiv(r, s, u, _Nux_.data(), _Nuy_.data());
    }
    else if (ncoeffs_u == 6){
      triangles_[id].quadderiv(r, s, u, _Nux_.data(), _Nuy_.data());
    }

    Matrix3d gradU_prev;
    gradU_prev <<
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, u_prev_block.begin(), 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, u_prev_block.begin(), 0.0),
      0.0,
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, &u_prev_block[ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, &u_prev_block[ncoeffs_u], 0.0),
      0.0,
      0.0,
      0.0,
      0.0;
    Matrix3d gradU_next;
    gradU_next << 
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, u_next_block.begin(), 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, u_next_block.begin(), 0.0),
      0.0,
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, &u_next_block[ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, &u_next_block[ncoeffs_u], 0.0),
      0.0,
      0.0,
      0.0,
      0.0;

    fields.gradU = _alpha_t * gradU_next + (1-_alpha_t) * gradU_prev;
    fields.gradA = stamp_rate(gradU_next, gradU_prev, t_prev, t_next);
  }
}

#endif
