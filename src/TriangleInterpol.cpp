#ifdef USE_DOLFIN
#include "TriangleInterpol.hpp"
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
  ufc_cells_.resize(mesh->num_cells());
  coordinate_dofs_.resize(mesh->num_cells());
  cell2cells_.resize(mesh->num_cells());

  for (std::size_t i = 0; i < mesh->num_cells(); ++i)
  {
    dolfin::Cell dolfin_cell(*mesh, i);
    dolfin_cell.get_coordinate_dofs(coordinate_dofs_[i]);

    triangles_[i] = Triangle(dolfin_cell);
    dolfin_cells_[i] = dolfin_cell;
    dolfin_cell.get_cell_data(ufc_cells_[i]);
  }
  // Build cell neighbour list for lookup speed
  build_neighbor_list(cell2cells_, mesh, dolfin_cells_);

  /*
  for ( std::size_t i = 0; i < cell2cells_.size(); ++i ){
    std::cout << i << ": ";
    for (auto & cell : cell2cells_[i] ){
       std::cout << " " << cell;
    }
    std::cout << std::endl;
  }
  */
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
    exit(0);
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
      exit(0);
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
}

void TriangleInterpol::update(const double t)
{
  // TODO: swapping!

  StampPair sp = ts.get(t);
  // std::cout << sp.prev.filename << " " << sp.next.filename << std::endl;

  if (!is_initialized || t_prev != sp.prev.t || t_next != sp.next.t){
    std::cout << "Prev: Timestep = " << sp.prev.t << ", filename = " << sp.prev.filename << std::endl;
    dolfin::HDF5File prevfile(MPI_COMM_WORLD, get_folder() + "/" + sp.prev.filename, "r");
    prevfile.read(*u_prev_, "u");
    u_prev_->vector()->get_local(u_prev_data_);
    if (include_pressure){
      prevfile.read(*p_prev_, "p");
      p_prev_->vector()->get_local(p_prev_data_);
    }

    std::cout << "Next: Timestep = " << sp.next.t << ", filename = " << sp.next.filename << std::endl;
    dolfin::HDF5File nextfile(MPI_COMM_WORLD, get_folder() + "/" + sp.next.filename, "r");
    nextfile.read(*u_next_, "u");
    u_next_->vector()->get_local(u_next_data_);
    if (include_pressure){
      nextfile.read(*p_next_, "p");
      p_next_->vector()->get_local(p_next_data_);
    }

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  //alpha_t = sp.weight_next(t);
  t_update = t;
}

void TriangleInterpol::probe(const Vector3d &x, const double t)
{
  int id_prev = -1;
  probe(x, t, id_prev);
}

void TriangleInterpol::_modx(dolfin::Array<double>& x_loc, const Vector3d &x){
  for (std::size_t i=0; i<dim; ++i){
    if (periodic[i]){
      x_loc[i] = x_min[i] + modulox(x[i]-x_min[i], x_max[i]-x_min[i]);
    }
    else {
      x_loc[i] = x[i];
    }
  }
}

void TriangleInterpol::probe(const Vector3d &x, const double t, int& id_prev)
{
  assert(t <= t_next && t >= t_prev);
  alpha_t = (t-t_prev)/(t_next-t_prev);

  // std::cout << "t=" << t << " t_next=" << t_next << " t_prev=" << t_prev << " alpha_t=" << alpha_t << std::endl;

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
  
  bool found = false;

  unsigned int id = 0;
  // Search in neighborhood first
  if (id_prev >= 0){
    dolfin::Cell prev_cell(*mesh, id_prev);
    if (prev_cell.contains(point)){
      id = id_prev;
      inside = true;
      found = true;
      ++found_same;
    }
    else {
      for ( auto neigh_id : cell2cells_[id_prev]){
        dolfin::Cell neigh_cell(*mesh, neigh_id);
        if (neigh_cell.contains(point)){
          inside = true;
          found = true;
          id = neigh_id;
          ++found_nneigh;
          break;
        }
      }
    }
  }
  if (!found){
    id = mesh->bounding_box_tree()->compute_first_entity_collision(point);
    inside = (id != std::numeric_limits<unsigned int>::max());
    if (inside) {
      found = true;
      ++found_other;
    }
  }
  if (found){
    id_prev = id;
  }
  
  if (inside)
  {
    // Compute Pk-Pl basis at x
    double r, s, t;
    triangles_[id].xy2bary(x_loc[0], x_loc[1], r, s, t);

    if (ncoeffs_u == 3){
      triangles_[id].linearbasis(r, s, t, Nu_);
    }
    else if (ncoeffs_u == 6){
      triangles_[id].quadbasis(r, s, t, Nu_);
    }
    else {
      std::cout << "Unrecognized ncoeffs_u = " << ncoeffs_u << std::endl;
      exit(0);
    }
    if (include_pressure){
      if (ncoeffs_p == 3){
        triangles_[id].linearbasis(r, s, t, Np_);
      }
      else if (ncoeffs_p == 6){
        triangles_[id].quadbasis(r, s, t, Np_);
      }
      else {
        std::cout << "Unrecognized ncoeffs_p = " << ncoeffs_p << std::endl;
        exit(0);
      }
    }

    // Restrict solution to cell
    u_prev_->restrict(u_prev_coefficients_.data(), *u_space_->element(), dolfin_cells_[id],
                      coordinate_dofs_[id].data(), ufc_cells_[id]);
    u_next_->restrict(u_next_coefficients_.data(), *u_space_->element(), dolfin_cells_[id],
                      coordinate_dofs_[id].data(), ufc_cells_[id]);
    if (include_pressure){
      p_prev_->restrict(p_prev_coefficients_.data(), *p_space_->element(), dolfin_cells_[id],
                        coordinate_dofs_[id].data(), ufc_cells_[id]);
      p_next_->restrict(p_next_coefficients_.data(), *p_space_->element(), dolfin_cells_[id],
                        coordinate_dofs_[id].data(), ufc_cells_[id]);
    }

    // Evaluate

    Vector3d U_prev = {std::inner_product(Nu_.begin(), Nu_.end(), u_prev_coefficients_.begin(), 0.0),
                       std::inner_product(Nu_.begin(), Nu_.end(), &u_prev_coefficients_[ncoeffs_u], 0.0),
                       0.0};
    Vector3d U_next = {std::inner_product(Nu_.begin(), Nu_.end(), u_next_coefficients_.begin(), 0.0),
                       std::inner_product(Nu_.begin(), Nu_.end(), &u_next_coefficients_[ncoeffs_u], 0.0),
                       0.0 };

    // Update
    U = alpha_t * U_next + (1-alpha_t) * U_prev;
    A = (U_next-U_prev)/(t_next-t_prev);

    if (include_pressure){
      // Evaluate
      double P_prev = std::inner_product(Np_.begin(), Np_.end(), p_prev_coefficients_.begin(), 0.0);
      double P_next = std::inner_product(Np_.begin(), Np_.end(), p_next_coefficients_.begin(), 0.0);
      P = alpha_t * P_next + (1-alpha_t) * P_prev;
    }

    if (this->int_order > 1){
      if (ncoeffs_u == 3){
        triangles_[id].linearderiv(r, s, t, Nux_, Nuy_);
      }
      else if (ncoeffs_u == 6){
        triangles_[id].quadderiv(r, s, t, Nux_, Nuy_);
      }

      Matrix3d gradU_prev;
      gradU_prev <<
        std::inner_product(Nux_.begin(), Nux_.end(), u_prev_coefficients_.begin(), 0.0),
        std::inner_product(Nuy_.begin(), Nuy_.end(), u_prev_coefficients_.begin(), 0.0),
        0.0,
        std::inner_product(Nux_.begin(), Nux_.end(), &u_prev_coefficients_[ncoeffs_u], 0.0),
        std::inner_product(Nuy_.begin(), Nuy_.end(), &u_prev_coefficients_[ncoeffs_u], 0.0),
        0.0,
        0.0,
        0.0,
        0.0;
      Matrix3d gradU_next;
      gradU_next << 
        std::inner_product(Nux_.begin(), Nux_.end(), u_next_coefficients_.begin(), 0.0),
        std::inner_product(Nuy_.begin(), Nuy_.end(), u_next_coefficients_.begin(), 0.0),
        0.0,
        std::inner_product(Nux_.begin(), Nux_.end(), &u_next_coefficients_[ncoeffs_u], 0.0),
        std::inner_product(Nuy_.begin(), Nuy_.end(), &u_next_coefficients_[ncoeffs_u], 0.0),
        0.0,
        0.0,
        0.0,
        0.0;

      gradU = alpha_t * gradU_next + (1-alpha_t) * gradU_prev;
      gradA = (gradU_next-gradU_prev)/(t_next-t_prev);
    }
  }
}

bool TriangleInterpol::probe_light(const Vector3d &x, const double t, int& id_prev)
{
  assert(t <= t_next && t >= t_prev);
  
  // std::cout << "t=" << t << " t_next=" << t_next << " t_prev=" << t_prev << " alpha_t=" << alpha_t << std::endl;

  dolfin::Array<double> x_loc(dim);
  _modx(x_loc, x);

  // Index of cell containing point
  const dolfin::Point point(dim, x_loc.data());
  
  bool found = false;

  unsigned int id = 0;
  // Search in neighborhood first
  if (id_prev >= 0){
    dolfin::Cell prev_cell(*mesh, id_prev);
    if (prev_cell.contains(point)){
      id = id_prev;
      inside = true;
      found = true;
      ++found_same;
    }
    else {
      for ( auto neigh_id : cell2cells_[id_prev]){
        dolfin::Cell neigh_cell(*mesh, neigh_id);
        if (neigh_cell.contains(point)){
          inside = true;
          found = true;
          id = neigh_id;
          ++found_nneigh;
          break;
        }
      }
    }
  }
  if (!found){
    id = mesh->bounding_box_tree()->compute_first_entity_collision(point);
    inside = (id != std::numeric_limits<unsigned int>::max());
    if (inside) {
      found = true;
      ++found_other;
    }
  }
  if (found){
    id_prev = id;
  }
  return inside;
}

void TriangleInterpol::probe_heavy(const Vector3d &x, const double tin, const int id, PointValues& fields)
{
  dolfin::Array<double> x_loc(dim);
  _modx(x_loc, x);

  double _alpha_t = (tin-t_prev)/(t_next-t_prev);

  // Compute Pk-Pl basis at x
  double r, s, u;
  triangles_[id].xy2bary(x_loc[0], x_loc[1], r, s, u);

  std::vector<double> _Nu_(ncoeffs_u);
  std::vector<double> _Np_(ncoeffs_p);
  std::vector<double> _Nux_(ncoeffs_u);
  std::vector<double> _Nuy_(ncoeffs_u);

  if (ncoeffs_u == 3){
    triangles_[id].linearbasis(r, s, u, _Nu_);
  }
  else if (ncoeffs_u == 6){
    triangles_[id].quadbasis(r, s, u, _Nu_);
  }
  else {
    std::cout << "Unrecognized ncoeffs_u = " << ncoeffs_u << std::endl;
    exit(0);
  }
  if (include_pressure){
    if (ncoeffs_p == 3){
      triangles_[id].linearbasis(r, s, u, _Np_);
    }
    else if (ncoeffs_p == 6){
      triangles_[id].quadbasis(r, s, u, _Np_);
    }
    else {
      std::cout << "Unrecognized ncoeffs_p = " << ncoeffs_p << std::endl;
      exit(0);
    }
  }

  std::vector<double> u_prev_block(ncoeffs_u*dim);
  std::vector<double> u_next_block(ncoeffs_u*dim);
  std::vector<double> p_prev_block(ncoeffs_p);
  std::vector<double> p_next_block(ncoeffs_p);

  // Restrict solution to cell
  const dolfin::GenericDofMap& u_dofmap = *u_space_->dofmap();
  auto u_dofs = u_dofmap.cell_dofs(dolfin_cells_[id].index());

  for (std::size_t i = 0; i < u_dofs.size(); ++i){
      u_prev_block[i] = u_prev_data_[u_dofs[i]];
      u_next_block[i] = u_next_data_[u_dofs[i]];
  }
  if (include_pressure){
      const dolfin::GenericDofMap& p_dofmap = *p_space_->dofmap();
      auto p_dofs = p_dofmap.cell_dofs(dolfin_cells_[id].index());

      for (std::size_t i = 0; i < p_dofs.size(); ++i){
          p_prev_block[i] = p_prev_data_[p_dofs[i]];
          p_next_block[i] = p_next_data_[p_dofs[i]];
      }
  }

  // Evaluate
  Vector3d U_prev = {std::inner_product(_Nu_.begin(), _Nu_.end(), u_prev_block.begin(), 0.0),
                     std::inner_product(_Nu_.begin(), _Nu_.end(), &u_prev_block[ncoeffs_u], 0.0),
                     0.0};
  Vector3d U_next = {std::inner_product(_Nu_.begin(), _Nu_.end(), u_next_block.begin(), 0.0),
                     std::inner_product(_Nu_.begin(), _Nu_.end(), &u_next_block[ncoeffs_u], 0.0),
                     0.0 };


  // Update
  fields.U = _alpha_t * U_next + (1-_alpha_t) * U_prev;
  fields.A = (U_next-U_prev)/(t_next-t_prev);

  if (include_pressure){
    // Evaluate
    double P_prev = std::inner_product(_Np_.begin(), _Np_.end(), p_prev_block.begin(), 0.0);
    double P_next = std::inner_product(_Np_.begin(), _Np_.end(), p_next_block.begin(), 0.0);
    fields.P = _alpha_t * P_next + (1-_alpha_t) * P_prev;
  }

  if (this->int_order > 1){
    if (ncoeffs_u == 3){
      triangles_[id].linearderiv(r, s, u, _Nux_, _Nuy_);
    }
    else if (ncoeffs_u == 6){
      triangles_[id].quadderiv(r, s, u, _Nux_, _Nuy_);
    }

    Matrix3d gradU_prev;
    gradU_prev <<
      std::inner_product(_Nux_.begin(), _Nux_.end(), u_prev_block.begin(), 0.0),
      std::inner_product(_Nuy_.begin(), _Nuy_.end(), u_prev_block.begin(), 0.0),
      0.0,
      std::inner_product(_Nux_.begin(), _Nux_.end(), &u_prev_block[ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.begin(), _Nuy_.end(), &u_prev_block[ncoeffs_u], 0.0),
      0.0,
      0.0,
      0.0,
      0.0;
    Matrix3d gradU_next;
    gradU_next << 
      std::inner_product(_Nux_.begin(), _Nux_.end(), u_next_block.begin(), 0.0),
      std::inner_product(_Nuy_.begin(), _Nuy_.end(), u_next_block.begin(), 0.0),
      0.0,
      std::inner_product(_Nux_.begin(), _Nux_.end(), &u_next_block[ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.begin(), _Nuy_.end(), &u_next_block[ncoeffs_u], 0.0),
      0.0,
      0.0,
      0.0,
      0.0;

    fields.gradU = _alpha_t * gradU_next + (1-_alpha_t) * gradU_prev;
    fields.gradA = (gradU_next-gradU_prev)/(t_next-t_prev);
  }
}

#endif
