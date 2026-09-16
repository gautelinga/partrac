#ifdef USE_DOLFIN
#include "DolfInterpol.hpp"
#include "dolfin_helpers.hpp"
#include <omp.h>
#include <boost/algorithm/string.hpp>

DolfInterpol::DolfInterpol(const std::string& infilename) : Interpol(infilename) {
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

  auto constrained_domain = std::make_shared<PeriodicBC>(periodic, x_min, x_max, dim);

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
      exit(1);
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
      exit(1);
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
      exit(1);
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
      exit(1);
    }
    break;
  default:
    std::cout << "Unsupported dimensionality" << std::endl;
    exit(1);
  }

  // Per-cell data for locate
  const std::size_t ncells = mesh->num_cells();
  dolfin_cells_.resize(ncells);
  cell2cells_.resize(ncells);
  if (dim == 2) triangles_.resize(ncells); else tets_.resize(ncells);
  // Flat vertex coordinates per cell
  ncoords_ = (dim + 1) * dim;
  coordinate_dofs_.resize(ncells * ncoords_);
  std::vector<double> coords;
  const std::vector<std::uint32_t> order =
    cell_order(*u_space->dofmap(), ncells, dolfin_params["renumber_cells"], dolfin2local_);
  for (std::size_t l = 0; l < ncells; ++l){
    dolfin::Cell dolfin_cell(*mesh, order[l]);
    dolfin_cell.get_coordinate_dofs(coords);
    assert(coords.size() == ncoords_);
    for (Uint k = 0; k < ncoords_; ++k)
      coordinate_dofs_[l*ncoords_ + k] = coords[k];
    dolfin_cells_[l] = dolfin_cell;
    if (dim == 2) triangles_[l] = Triangle(dolfin_cell); else tets_[l] = Tet(dolfin_cell);
  }
  // Cell orientations (all evaluate needs of ufc::cell), in the cells' order
  cell_orientations_ = mesh->cell_orientations();
  if (!dolfin2local_.empty() && !cell_orientations_.empty()){
    std::vector<int> reordered(ncells);
    for (std::size_t l = 0; l < ncells; ++l) reordered[l] = cell_orientations_[order[l]];
    cell_orientations_.swap(reordered);
  }
  build_neighbor_list(cell2cells_, mesh, dolfin_cells_,
                      dolfin2local_.empty() ? nullptr : &dolfin2local_);
  found_.resize(omp_get_max_threads());
  //std::cout << x_min << std::endl;
  //std::cout << x_max << std::endl;
  //this->Lx = x_max[0]-x_min[0];
  //this->Ly = x_max[1]-x_min[1];
  //this->Lz = x_max[2]-x_min[2];

  u_prev_ = std::make_shared<dolfin::Function>(u_space);
  u_next_ = std::make_shared<dolfin::Function>(u_space);
  p_prev_ = std::make_shared<dolfin::Function>(p_space);
  p_next_ = std::make_shared<dolfin::Function>(p_space);
  u_element_ = u_space->element();
  p_element_ = p_space->element();
  u_dim_ = u_element_->space_dimension();
  p_dim_ = p_element_->space_dimension();
  // Basis buffers in evaluate are sized by value size
  const Uint u_value_size = u_element_->value_rank() == 0 ? 1 : u_element_->value_dimension(0);
  const Uint p_value_size = p_element_->value_rank() == 0 ? 1 : p_element_->value_dimension(0);
  if (u_value_size != dim || p_value_size != 1){
    std::cout << "DolfInterpol: velocity has value size " << u_value_size
              << " and pressure " << p_value_size
              << ", against " << dim << " and 1" << std::endl;
    exit(1);
  }
  u_dofs_.build(*u_space->dofmap(), dolfin_cells_, "DolfInterpol");
  p_dofs_.build(*p_space->dofmap(), dolfin_cells_, "DolfInterpol");

  //std::cout << "GOT THIS FAR" << std::endl;
}

void DolfInterpol::update(const double t){
  StampPair sp = ts.get(t);
  // std::cout << sp.prev.filename << " " << sp.next.filename << std::endl;

  if (!is_initialized || t_prev != sp.prev.t || t_next != sp.next.t){
    // Swap if possible
    if (is_initialized && t_next == sp.prev.t){
      std::cout << "Prev: Timestep = " << sp.prev.t << ", swapping... " << std::endl;
      u_prev_data_.swap(u_next_data_);
      p_prev_data_.swap(p_next_data_);
    }
    else {
      std::cout << "Prev: Timestep = " << sp.prev.t << ", filename = " << sp.prev.filename << std::endl;
      dolfin::HDF5File prevfile(MPI_COMM_WORLD, get_folder() + "/" + sp.prev.filename, "r");
      prevfile.read(*u_prev_, "u");
      prevfile.read(*p_prev_, "p");
      u_prev_->vector()->get_local(u_prev_data_);
      p_prev_->vector()->get_local(p_prev_data_);
    }

    std::cout << "Next: Timestep = " << sp.next.t << ", filename = " << sp.next.filename << std::endl;
    // Single stamp: copy prev
    if (sp.next.filename == sp.prev.filename){
      u_next_data_ = u_prev_data_;
      p_next_data_ = p_prev_data_;
    }
    else {
      dolfin::HDF5File nextfile(MPI_COMM_WORLD, get_folder() + "/" + sp.next.filename, "r");
      nextfile.read(*u_next_, "u");
      nextfile.read(*p_next_, "p");
      u_next_->vector()->get_local(u_next_data_);
      p_next_->vector()->get_local(p_next_data_);
    }

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  t_update = t;
}



Vector3d DolfInterpol::_modx(const Vector3d &x){
  Vector3d x_loc = x;
  for (std::size_t i=0; i<dim; ++i)
    if (periodic[i])
      x_loc[i] = x_min[i] + modulox(x[i]-x_min[i], x_max[i]-x_min[i]);
  return x_loc;
}

bool DolfInterpol::locate(const Vector3d &x, const double t, CellPos& pos){
  assert(t <= t_next && t >= t_prev);
  const Vector3d xx = _modx(x);
  if (dim == 2)
    return locate_in_cells(triangles_, cell2cells_, *mesh, dim, xx, pos,
                           found_, dolfin2local_.empty() ? nullptr : &dolfin2local_);
  return locate_in_cells(tets_, cell2cells_, *mesh, dim, xx, pos,
                         found_, dolfin2local_.empty() ? nullptr : &dolfin2local_);
}
// Basis from dolfin at x; bary unused
void DolfInterpol::evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  const int id = pos.id;
  assert(t <= t_next && t >= t_prev);
  const double alpha_t = stamp_weight(t, t_prev, t_next);
  const Vector3d x_loc = _modx(x);

  const int orientation = cell_orientations_.empty() ? -1 : cell_orientations_[id];
  const double* coordinate_dofs = coordinate_dofs_.data() + id*ncoords_;
  const dolfin::FiniteElement& u_element = *u_element_;
  const dolfin::FiniteElement& p_element = *p_element_;

  const std::uint32_t* u_dofs = u_dofs_[id];
  const std::uint32_t* p_dofs = p_dofs_[id];

  // At most three components
  double u_basis[3], gradu_basis[9], p_basis;

  Vector3d U_prev = {0., 0., 0.};
  Vector3d U_next = {0., 0., 0.};
  double P_prev = 0.;
  double P_next = 0.;
  Matrix3d gradU_prev = Matrix3d::Zero();
  Matrix3d gradU_next = Matrix3d::Zero();

  const double* _x = x_loc.data();

  for (Uint i=0; i<u_dim_; ++i){
    u_element.evaluate_basis(i, u_basis, _x, coordinate_dofs, orientation);
    for (Uint j=0; j<dim; ++j){
      U_prev[j] += u_prev_data_[u_dofs[i]]*u_basis[j];
      U_next[j] += u_next_data_[u_dofs[i]]*u_basis[j];
    }
  }
  for (Uint i=0; i<p_dim_; ++i){
    p_element.evaluate_basis(i, &p_basis, _x, coordinate_dofs, orientation);
    P_prev += p_prev_data_[p_dofs[i]]*p_basis;
    P_next += p_next_data_[p_dofs[i]]*p_basis;
  }

  if (wants_gradient()){
    for (Uint i=0; i<u_dim_; ++i){
      u_element.evaluate_basis_derivatives(i, 1, gradu_basis, _x, coordinate_dofs, orientation);
      for (Uint j=0; j<dim; ++j){
        for (Uint k=0; k<dim; ++k){
          gradU_prev(j, k) += u_prev_data_[u_dofs[i]]*gradu_basis[dim*j+k];
          gradU_next(j, k) += u_next_data_[u_dofs[i]]*gradu_basis[dim*j+k];
        }
      }
    }
  }

  fields.U = alpha_t * U_next + (1-alpha_t) * U_prev;
  fields.A = stamp_rate(U_next, U_prev, t_prev, t_next);
  fields.P = alpha_t * P_next + (1-alpha_t) * P_prev;
  if (wants_gradient()){
    fields.gradU = alpha_t * gradU_next + (1-alpha_t) * gradU_prev;
    fields.gradA = stamp_rate(gradU_next, gradU_prev, t_prev, t_next);
  }
}

#endif
