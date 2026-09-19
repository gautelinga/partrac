#ifdef USE_DOLFIN
#include "Error.hpp"
#include "DolfInterpol.hpp"
#include "loader_params.hpp"
#include "PeriodicBC.hpp"
#include "Params.hpp"
#include "H5Cpp.h"
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
#include <cassert>

namespace {

// A Lagrange space P1-P3 by name, vector or scalar, of dimension D
template<int D, bool Vector>
std::shared_ptr<dolfin::FunctionSpace> lagrange_space(const std::string& el,
                                                      std::shared_ptr<dolfin::Mesh> mesh,
                                                      std::shared_ptr<const dolfin::SubDomain> cd,
                                                      const char* what){
  if constexpr (D == 2 && Vector){
    if (el == "P1") return std::make_shared<vP1_2::FunctionSpace>(mesh, cd);
    if (el == "P2") return std::make_shared<vP2_2::FunctionSpace>(mesh, cd);
    if (el == "P3") return std::make_shared<vP3_2::FunctionSpace>(mesh, cd);
  }
  else if constexpr (D == 2){
    if (el == "P1") return std::make_shared<P1_2::FunctionSpace>(mesh, cd);
    if (el == "P2") return std::make_shared<P2_2::FunctionSpace>(mesh, cd);
    if (el == "P3") return std::make_shared<P3_2::FunctionSpace>(mesh, cd);
  }
  else if constexpr (Vector){
    if (el == "P1") return std::make_shared<vP1_3::FunctionSpace>(mesh, cd);
    if (el == "P2") return std::make_shared<vP2_3::FunctionSpace>(mesh, cd);
    if (el == "P3") return std::make_shared<vP3_3::FunctionSpace>(mesh, cd);
  }
  else {
    if (el == "P1") return std::make_shared<P1_3::FunctionSpace>(mesh, cd);
    if (el == "P2") return std::make_shared<P2_3::FunctionSpace>(mesh, cd);
    if (el == "P3") return std::make_shared<P3_3::FunctionSpace>(mesh, cd);
  }
  partrac::fail("unrecognized ", what, " element: ", el);
}

}  // namespace

Uint dolfin_mesh_dim(const std::string& infilename){
  const std::string mesh_key = partrac::peek_file(infilename, "mesh");
  if (mesh_key.empty()) return 0;
  const std::string meshfilename = infilename.substr(0, infilename.find_last_of("/")) + "/" + mesh_key;
  try {
    H5::Exception::dontPrint();
    H5::H5File f(meshfilename, H5F_ACC_RDONLY);
    H5::DataSpace space = f.openDataSet("mesh/coordinates").getSpace();
    hsize_t dims[2] = {0, 0};
    if (space.getSimpleExtentNdims() != 2) return 0;
    space.getSimpleExtentDims(dims);
    return Uint(dims[1]);
  } catch (const H5::Exception&) {
    return 0;
  }
}

template<typename Cell>
DolfInterpol<Cell>::DolfInterpol(const std::string& infilename)
  : MeshInterpol<Cell>(infilename)
{
  dolfin_params = partrac::parse_file_or_exit(dolfin_h5_schema("fenics"), infilename);

  std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  ts.initialize(get_folder() + "/" + dolfin_params.template get<std::string>("timestamps"));

  read_mesh_params();

  std::string meshfilename = get_folder() + "/" + dolfin_params.template get<std::string>("mesh");
  dolfin::HDF5File meshfile(MPI_COMM_WORLD, meshfilename, "r");

  dolfin::Mesh mesh_in;
  meshfile.read(mesh_in, "mesh", false);

  mesh = std::make_shared<dolfin::Mesh>(mesh_in);
  init_mesh_geometry();

  auto constrained_domain = std::make_shared<PeriodicBC>(periodic, x_min, x_max, dim);

  const std::string u_el = dolfin_params.template get<std::string>("velocity_space");
  const std::string p_el = dolfin_params.template get<std::string>("pressure_space");
  u_space_ = lagrange_space<D, true>(u_el, mesh, constrained_domain, "velocity");
  if (include_pressure)
    p_space_ = lagrange_space<D, false>(p_el, mesh, constrained_domain, "pressure");
  else
    std::cout << "Note: Ignoring pressure." << std::endl;

  // Cells in the velocity's dof order, with periodic neighbours
  build_cells(*u_space_->dofmap());

  // Flat vertex coordinates and orientations per cell, in the cells' order
  const std::size_t ncells = dolfin_cells_.size();
  ncoords_ = (dim + 1) * dim;
  coordinate_dofs_.resize(ncells * ncoords_);
  std::vector<double> coords;
  for (std::size_t l = 0; l < ncells; ++l){
    dolfin_cells_[l].get_coordinate_dofs(coords);
    assert(coords.size() == ncoords_);
    for (Uint k = 0; k < ncoords_; ++k)
      coordinate_dofs_[l*ncoords_ + k] = coords[k];
  }
  const std::vector<int>& orientations = mesh->cell_orientations();
  if (!orientations.empty()){
    cell_orientations_.resize(ncells);
    for (std::size_t l = 0; l < ncells; ++l)
      cell_orientations_[l] = orientations[dolfin_cells_[l].index()];
  }

  u_prev_ = std::make_shared<dolfin::Function>(u_space_);
  u_next_ = std::make_shared<dolfin::Function>(u_space_);
  u_element_ = u_space_->element();
  u_dim_ = u_element_->space_dimension();
  // Basis buffers in evaluate are sized by value size
  const Uint u_value_size = u_element_->value_rank() == 0 ? 1 : u_element_->value_dimension(0);
  Uint p_value_size = 1;
  if (include_pressure){
    p_prev_ = std::make_shared<dolfin::Function>(p_space_);
    p_next_ = std::make_shared<dolfin::Function>(p_space_);
    p_element_ = p_space_->element();
    p_dim_ = p_element_->space_dimension();
    p_value_size = p_element_->value_rank() == 0 ? 1 : p_element_->value_dimension(0);
  }
  if (u_value_size != dim || p_value_size != 1){
    partrac::fail("DolfInterpol: velocity has value size ", u_value_size, " and pressure ",
                  p_value_size, ", against ", dim, " and 1");
  }
  u_dofs_.build(*u_space_->dofmap(), dolfin_cells_, "DolfInterpol");
  if (include_pressure)
    p_dofs_.build(*p_space_->dofmap(), dolfin_cells_, "DolfInterpol");
}

template<typename Cell>
void DolfInterpol<Cell>::update(const double t){
  StampPair sp = ts.get(t);

  // Always load once; keep last bracket past t_max
  if ( !is_initialized || ((t_prev != sp.prev.t || t_next != sp.next.t) && t < ts.get_t_max()) ){
    const std::string u_field = dolfin_params.template get<std::string>("velocity_field");
    const std::string p_field = dolfin_params.template get<std::string>("pressure_field");

    // Swap if possible
    if (is_initialized && t_next == sp.prev.t){
      std::cout << "Prev: Timestep = " << sp.prev.t << ", swapping... " << std::endl;
      u_prev_data_.swap(u_next_data_);
      p_prev_data_.swap(p_next_data_);
    }
    else {
      std::cout << "Prev: Timestep = " << sp.prev.t << ", filename = " << sp.prev.filename << std::endl;
      dolfin::HDF5File prevfile(MPI_COMM_WORLD, get_folder() + "/" + sp.prev.filename, "r");
      prevfile.read(*u_prev_, u_field);
      u_prev_->vector()->get_local(u_prev_data_);
      if (include_pressure){
        prevfile.read(*p_prev_, p_field);
        p_prev_->vector()->get_local(p_prev_data_);
      }
    }

    std::cout << "Next: Timestep = " << sp.next.t << ", filename = " << sp.next.filename << std::endl;
    // Single stamp: copy prev
    if (sp.next.filename == sp.prev.filename){
      u_next_data_ = u_prev_data_;
      p_next_data_ = p_prev_data_;
    }
    else {
      dolfin::HDF5File nextfile(MPI_COMM_WORLD, get_folder() + "/" + sp.next.filename, "r");
      nextfile.read(*u_next_, u_field);
      u_next_->vector()->get_local(u_next_data_);
      if (include_pressure){
        nextfile.read(*p_next_, p_field);
        p_next_->vector()->get_local(p_next_data_);
      }
    }

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  t_update = t;
}

// Basis from dolfin at x; bary unused
template<typename Cell>
void DolfInterpol<Cell>::evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  const int id = pos.id;
  assert(t <= t_next && t >= t_prev);
  const double alpha_t = stamp_weight(t, t_prev, t_next);
  const Vector3d x_loc = _modx(x);

  const int orientation = cell_orientations_.empty() ? -1 : cell_orientations_[id];
  const double* coordinate_dofs = coordinate_dofs_.data() + id*ncoords_;
  const dolfin::FiniteElement& u_element = *u_element_;

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
    p_element_->evaluate_basis(i, &p_basis, _x, coordinate_dofs, orientation);
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

template class DolfInterpol<Triangle>;
template class DolfInterpol<Tet>;

#endif
