#ifdef USE_DOLFIN
#include "Error.hpp"
#include "DolfInterpol.hpp"
#include "loader_params.hpp"
#include "dolfin_ref.hpp"
#include "Params.hpp"
#include "phase_timing.hpp"
#include "H5Cpp.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <map>


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
void DolfInterpol<Cell>::init_mesh_geometry(){
  dim = mesh->geometry().dim();
  mesh->init();
  partrac::phase("mesh init");
  mesh->bounding_box_tree();
  partrac::phase("tree");

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
  hmin_ = mesh->hmin();
  partrac::phase("bounds");
}

template<typename Cell>
void DolfInterpol<Cell>::build_cells(const dolfin::GenericDofMap& dofmap){
  const std::size_t ncells = mesh->num_cells();
  cells_.resize(ncells);
  dolfin_cells_.resize(ncells);

  const std::vector<std::uint32_t> order =
    cell_order(dofmap, ncells, dolfin_params.template get<std::string>("renumber_cells"), dolfin2local_);
  partrac::phase("cell order");
  for (std::size_t l = 0; l < ncells; ++l)
  {
    dolfin::Cell dolfin_cell(*mesh, order[l]);
    cells_[l] = Cell(dolfin_cell);
    dolfin_cells_[l] = dolfin_cell;
  }
  partrac::phase("build cells");
  build_facet_table();
  partrac::phase("facet table");
}

template<typename Cell>
void DolfInterpol<Cell>::build_facet_table()
{
  build_facet_neighbours(facet_neigh_, mesh, dolfin_cells_,
                         dolfin2local_.empty() ? nullptr : &dolfin2local_,
                         periodic, x_min, x_max, dim, periodic_tol);
  set_period();
}

template<typename Cell>
bool DolfInterpol<Cell>::locate_tree(const Vector3d& xx, CellPos& pos)
{
  return tree_to_cell(cells_, *mesh, dim, xx, pos,
                      dolfin2local_.empty() ? nullptr : &dolfin2local_);
}

namespace {
// An element's value components
std::size_t components(const dolfin::FiniteElement& element){
  return element.value_rank() == 0 ? 1 : element.value_dimension(0);
}
}

template<typename Cell>
std::vector<std::array<std::uint32_t, 2>> DolfInterpol<Cell>::seam_pairs(const dolfin::FunctionSpace& space) const
{
  const dolfin::FiniteElement& element = *space.element();
  const dolfin::GenericDofMap& dofmap = *space.dofmap();
  const std::size_t ndofs = element.space_dimension();
  const std::size_t per_component = ndofs / components(element);
  const double tol = 1e-6 * hmin_;
  std::map<std::array<long long, 4>, std::uint32_t> first;
  std::vector<std::array<std::uint32_t, 2>> pairs;
  boost::multi_array<double, 2> X;
  for (std::size_t l = 0; l < dolfin_cells_.size(); ++l){
    const double* v = coordinate_dofs_.data() + l*ncoords_;
    bool on_seam = false;
    for (Uint k = 0; k < ncoords_; ++k){
      const Uint d = k % dim;
      on_seam |= periodic[d] && (std::abs(v[k] - x_min[d]) < tol || std::abs(v[k] - x_max[d]) < tol);
    }
    if (!on_seam) continue;
    element.tabulate_dof_coordinates(X, std::vector<double>(v, v + ncoords_), dolfin_cells_[l]);
    const auto dofs = dofmap.cell_dofs(dolfin_cells_[l].index());
    for (std::size_t i = 0; i < ndofs; ++i){
      // Folded onto the min faces, by component
      std::array<long long, 4> key = {0, 0, 0, (long long)(i / per_component)};
      bool seam = false;
      for (Uint d = 0; d < dim; ++d){
        double x = X[i][d];
        if (periodic[d] && std::abs(x - x_max[d]) < tol){ x -= x_max[d] - x_min[d]; seam = true; }
        else if (periodic[d] && std::abs(x - x_min[d]) < tol) seam = true;
        key[d] = std::llround((x - x_min[d]) / tol);
      }
      if (!seam) continue;
      const auto it = first.emplace(key, std::uint32_t(dofs[i])).first;
      if (it->second != std::uint32_t(dofs[i]))
        pairs.push_back({it->second, std::uint32_t(dofs[i])});
    }
  }
  std::sort(pairs.begin(), pairs.end());
  pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
  return pairs;
}

template<typename Cell>
void DolfInterpol<Cell>::check_seam(const dolfin::Function& f, const std::vector<double>& data,
                                    const std::vector<std::array<std::uint32_t, 2>>& pairs,
                                    const std::string& filename, const char* what) const
{
  if (pairs.empty()) return;
  const dolfin::FiniteElement& element = *f.function_space()->element();
  const dolfin::GenericDofMap& dofmap = *f.function_space()->dofmap();
  const std::size_t ncomp = components(element);
  const std::size_t per_component = element.space_dimension() / ncomp;
  // Each dof's component, and each component's range
  std::vector<std::uint8_t> comp(data.size(), 0);
  std::vector<double> lo(ncomp, std::numeric_limits<double>::infinity());
  std::vector<double> hi(ncomp, -std::numeric_limits<double>::infinity());
  for (std::size_t c = 0; c < mesh->num_cells(); ++c){
    const auto dofs = dofmap.cell_dofs(c);
    for (Eigen::Index i = 0; i < dofs.size(); ++i){
      const std::size_t k = std::size_t(i) / per_component;
      const double x = data[dofs[i]];
      comp[dofs[i]] = std::uint8_t(k);
      lo[k] = std::min(lo[k], x);
      hi[k] = std::max(hi[k], x);
    }
  }
  // Past 1e-6 of the range, and past the values' round-off
  double jump = 0., range = 0., excess = 0.;
  for (const auto& q : pairs){
    const std::size_t k = comp[q[0]];
    const double j = std::abs(data[q[0]] - data[q[1]]);
    const double tol = 1e-6*(hi[k] - lo[k]) + 1e-14*std::max(std::abs(lo[k]), std::abs(hi[k]));
    if (j - tol > excess){
      excess = j - tol;
      jump = j;
      range = hi[k] - lo[k];
    }
  }
  if (excess > 0.){
    partrac::fail("DolfInterpol: in ", filename, " the ", what, " differs between periodic images by ",
                  jump, " (", jump / range, " of its range, more than 1e-6). dolfin's periodic "
                  "P3 space is wrong at the seam on a mesh whose vertex numbering reverses a seam edge "
                  "against its image: number the vertices so that images keep their orientation, or "
                  "write the field from a space without the periodic constraint.");
  }
}

namespace {
// A one-dimensional dataset whole, or false if there is none
template<typename T>
bool read_1d(const H5::H5File& f, const std::string& path, const H5::PredType& type, std::vector<T>& out){
  if (!f.nameExists(path)) return false;
  const H5::DataSet set = f.openDataSet(path);
  const H5::DataSpace space = set.getSpace();
  if (space.getSimpleExtentNdims() != 1) return false;
  hsize_t n = 0;
  space.getSimpleExtentDims(&n);
  out.resize(n);
  set.read(out.data(), type);
  return true;
}
}

template<typename Cell>
void DolfInterpol<Cell>::read_by_cells(const std::string& filename, const std::string& name,
                                       const dolfin::Function& f, std::vector<double>& data) const
{
  std::vector<std::uint64_t> cells, x_cell_dofs;
  std::vector<std::int64_t> cell_dofs;
  std::vector<double> values;
  try {
    H5::Exception::dontPrint();
    const H5::H5File h(get_folder() + "/" + filename, H5F_ACC_RDONLY);
    // As dolfin's read: a group's vector_0, else its vector; or a vector named, in its group
    std::string group = name, vector = name + "/vector_0";
    if (!(h.nameExists(name) && h.childObjType(name) == H5O_TYPE_GROUP)){
      group = name.substr(0, name.rfind('/'));
      vector = name;
    }
    const std::size_t v0 = vector.rfind("/vector_0");
    if (!(read_1d(h, group + "/cells", H5::PredType::NATIVE_UINT64, cells) &&
          read_1d(h, group + "/x_cell_dofs", H5::PredType::NATIVE_UINT64, x_cell_dofs) &&
          read_1d(h, group + "/cell_dofs", H5::PredType::NATIVE_INT64, cell_dofs) &&
          (read_1d(h, vector, H5::PredType::NATIVE_DOUBLE, values) ||
           (v0 != std::string::npos &&
            read_1d(h, vector.substr(0, v0) + "/vector", H5::PredType::NATIVE_DOUBLE, values))))){
      partrac::fail("DolfInterpol: ", filename, " holds no function ", name);
    }
  } catch (const H5::Exception& e) {
    partrac::fail("DolfInterpol: reading ", name, " from ", filename, ": ", e.getDetailMsg());
  }
  const std::size_t nrows = cells.size();
  if (nrows != mesh->num_cells() || x_cell_dofs.size() != nrows + 1 ||
      x_cell_dofs.back() - x_cell_dofs.front() != cell_dofs.size()){
    partrac::fail("DolfInterpol: ", name, " in ", filename, " is not a function on the mesh's ",
                  mesh->num_cells(), " cells");
  }
  std::vector<std::int64_t> local_of(mesh->num_cells(), -1);
  for (dolfin::CellIterator c(*mesh); !c.end(); ++c){
    const std::size_t g = std::size_t(c->global_index());
    if (g < local_of.size()) local_of[g] = c->index();
  }
  const dolfin::GenericDofMap& dofmap = *f.function_space()->dofmap();
  data.assign(f.vector()->local_size(), 0.);
  for (std::size_t k = 0; k < nrows; ++k){
    if (cells[k] >= local_of.size() || local_of[cells[k]] < 0)
      partrac::fail("DolfInterpol: ", name, " names cell ", cells[k], ", not in the mesh");
    const auto dofs = dofmap.cell_dofs(local_of[cells[k]]);
    if (std::size_t(dofs.size()) != x_cell_dofs[k + 1] - x_cell_dofs[k])
      partrac::fail("DolfInterpol: ", name, " has ", x_cell_dofs[k + 1] - x_cell_dofs[k],
                    " dofs a cell, the space ", dofs.size());
    for (Eigen::Index i = 0; i < dofs.size(); ++i){
      const std::int64_t j = cell_dofs[x_cell_dofs[k] - x_cell_dofs.front() + i];
      if (j < 0 || std::size_t(j) >= values.size())
        partrac::fail("DolfInterpol: ", name, " has dof ", j, ", past its ", values.size(), " values");
      data[dofs[i]] = values[j];
    }
  }
}

template<typename Cell>
void DolfInterpol<Cell>::read_field(dolfin::HDF5File& file, const std::string& filename,
                                    const std::string& name, dolfin::Function& f,
                                    std::vector<double>& data, const bool by_cells,
                                    const std::vector<std::array<std::uint32_t, 2>>& seam,
                                    const char* what) const
{
  if (!by_cells){
    file.read(f, name);
    f.vector()->get_local(data);
    return;
  }
  read_by_cells(filename, name, f, data);
  check_seam(f, data, seam, filename, what);
}

template<typename Cell>
DolfInterpol<Cell>::DolfInterpol(const std::string& infilename)
  : MeshCore<Cell>(infilename)
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
  // P3 seam edge dofs pair by edge orientation under the constraint: read unconstrained
  const bool any_periodic = std::any_of(periodic.begin(), periodic.begin() + dim, [](bool b){ return b; });
  const auto domain_of = [&](const std::string& el){
    return any_periodic && el == "P3" ? nullptr : constrained_domain;
  };
  u_space_ = lagrange_space<D, true>(u_el, mesh, domain_of(u_el), "velocity");
  if (include_pressure)
    p_space_ = lagrange_space<D, false>(p_el, mesh, domain_of(p_el), "pressure");
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
  const Uint u_value_size = components(*u_element_);
  Uint p_value_size = 1;
  if (include_pressure){
    p_prev_ = std::make_shared<dolfin::Function>(p_space_);
    p_next_ = std::make_shared<dolfin::Function>(p_space_);
    p_element_ = p_space_->element();
    p_dim_ = p_element_->space_dimension();
    p_value_size = components(*p_element_);
  }
  if (u_value_size != dim || p_value_size != 1){
    partrac::fail("DolfInterpol: velocity has value size ", u_value_size, " and pressure ",
                  p_value_size, ", against ", dim, " and 1");
  }
  u_dofs_.build(*u_space_->dofmap(), dolfin_cells_, "DolfInterpol");
  if (include_pressure)
    p_dofs_.build(*p_space_->dofmap(), dolfin_cells_, "DolfInterpol");
  u_by_cells_ = !domain_of(u_el);
  p_by_cells_ = include_pressure && !domain_of(p_el);
  if (u_by_cells_) u_seam_ = seam_pairs(*u_space_);
  if (p_by_cells_) p_seam_ = seam_pairs(*p_space_);
}

template<typename Cell>
void DolfInterpol<Cell>::update(const double t){
  StampPair sp = ts.get(t);

  // Always load once; keep last bracket past t_max
  if ( !is_initialized || ((t_prev != sp.prev.t || t_next != sp.next.t) && t < ts.get_t_max()) ){
    const std::string u_field = dolfin_params.template get<std::string>("velocity_field");
    const std::string p_field = dolfin_params.template get<std::string>("pressure_field");
    // A stamp's velocity and pressure
    const auto read_stamp = [&](const std::string& filename, const std::shared_ptr<dolfin::Function>& u,
                                std::vector<double>& u_data, const std::shared_ptr<dolfin::Function>& p,
                                std::vector<double>& p_data){
      dolfin::HDF5File file(MPI_COMM_WORLD, get_folder() + "/" + filename, "r");
      read_field(file, filename, u_field, *u, u_data, u_by_cells_, u_seam_, "velocity");
      if (include_pressure)
        read_field(file, filename, p_field, *p, p_data, p_by_cells_, p_seam_, "pressure");
    };

    // Swap if possible
    if (is_initialized && t_next == sp.prev.t){
      std::cout << "Prev: Timestep = " << sp.prev.t << ", swapping... " << std::endl;
      u_prev_data_.swap(u_next_data_);
      p_prev_data_.swap(p_next_data_);
    }
    else {
      std::cout << "Prev: Timestep = " << sp.prev.t << ", filename = " << sp.prev.filename << std::endl;
      read_stamp(sp.prev.filename, u_prev_, u_prev_data_, p_prev_, p_prev_data_);
    }

    std::cout << "Next: Timestep = " << sp.next.t << ", filename = " << sp.next.filename << std::endl;
    // Single stamp: copy prev
    if (sp.next.filename == sp.prev.filename){
      u_next_data_ = u_prev_data_;
      p_next_data_ = p_prev_data_;
    }
    else {
      read_stamp(sp.next.filename, u_next_, u_next_data_, p_next_, p_next_data_);
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
  assert(in_bracket(t, t_prev, t_next, this->stamp_snap));
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
