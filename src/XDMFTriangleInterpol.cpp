#ifdef USE_DOLFIN
#include "XDMFTriangleInterpol.hpp"
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

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "dolfin/io/HDF5Interface.h"

namespace pt = boost::property_tree;

//using namespace H5;

void read_dataset_scalar(std::string& h5filename, std::string& field, std::vector<double>& data){
    H5::H5File h5file(h5filename, H5F_ACC_RDONLY);
    H5::DataSet dataset = h5file.openDataSet(field.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    // Move out, probably
    int rank = 2; // dataspace.getSimpleExtentNdims();
    std::vector<hsize_t> shape(rank);
    int ndims = dataspace.getSimpleExtentDims( shape.data(), NULL);
    shape[1] = 1;
    std::pair<std::int64_t, std::int64_t> range = dolfin::MPI::local_range(MPI_COMM_WORLD, shape[0]);

    // Hyperslab selection
    std::vector<hsize_t> offset(rank, 0);
    std::vector<hsize_t> count = shape;

    offset[0] = range.first;
    count[0] = range.second - range.first;

    // Allocate data for shape
    dataspace.selectHyperslab( H5S_SELECT_SET, count.data(), offset.data() );
    H5::DataSpace memspace( rank, count.data() );
    // memspace.selectHyperslab( H5S_SELECT_SET, count.data(), offset.data() );

    std::size_t data_size = 1;
    for (std::size_t i = 0; i < count.size(); ++i)
    {
      data_size *= count[i];
    }
    data.resize(data_size);

    dataset.read( data.data(), H5::PredType::NATIVE_DOUBLE, memspace, dataspace );
}

void shuffle(std::vector<double>& data_, const std::vector<double>& xdata, const std::vector<Uint>& j2i, const int dim){
  for ( Uint j=0; j < j2i.size(); ++j ){
    Uint i = j2i[j];
    for ( Uint k=0; k < dim; ++k)
      data_[dim*j+k] = xdata[dim*i+k];
  }
}

void read_dataset_vector(std::string& h5filename_u, std::string& field, std::vector<double>& data){
    H5::H5File h5file_u(h5filename_u, H5F_ACC_RDONLY);
    H5::DataSet dataset = h5file_u.openDataSet(field.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    // Move out, probably
    int rank = 2; // dataspace.getSimpleExtentNdims();
    std::vector<hsize_t> shape(rank);
    int ndims = dataspace.getSimpleExtentDims( shape.data(), NULL);
    shape[1] = 2;
    std::pair<std::int64_t, std::int64_t> range = dolfin::MPI::local_range(MPI_COMM_WORLD, shape[0]);

    // Hyperslab selection
    std::vector<hsize_t> offset(rank, 0);
    std::vector<hsize_t> count = shape;

    offset[0] = range.first;
    count[0] = range.second - range.first;

    // Allocate data for shape
    dataspace.selectHyperslab( H5S_SELECT_SET, count.data(), offset.data() );
    H5::DataSpace memspace( rank, count.data() );
    // memspace.selectHyperslab( H5S_SELECT_SET, count.data(), offset.data() );

    std::size_t data_size = 1;
    for (std::size_t i = 0; i < count.size(); ++i)
    {
      data_size *= count[i];
    }
    data.resize(data_size);

    dataset.read( data.data(), H5::PredType::NATIVE_DOUBLE, memspace, dataspace );
}

std::vector<std::pair<double, std::string>> parse_xdmf(const std::string& xdmffilename, 
                                                       std::string& h5filename_u, 
                                                       std::string& topology_path, 
                                                       std::string& geometry_path)
{
  // Create empty property tree object
  pt::ptree tree;

  // Parse the XML into the property tree.
  pt::read_xml(xdmffilename, tree);

  topology_path = tree.get<std::string>("Xdmf.Domain.Grid.Grid.Topology.DataItem");
  geometry_path = tree.get<std::string>("Xdmf.Domain.Grid.Grid.Geometry.DataItem");

  auto topology_pos = topology_path.find(":");
  h5filename_u = topology_path.substr(0, topology_pos);
  topology_path.erase(0, topology_pos + 1);
  if (h5filename_u != geometry_path.substr(0, topology_pos)){
    std::cout << "XDMF error: Not matching filenames." << std::endl;
    exit(0);
  }
  geometry_path.erase(0, topology_pos + 1);

  std::vector<std::pair<double, std::string>> titems;
  for (auto & p : tree.get_child("Xdmf.Domain.Grid")) {
    //std :: cout << "[" << p.first << "]" << std :: endl;
    if (p.first == "Grid"){
      double time;
      std::string location;
      for (auto & pp : p.second) {
        if ( pp.first == "Time" ){
          time = pp.second.get<double>("<xmlattr>.Value");
        }
        else if ( pp.first == "Attribute")
        {
          location = pp.second.get<std::string>("DataItem");
          auto pos = location.find(":");
          std::string filename = location.substr(0, pos);
          location.erase(0, pos + 1);
          // Check filename == h5filename_u too
        }
      }
      //std::cout << " " << time << " " << location << std::endl;
      titems.push_back({time, location});
    }
  }
  return titems;
}

XDMFTriangleInterpol::XDMFTriangleInterpol(const std::string& infilename)
  : Interpol(infilename)
{

  // Input file (e.g. dolfin_params.dat)
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

  // Set periodicity
  if (dolfin_params["periodic_x"] == "true"){
    periodic[0] = true;
  }
  if (dolfin_params["periodic_y"] == "true"){
    periodic[1] = true;
  }
  include_pressure = !(dolfin_params["ignore_pressure"] == "true");

  // Make the mesh
  dolfin::Mesh mesh_in;

  // find velocity file
  std::string xdmffilename_u = get_folder() + "/" + dolfin_params["u"];
  std::string xdmffilename_p = get_folder() + "/" + dolfin_params["p"];

  std::string topology_path, geometry_path;
  auto titems = parse_xdmf(xdmffilename_u, h5filename_u, topology_path, geometry_path);
  std::cout << "mesh: " << h5filename_u << " -- " << topology_path << " " << geometry_path << std::endl; 
  if (include_pressure){
    std::string _dummy1, _dummy2;
    auto titems2 = parse_xdmf(xdmffilename_p, h5filename_p, _dummy1, _dummy2);
    std::cout << "h5filename_p = " << h5filename_p << std::endl;
  }

  ts.initialize(titems);

  dolfin::HDF5File meshfile(MPI_COMM_WORLD, h5filename_u, "r");

  std::string cell_type_str = "triangle";
  std::unique_ptr<dolfin::CellType> cell_type(dolfin::CellType::create(cell_type_str));

  std::vector<std::int64_t> coords_shape = dolfin::HDF5Interface::get_dataset_shape(meshfile.h5_id(), geometry_path);

  Uint gdim = 2;

  meshfile.read(mesh_in, topology_path, geometry_path, gdim, *cell_type, -1, coords_shape[0], false);

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

  std::cout << "Built neighbour list" << std::endl;
  
  auto constrained_domain = std::make_shared<PeriodicBC>(periodic, x_min, x_max, dim);
  std::cout << "Made periodic domain." << std::endl;

  std::string u_el = dolfin_params["velocity_space"];
  std::string p_el = dolfin_params["pressure_space"];
  
  // Velocity
  u_space_ = std::make_shared<vP1_2::FunctionSpace>(mesh, constrained_domain);
  ncoeffs_u = 3;

  // Pressure
  p_space_ = std::make_shared<P1_2::FunctionSpace>(mesh, constrained_domain);
  ncoeffs_p = 3;

  // const dolfin::GenericDofMap& u_dofmap = *u_space_->dofmap();
  auto xdof = p_space_->tabulate_dof_coordinates();

  std::map<std::tuple<double, double>, Uint> xmap;
  for ( Uint j=0; j < xdof.size()/2; ++j ){
    //std::cout << " " << x_dof[i] << " " << x_dof[i+1] << std::endl;
    xmap[{xdof[2*j], xdof[2*j+1]}] = j;
  }

  std::vector<double> xdata;
  read_dataset_vector(h5filename_u, geometry_path, xdata);

  j2i.resize(xdof.size()/2);
  for ( Uint i=0; i < xdata.size()/2; ++i ){
    Uint j = xmap[{xdata[2*i], xdata[2*i+1]}];
    j2i[j] = i;
  }

  u_prev_data_.resize(xdof.size());
  u_next_data_.resize(xdof.size());
  p_prev_data_.resize(xdof.size()/2);
  p_next_data_.resize(xdof.size()/2);

  u_prev_coefficients_.resize(dim*ncoeffs_u);
  u_next_coefficients_.resize(dim*ncoeffs_u);

  Nu_.resize(ncoeffs_u);
  Nux_.resize(ncoeffs_u);
  Nuy_.resize(ncoeffs_u);

  if (include_pressure){
    p_prev_coefficients_.resize(ncoeffs_p);
    p_next_coefficients_.resize(ncoeffs_p);
    
    Np_.resize(ncoeffs_p);
  }

  std::cout << "Setting max threads: " << omp_get_max_threads() << std::endl;

  found_same_.resize(omp_get_max_threads());
  found_nneigh_.resize(omp_get_max_threads());
  found_other_.resize(omp_get_max_threads());
}

void XDMFTriangleInterpol::update(const double t)
{
  // TODO: swapping!

  StampPair sp = ts.get(t);
  // std::cout << sp.prev.filename << " " << sp.next.filename << std::endl;

  if (( !is_initialized || t_prev != sp.prev.t || t_next != sp.next.t) && t < ts.get_t_max() )
  {
    std::vector<double> data_;
    // Swap if possible
    if (is_initialized && t_next == sp.prev.t)
    {
      std::cout << "Prev: Timestep = " << sp.prev.t << ", swapping... "<< std::endl;
      u_prev_data_.swap(u_next_data_);
      if (include_pressure)
        p_prev_data_.swap(p_next_data_);
    }
    else 
    {
      std::cout << "Prev: Timestep = " << sp.prev.t << ", file = " << h5filename_u << ":" << sp.prev.filename << std::endl;
      read_dataset_vector(h5filename_u, sp.prev.filename, data_);
      shuffle(u_prev_data_, data_, j2i, 2);

      if (include_pressure){
        read_dataset_scalar(h5filename_p, sp.prev.filename, data_);
        shuffle(p_prev_data_, data_, j2i, 1);
      }
    }

    std::cout << "Next: Timestep = " << sp.next.t << ", file = " << h5filename_u << ":" << sp.next.filename << std::endl;
    read_dataset_vector(h5filename_u, sp.next.filename, data_);
    shuffle(u_next_data_, data_, j2i, 2);

    if (include_pressure){
      read_dataset_scalar(h5filename_p, sp.prev.filename, data_);
      shuffle(p_next_data_, data_, j2i, 1);
    }

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  t_update = t;
}

void XDMFTriangleInterpol::probe(const Vector3d &x, const double t)
{
  int id_prev = -1;
  probe(x, t, id_prev);
}

void XDMFTriangleInterpol::_modx(dolfin::Array<double>& x_loc, const Vector3d &x){
  for (std::size_t i=0; i<dim; ++i){
    if (periodic[i]){
      x_loc[i] = x_min[i] + modulox(x[i]-x_min[i], x_max[i]-x_min[i]);
    }
    else {
      x_loc[i] = x[i];
    }
  }
}

void XDMFTriangleInterpol::probe(const Vector3d &x, const double t, int& id_prev)
{
  // Not good for parallelization
  inside = probe_light(x, t, id_prev);

  if (inside)
  {
    PointValues fields(U0);
    probe_heavy(x, t, id_prev, fields);

    // Update
    U = fields.U;
    A = fields.A;

    if (include_pressure){
      // Evaluate
      P = fields.P;
    }

    if (this->int_order > 1){
      gradU = fields.gradU;
      gradA = fields.gradA;
    }
  }
}

bool XDMFTriangleInterpol::probe_light(const Vector3d &x, const double t, int& id_prev)
{
  // TODO: CHECK if thread safe
  assert(t <= t_next && t >= t_prev);
  
  dolfin::Array<double> x_loc(dim);
  _modx(x_loc, x);

  // Index of cell containing point
  const dolfin::Point point(dim, x_loc.data());
  
  bool found = false;
  bool inside_loc = false;
  unsigned int id = 0;

  // Search in neighborhood first
  if (id_prev >= 0){
    dolfin::Cell prev_cell(*mesh, id_prev);
    if (prev_cell.contains(point)){
      id = id_prev;
      inside_loc = true;
      found = true;
      found_same_[omp_get_thread_num()]++;
    }
    else {
      for ( auto neigh_id : cell2cells_[id_prev]){
        dolfin::Cell neigh_cell(*mesh, neigh_id);
        if (neigh_cell.contains(point)){
          inside_loc = true;
          found = true;
          id = neigh_id;
          found_nneigh_[omp_get_thread_num()]++;
          break;
        }
      }
    }
  }
  if (!found){
    id = mesh->bounding_box_tree()->compute_first_entity_collision(point);
    inside_loc = (id != std::numeric_limits<unsigned int>::max());
    if (inside_loc) {
      found = true;
      found_other_[omp_get_thread_num()]++;
    }
  }
  if (found){
    id_prev = id;
  }
  return inside_loc;
}

void XDMFTriangleInterpol::probe_heavy(const Vector3d &x, const double tin, const int id, PointValues& fields)
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

  triangles_[id].linearbasis(r, s, u, _Nu_);
  if (include_pressure){
    triangles_[id].linearbasis(r, s, u, _Np_);
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
    triangles_[id].linearderiv(r, s, u, _Nux_, _Nuy_);

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
