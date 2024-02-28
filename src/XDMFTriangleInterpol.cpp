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

std::vector<std::pair<double, std::vector<std::string>>> parse_xdmf(const std::string& xdmffilename, 
                                                                    std::string& h5filename,
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
  h5filename = topology_path.substr(0, topology_pos);
  topology_path.erase(0, topology_pos + 1);
  if (h5filename != geometry_path.substr(0, topology_pos)){
    std::cout << "XDMF error: Not matching filenames." << std::endl;
    exit(0);
  }
  geometry_path.erase(0, topology_pos + 1);

  std::vector<std::pair<double, std::vector<std::string>>> titems;
  for (auto & p : tree.get_child("Xdmf.Domain.Grid")) {
    //std :: cout << "[" << p.first << "]" << std :: endl;    
    if (p.first == "Grid"){
      double time;
      std::string filename;
      std::string location;

      for (auto & pp : p.second) {
        if ( pp.first == "Time" ){
          time = pp.second.get<double>("<xmlattr>.Value");
        }
        else if ( pp.first == "Attribute")
        {
          location = pp.second.get<std::string>("DataItem");
          auto pos = location.find(":");
          filename = location.substr(0, pos);
          location.erase(0, pos + 1);
          // Check filename == h5filename too
        }
      }
      //std::cout << " " << time << " " << location << std::endl;
      //titems.push_back({time, {filename, location}});
      std::vector<std::string> path = {filename, location};
      titems.push_back({time, path});
    }
  }
  return titems;
}

std::vector<std::pair<double, std::vector<std::string>>> parse_xdmf(const std::string& xdmffilename)
{
  std::string _dummy0, _dummy1, _dummy2;
  return parse_xdmf(xdmffilename, _dummy0, _dummy1, _dummy2);
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

  // Default params

  dolfin_params["include_pf"] = "false";

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
  include_phi = (dolfin_params["include_phi"] == "true");

  // Make the mesh
  dolfin::Mesh mesh_in;

  // find velocity file
  std::string xdmffilename_u = get_folder() + "/" + dolfin_params["u"];
  std::string xdmffilename_p = get_folder() + "/" + dolfin_params["p"];
  std::string xdmffilename_phi = get_folder() + "/" + dolfin_params["phi"];
  //if (include_phi)
  //std::string xdmffilename_phi = get_folder() + "/" + dolfin_params["phi"];

  std::vector<std::pair<double, std::vector<std::string>>> titems_u, titems_p, titems_phi;

  std::string topology_path, geometry_path;
  titems_u = parse_xdmf(xdmffilename_u, h5filename_u, topology_path, geometry_path);
  std::cout << "mesh: " << h5filename_u << ": " << topology_path << " " << geometry_path << std::endl; 
  
  ts.initialize(titems_u);

  if (include_pressure){
    titems_p = parse_xdmf(xdmffilename_p);
    // std::cout << "h5filename_p = " << h5filename_p << std::endl;
    ts.add("p", titems_p);
  }
  if (include_phi){
    titems_phi = parse_xdmf(xdmffilename_phi);
    // std::cout << "h5filename_phi = " << h5filename_phi << std::endl;
    ts.add("phi", titems_phi);
  }

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
  //ufc_cells_.resize(mesh->num_cells());
  coordinate_dofs_.resize(mesh->num_cells());
  cell2cells_.resize(mesh->num_cells());

  for (std::size_t i = 0; i < mesh->num_cells(); ++i)
  {
    dolfin::Cell dolfin_cell(*mesh, i);
    dolfin_cell.get_coordinate_dofs(coordinate_dofs_[i]);

    triangles_[i] = Triangle(dolfin_cell);
    dolfin_cells_[i] = dolfin_cell;
    //dolfin_cell.get_cell_data(ufc_cells_[i]);
  }
  // Build cell neighbour list for lookup speed
  build_neighbor_list(cell2cells_, mesh, dolfin_cells_);

  // Identify edge cells
  cell_type_.resize(mesh->num_cells());
  //label_cell_type(cell_type_, cell2cells_, dim);

  double tol = 1e-12; // heuristic

  apply_periodic_boundaries(cell2cells_, periodic, x_min, x_max, mesh, dolfin_cells_, dim, tol);

  label_cell_type(cell_type_, cell2cells_, dim);

  //compute_normals(cell_type_, triangles_);
  cell_normal_.resize(mesh->num_cells());
  cell_facet_midpoint_.resize(mesh->num_cells());
  perm_.resize(mesh->num_cells());

  for ( Uint i=1; i < mesh->num_cells(); ++i)
  {
    if (cell_type_[i] == 1)
    {
      auto facets = dolfin_cells_[i].entities(dim-1);
      for ( std::size_t j = 0; j < dolfin_cells_[i].num_entities(dim-1); ++j ){
        dolfin::Facet dolfin_facet(*mesh, facets[j]);

        if (dolfin_facet.exterior()){
          Vector3d pt(dolfin_facet.midpoint().coordinates());

          bool periodic_facet = false;
          for ( Uint k=0; k < dim; ++k)
          {
            if (periodic[k] && (pt[k] < x_min[k] + tol || pt[k] > x_max[k] - tol))
            {
              periodic_facet = true;
              break;
            }
          }
          if (!periodic_facet){
            Vector3d n_loc (dolfin_facet.normal().coordinates());
            cell_normal_[i] = n_loc;
            cell_facet_midpoint_[i] = pt;
          }
        }
      }
    }
  }

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
  phi_prev_data_.resize(xdof.size()/2);
  phi_next_data_.resize(xdof.size()/2);

  //u_prev_coefficients_.resize(dim*ncoeffs_u);
  //u_next_coefficients_.resize(dim*ncoeffs_u);

  //Nu_.resize(ncoeffs_u);
  //Nux_.resize(ncoeffs_u);
  //Nuy_.resize(ncoeffs_u);

  //if (include_pressure){
    //p_prev_coefficients_.resize(ncoeffs_p); // not needed?
    //p_next_coefficients_.resize(ncoeffs_p);
    
    //Np_.resize(ncoeffs_p);
  //}

  std::cout << "Setting max threads: " << omp_get_max_threads() << std::endl;

  found_same_.resize(omp_get_max_threads());
  found_nneigh_.resize(omp_get_max_threads());
  found_other_.resize(omp_get_max_threads());

  // Precomputing dofs
  u_dofs_.resize(mesh->num_cells());
  p_dofs_.resize(mesh->num_cells());

  const dolfin::GenericDofMap& u_dofmap = *u_space_->dofmap();
  const dolfin::GenericDofMap& p_dofmap = *p_space_->dofmap();
  
  for (Uint id = 0; id < mesh->num_cells(); ++id)
  {
    auto u_dofs = u_dofmap.cell_dofs(dolfin_cells_[id].index());
    u_dofs_[id].resize(u_dofs.size());
    for (std::size_t i = 0; i < u_dofs.size(); ++i){
      u_dofs_[id][i] = u_dofs[i];
    }

    auto p_dofs = p_dofmap.cell_dofs(dolfin_cells_[id].index());
    p_dofs_[id].resize(p_dofs.size());
    for (std::size_t i = 0; i < p_dofs.size(); ++i){
      p_dofs_[id][i] = p_dofs[i];
    }
  }
}

void XDMFTriangleInterpol::update(const double t)
{
  MultiStampPair sp = ts.get(t);

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
      if (include_phi)
        phi_prev_data_.swap(phi_next_data_);
    }
    else 
    {
      auto u_path_prev = ts.get_path("u", sp.prev.it);
      std::cout << "Prev: Timestep = " << sp.prev.t << ", file = " << u_path_prev[0] << ":" << u_path_prev[1] << std::endl;
      read_dataset_vector(u_path_prev[0], u_path_prev[1], data_);
      shuffle(u_prev_data_, data_, j2i, 2);

      if (include_pressure){
        auto p_path_prev = ts.get_path("p", sp.prev.it);
        read_dataset_scalar(p_path_prev[0], p_path_prev[1], data_);
        shuffle(p_prev_data_, data_, j2i, 1);
      }

      if (include_phi){
        auto phi_path_prev = ts.get_path("phi", sp.prev.it);
        read_dataset_scalar(phi_path_prev[0], phi_path_prev[1], data_);
        shuffle(phi_prev_data_, data_, j2i, 1);
      }
    }

    auto u_path_next = ts.get_path("u", sp.next.it);
    std::cout << "Next: Timestep = " << sp.next.t << ", file = " << u_path_next[0] << ":" << u_path_next[1] << std::endl;
    read_dataset_vector(u_path_next[0], u_path_next[1], data_);
    shuffle(u_next_data_, data_, j2i, 2);

    if (include_pressure){
      auto p_path_next = ts.get_path("p", sp.next.it);
      read_dataset_scalar(p_path_next[0], p_path_next[1], data_);
      shuffle(p_next_data_, data_, j2i, 1);
    }

    if (include_phi){
      auto phi_path_next = ts.get_path("phi", sp.next.it);
      read_dataset_scalar(phi_path_next[0], phi_path_next[1], data_);
      shuffle(phi_next_data_, data_, j2i, 1);
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

Vector3d XDMFTriangleInterpol::_modx(const Vector3d &x){
  Vector3d x_loc;
  for (std::size_t i=0; i<dim; ++i){
    if (periodic[i]){
      x_loc[i] = x_min[i] + modulox(x[i]-x_min[i], x_max[i]-x_min[i]);
    }
    else {
      x_loc[i] = x[i];
    }
  }
  return x_loc;
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

  auto xx_loc = _modx(x);
  
  bool found = false;
  bool inside_loc = false;
  unsigned int id = 0;

  // Search in neighborhood first
  if (id_prev >= 0){
    if (triangles_[id_prev].contains(xx_loc))
    {
      id = id_prev;
      inside_loc = true;
      found = true;
      found_same_[omp_get_thread_num()]++;
    }
    else {
      for ( auto neigh_id : cell2cells_[id_prev]){
        if (triangles_[neigh_id].contains(xx_loc))
        {
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
    // Index of cell containing point
    dolfin::Array<double> x_loc(dim);
    _modx(x_loc, x);
    const dolfin::Point point(dim, x_loc.data());

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
  double r1, r2, r3;
  triangles_[id].xy2bary(x_loc[0], x_loc[1], r1, r2, r3);

  std::vector<double> _Nu_(ncoeffs_u);
  std::vector<double> _Np_(ncoeffs_p);
  std::vector<double> _Nux_(ncoeffs_u);
  std::vector<double> _Nuy_(ncoeffs_u);

  triangles_[id].linearbasis(r1, r2, r3, _Nu_);

  if (include_pressure){
    triangles_[id].linearbasis(r1, r2, r3, _Np_);
  }

  std::vector<double> u_prev_block(ncoeffs_u*dim);
  std::vector<double> u_next_block(ncoeffs_u*dim);

  // Restrict solution to cell
  //const dolfin::GenericDofMap& u_dofmap = *u_space_->dofmap();
  //auto u_dofs = u_dofmap.cell_dofs(dolfin_cells_[id].index());

  for (std::size_t i = 0; i < u_dofs_[id].size(); ++i){
      u_prev_block[i] = u_prev_data_[u_dofs_[id][i]];
      u_next_block[i] = u_next_data_[u_dofs_[id][i]];
  }

  // Evaluate
  Vector3d U_prev = {std::inner_product(_Nu_.begin(), _Nu_.end(), u_prev_block.begin(), 0.0),
                     std::inner_product(_Nu_.begin(), _Nu_.end(), &u_prev_block[ncoeffs_u], 0.0),
                     0.0};
  Vector3d U_next = {std::inner_product(_Nu_.begin(), _Nu_.end(), u_next_block.begin(), 0.0),
                     std::inner_product(_Nu_.begin(), _Nu_.end(), &u_next_block[ncoeffs_u], 0.0),
                     0.0 };

  Matrix3d gradU_prev, gradU_next;
  if (this->int_order > 1){
    triangles_[id].linearderiv(r1, r2, r3, _Nux_, _Nuy_);

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
  }

  if (cell_type_[id] == 1 && true){
    const Uint ncoeffs_u_2 = 6;
    std::vector<double> _Nu2_(ncoeffs_u_2);

    triangles_[id].quadbasis(r1, r2, r3, _Nu2_);
    std::vector<double> u_prev_block_2(ncoeffs_u_2*dim);
    std::vector<double> u_next_block_2(ncoeffs_u_2*dim);
    
    // Shouldn't do this every time...

    Vector3d uu_prev = {0., 0., 0.};
    Vector3d uu_next = {0., 0., 0.};

    Uint id_uu = 0;
    Uint count_uu = 0;

    for (Uint i=0; i < ncoeffs_u; ++i)
    {
      Vector3d uu_prev_loc = {0., 0., 0.};
      Vector3d uu_next_loc = {0., 0., 0.};
      for (Uint j=0; j < dim; ++j)
      {
        uu_prev_loc[j] = u_prev_block[j*ncoeffs_u + i];
        uu_next_loc[j] = u_next_block[j*ncoeffs_u + i];
      }
      //std::cout << uu << std::endl;
      if (uu_prev_loc.norm() > 1e-14 || uu_next_loc.norm() > 1e-14){
        uu_prev = uu_prev_loc;
        uu_next = uu_next_loc;
        id_uu = i;
        ++count_uu;
      }
    }

    //std::cout << "cell_id=" << id << " id_uu=" << id_uu << " count=" << count_uu << " norm=" << uu.norm() << std::endl;

    // On-the-fly brute force calculation of permutations
    if (perm_[id].size() == 0){
      std::vector<double> u_block = {1.0, 2.0, 3.0}; // dummy data
      double U_1 = std::inner_product(_Nu_.begin(), _Nu_.end(), u_block.begin(), 0.0);

      std::vector<double> u_block_2(ncoeffs_u_2);

      // Corners always work
      for ( Uint i=0; i < ncoeffs_u; ++i)
      {
        u_block_2[i] = u_block[i];
      }

      const std::vector<std::vector<int>> perm_loc = {{3, 4, 5}, {3, 5, 4}, {4, 3, 5}, {4, 5, 3}, {5, 3, 4}, {5, 4, 3}};

      for (int iperm = 0; iperm < perm_loc.size(); ++iperm){
        u_block_2[perm_loc[iperm][0]] = 0.5*(u_block[0] + u_block[1]);
        u_block_2[perm_loc[iperm][1]] = 0.5*(u_block[0] + u_block[2]);
        u_block_2[perm_loc[iperm][2]] = 0.5*(u_block[1] + u_block[2]);

        double U_2 = std::inner_product(_Nu2_.begin(), _Nu2_.end(), u_block_2.begin(), 0.0);
        if (abs(U_2-U_1) < 1e-12){
          perm_[id].insert(perm_[id].end(), perm_loc[iperm].begin(), perm_loc[iperm].end());
          break;
        }
      }
    }
    if (perm_[id].size() == 0){
      std::cout << "ERROR: No permutations worked" << std::endl;
      exit(0);
    }

    if (count_uu == 1 && true)
    {
      const std::vector<std::vector<int>> neigh_dof_loc = {{perm_[id][0], perm_[id][1]}, 
                                                           {perm_[id][0], perm_[id][2]},
                                                           {perm_[id][1], perm_[id][2]}};

      double un_prev = uu_prev.dot(cell_normal_[id]);
      double un_next = uu_next.dot(cell_normal_[id]);

      for ( Uint j=0; j < dim; ++j){
        u_prev_block_2[j*ncoeffs_u_2 + id_uu] = u_prev_block[j*ncoeffs_u + id_uu];
        u_next_block_2[j*ncoeffs_u_2 + id_uu] = u_next_block[j*ncoeffs_u + id_uu];
      }
      for ( auto & i_loc : neigh_dof_loc[id_uu] ){
        for ( Uint j=0; j < dim; ++j){
          u_prev_block_2[j*ncoeffs_u_2 + i_loc] += 0.5*u_prev_block[j*ncoeffs_u + id_uu];
          u_next_block_2[j*ncoeffs_u_2 + i_loc] += 0.5*u_next_block[j*ncoeffs_u + id_uu];

          // Normal correction
          u_prev_block_2[j*ncoeffs_u_2 + i_loc] -= 0.25 * un_prev * cell_normal_[id][j];
          u_next_block_2[j*ncoeffs_u_2 + i_loc] -= 0.25 * un_next * cell_normal_[id][j];
        }
      }
      /*
      double n_dot_grad_g_uu = triangles_[id].dot_grad_gi(cell_normal_[id][0], cell_normal_[id][1], id_uu);
      std::vector<double> tau = {cell_normal_[id][1], -cell_normal_[id][0]};
      std::vector<double> t_dot_grad_g = {triangles_[id].dot_grad_gi(tau[0], tau[1], (id_uu + 1) % (dim+1)),
                                          triangles_[id].dot_grad_gi(tau[0], tau[1], (id_uu + 2) % (dim+1))};

      for (Uint i=0; i<dim+1; ++i){
        for (Uint j=0; j<dim; ++j){
          // Tangential correction
          if (t_dot_grad_g[i] > 1e-12){
            // std::cout << "Larger" << std::endl;
            u_prev_block_2[j*ncoeffs_u_2 + neigh_dof_loc[id_uu][i]] += -0.25 * n_dot_grad_g_uu / t_dot_grad_g[i] * tau[j] * un_prev;
            u_next_block_2[j*ncoeffs_u_2 + neigh_dof_loc[id_uu][i]] += -0.25 * n_dot_grad_g_uu / t_dot_grad_g[i] * tau[j] * un_next;
          }
        }
      }
      */
    }
    else {
      // Linear combination
      for ( Uint j=0; j < dim; ++j){
        for ( Uint i=0; i < ncoeffs_u; ++i)
        {
          u_prev_block_2[j*ncoeffs_u_2 + i] = u_prev_block[j*ncoeffs_u + i];
          u_next_block_2[j*ncoeffs_u_2 + i] = u_next_block[j*ncoeffs_u + i];
        }

        u_prev_block_2[j*ncoeffs_u_2 + perm_[id][0]] = 0.5*(u_prev_block[j*ncoeffs_u + 0] + u_prev_block[j*ncoeffs_u + 1]);
        u_prev_block_2[j*ncoeffs_u_2 + perm_[id][1]] = 0.5*(u_prev_block[j*ncoeffs_u + 0] + u_prev_block[j*ncoeffs_u + 2]);
        u_prev_block_2[j*ncoeffs_u_2 + perm_[id][2]] = 0.5*(u_prev_block[j*ncoeffs_u + 1] + u_prev_block[j*ncoeffs_u + 2]);

        u_next_block_2[j*ncoeffs_u_2 + perm_[id][0]] = 0.5*(u_next_block[j*ncoeffs_u + 0] + u_next_block[j*ncoeffs_u + 1]);
        u_next_block_2[j*ncoeffs_u_2 + perm_[id][1]] = 0.5*(u_next_block[j*ncoeffs_u + 0] + u_next_block[j*ncoeffs_u + 2]);
        u_next_block_2[j*ncoeffs_u_2 + perm_[id][2]] = 0.5*(u_next_block[j*ncoeffs_u + 1] + u_next_block[j*ncoeffs_u + 2]);
      }
    }

    Vector3d U_prev2, U_next2;
    U_prev2 = {std::inner_product(_Nu2_.begin(), _Nu2_.end(), u_prev_block_2.begin(), 0.0),
               std::inner_product(_Nu2_.begin(), _Nu2_.end(), &u_prev_block_2[ncoeffs_u_2], 0.0),
               0.0};
    U_next2 = {std::inner_product(_Nu2_.begin(), _Nu2_.end(), u_next_block_2.begin(), 0.0),
               std::inner_product(_Nu2_.begin(), _Nu2_.end(), &u_next_block_2[ncoeffs_u_2], 0.0),
               0.0 };

    U_prev = U_prev2;
    U_next = U_next2;
    // = {0., 0., 0.}; //
  
    if (this->int_order > 1){
      Matrix3d gradU_prev2, gradU_next2;

      std::vector<double> _Nu2x_(ncoeffs_u_2);
      std::vector<double> _Nu2y_(ncoeffs_u_2);

      triangles_[id].quadderiv(r1, r2, r3, _Nu2x_, _Nu2y_);

      gradU_prev2 <<
        std::inner_product(_Nu2x_.begin(), _Nu2x_.end(), u_prev_block_2.begin(), 0.0),
        std::inner_product(_Nu2y_.begin(), _Nu2y_.end(), u_prev_block_2.begin(), 0.0),
        0.0,
        std::inner_product(_Nu2x_.begin(), _Nu2x_.end(), &u_prev_block_2[ncoeffs_u_2], 0.0),
        std::inner_product(_Nu2y_.begin(), _Nu2y_.end(), &u_prev_block_2[ncoeffs_u_2], 0.0),
        0.0,
        0.0,
        0.0,
        0.0;
      gradU_next2 << 
        std::inner_product(_Nu2x_.begin(), _Nu2x_.end(), u_next_block_2.begin(), 0.0),
        std::inner_product(_Nu2y_.begin(), _Nu2y_.end(), u_next_block_2.begin(), 0.0),
        0.0,
        std::inner_product(_Nu2x_.begin(), _Nu2x_.end(), &u_next_block_2[ncoeffs_u_2], 0.0),
        std::inner_product(_Nu2y_.begin(), _Nu2y_.end(), &u_next_block_2[ncoeffs_u_2], 0.0),
        0.0,
        0.0,
        0.0,
        0.0;

      gradU_prev = gradU_prev2;
      gradU_next = gradU_next2;
    }
  }

  // Update
  fields.U = _alpha_t * U_next + (1-_alpha_t) * U_prev;
  fields.A = (U_next-U_prev)/(t_next-t_prev);

  if (int_order > 1){
    fields.gradU = _alpha_t * gradU_next + (1-_alpha_t) * gradU_prev;
    fields.gradA = (gradU_next-gradU_prev)/(t_next-t_prev);
  }

  if (include_pressure){
    std::vector<double> p_prev_block(ncoeffs_p);
    std::vector<double> p_next_block(ncoeffs_p);
    
    for (std::size_t i = 0; i < p_dofs_[id].size(); ++i){
        p_prev_block[i] = p_prev_data_[p_dofs_[id][i]];
        p_next_block[i] = p_next_data_[p_dofs_[id][i]];
    }
    
    // Evaluate
    double P_prev = std::inner_product(_Np_.begin(), _Np_.end(), p_prev_block.begin(), 0.0);
    double P_next = std::inner_product(_Np_.begin(), _Np_.end(), p_next_block.begin(), 0.0);
    fields.P = _alpha_t * P_next + (1-_alpha_t) * P_prev;
  }

  if (include_phi){
    std::vector<double> phi_prev_block(ncoeffs_p);
    std::vector<double> phi_next_block(ncoeffs_p);
    
    for (std::size_t i = 0; i < p_dofs_[id].size(); ++i){
        phi_prev_block[i] = phi_prev_data_[p_dofs_[id][i]];
        phi_next_block[i] = phi_next_data_[p_dofs_[id][i]];
    }

    double Phi_prev = std::inner_product(_Np_.begin(), _Np_.end(), phi_prev_block.begin(), 0.0);
    double Phi_next = std::inner_product(_Np_.begin(), _Np_.end(), phi_next_block.begin(), 0.0);
    fields.Rho = _alpha_t * Phi_next + (1-_alpha_t) * Phi_prev;
  }

  // cell_type
  fields.cell_type = cell_type_[id];
}

#endif
