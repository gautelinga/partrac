#ifdef USE_DOLFIN
#include "XDMFTetInterpol.hpp"
#include <array>
#include <numeric>
#include "Timestamps.hpp"
//#include "H5Cpp.h"
#include <boost/algorithm/string.hpp>
#include <cassert>

#include "dolfin_elements/P1_3.h"
//#include "dolfin_elements/P2_3.h"
#include "dolfin_elements/vP1_3.h"
//#include "dolfin_elements/vP2_3.h"
#include "PeriodicBC.hpp"
#include "dolfin_helpers.hpp"
#include "xdmf_helpers.hpp"

//using namespace H5;

XDMFTetInterpol::XDMFTetInterpol(const std::string& infilename)
  : Interpol(infilename)
{

  // Input file (e.g. dolfin_params.dat)
  std::ifstream input(infilename);
  if (!input){
    std::cout << "File " << infilename <<" doesn't exist." << std::endl;
    exit(1);
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
  if (dolfin_params["periodic_z"] == "true"){
    periodic[2] = true;
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

  std::string cell_type_str = "tetrahedron";
  Uint gdim = 3;

  std::unique_ptr<dolfin::CellType> cell_type(dolfin::CellType::create(cell_type_str));

  std::vector<std::int64_t> coords_shape = dolfin::HDF5Interface::get_dataset_shape(meshfile.h5_id(), geometry_path);

  meshfile.read(mesh_in, topology_path, geometry_path, gdim, *cell_type, -1, coords_shape[0], false);

  mesh = std::make_shared<dolfin::Mesh>(mesh_in);
  dim = mesh->geometry().dim();
  assert(gdim == dim);
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

  // Precompute all tets Taylor-Hood P2-P1
  // FIXME compute on the fly and save
  tets_.resize(mesh->num_cells());
  dolfin_cells_.resize(mesh->num_cells());
  cell2cells_.resize(mesh->num_cells());

  for (std::size_t i = 0; i < mesh->num_cells(); ++i)
  {
    dolfin::Cell dolfin_cell(*mesh, i);
    tets_[i] = Tet(dolfin_cell);
    dolfin_cells_[i] = dolfin_cell;
  }
  // Build cell neighbour list for lookup speed
  build_neighbor_list(cell2cells_, mesh, dolfin_cells_);

  // Identify edge cells
  cell_type_.resize(mesh->num_cells());

  double tol = 1e-4; // 1e-12; // heuristic

  apply_periodic_boundaries(cell2cells_, periodic, x_min, x_max, mesh, dolfin_cells_, dim, tol);

  label_cell_type(cell_type_, cell2cells_, dim);

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
          for ( Uint k=0; k < static_cast<Uint>(dim); ++k)
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
  
  // Velocity
  u_space_ = std::make_shared<vP1_3::FunctionSpace>(mesh, constrained_domain);
  //u_space_ = std::make_shared<vP1_3::FunctionSpace>(mesh);
  ncoeffs_u = 4;

  // Pressure
  p_space_ = std::make_shared<P1_3::FunctionSpace>(mesh, constrained_domain);
  //p_space_ = std::make_shared<P1_3::FunctionSpace>(mesh);
  ncoeffs_p = 4;

  // const dolfin::GenericDofMap& u_dofmap = *u_space_->dofmap();
  auto xdof = p_space_->tabulate_dof_coordinates();

  std::map<std::tuple<double, double, double>, Uint> xmap;
  for ( Uint j=0; j < xdof.size()/dim; ++j ){
    xmap[{xdof[dim*j], xdof[dim*j+1], xdof[dim*j+2]}] = j;
  }

  std::vector<double> xdata;
  read_dataset_vector(h5filename_u, geometry_path, xdata, dim);

  j2i.resize(xdof.size()/dim);
  for ( Uint i=0; i < xdata.size()/dim; ++i ){
    Uint j = xmap[{xdata[dim*i], xdata[dim*i+1], xdata[dim*i+2]}];
    j2i[j] = i;
  }

  u_prev_data_.resize(xdof.size());
  u_next_data_.resize(xdof.size());
  p_prev_data_.resize(xdof.size()/dim);
  p_next_data_.resize(xdof.size()/dim);
  phi_prev_data_.resize(xdof.size()/dim);
  phi_next_data_.resize(xdof.size()/dim);

  check_dofs_fit(ncoeffs_u, ncoeffs_p, Tet::n_dofs_max, "XDMFTetInterpol");

  std::cout << "Setting max threads: " << omp_get_max_threads() << std::endl;

  found_.resize(omp_get_max_threads());

  // Precomputing dofs

  const dolfin::GenericDofMap& u_dofmap = *u_space_->dofmap();
  const dolfin::GenericDofMap& p_dofmap = *p_space_->dofmap();
  
  u_dofs_.build(u_dofmap, dolfin_cells_, "XDMFTetInterpol");
  u_dofs_.check_stride(3*ncoeffs_u, "XDMFTetInterpol");

  p_dofs_.build(p_dofmap, dolfin_cells_, "XDMFTetInterpol");
  p_dofs_.check_stride(ncoeffs_p, "XDMFTetInterpol");
}

void XDMFTetInterpol::update(const double t)
{
  MultiStampPair sp = ts.get(t);

  // Always load once; keep last bracket past t_max
  if ( !is_initialized || ((t_prev != sp.prev.t || t_next != sp.next.t) && t < ts.get_t_max()) )
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
      read_dataset_vector(u_path_prev[0], u_path_prev[1], data_, dim);
      reorder_indices(u_prev_data_, data_, j2i, dim);

      if (include_pressure){
        auto p_path_prev = ts.get_path("p", sp.prev.it);
        read_dataset_scalar(p_path_prev[0], p_path_prev[1], data_);
        reorder_indices(p_prev_data_, data_, j2i, 1);
      }

      if (include_phi){
        auto phi_path_prev = ts.get_path("phi", sp.prev.it);
        read_dataset_scalar(phi_path_prev[0], phi_path_prev[1], data_);
        reorder_indices(phi_prev_data_, data_, j2i, 1);
      }
    }

    auto u_path_next = ts.get_path("u", sp.next.it);
    std::cout << "Next: Timestep = " << sp.next.t << ", file = " << u_path_next[0] << ":" << u_path_next[1] << std::endl;
    read_dataset_vector(u_path_next[0], u_path_next[1], data_, dim);
    reorder_indices(u_next_data_, data_, j2i, dim);

    if (include_pressure){
      auto p_path_next = ts.get_path("p", sp.next.it);
      read_dataset_scalar(p_path_next[0], p_path_next[1], data_);
      reorder_indices(p_next_data_, data_, j2i, 1);
    }

    if (include_phi){
      auto phi_path_next = ts.get_path("phi", sp.next.it);
      read_dataset_scalar(phi_path_next[0], phi_path_next[1], data_);
      reorder_indices(phi_next_data_, data_, j2i, 1);
    }

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  t_update = t;

}



Vector3d XDMFTetInterpol::_modx(const Vector3d &x){
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


bool XDMFTetInterpol::locate(const Vector3d &x, const double t, CellPos& pos)
{
  assert(t <= t_next && t >= t_prev);
  const Vector3d xx = _modx(x);
  return locate_in_cells(tets_, cell2cells_, *mesh, dim, xx, pos,
                         found_);
}

void XDMFTetInterpol::evaluate(const Vector3d &x, const double tin, const CellPos& pos, PointValues& fields)
{
  const double _alpha_t = stamp_weight(tin, t_prev, t_next);

  // Compute Pk-Pl basis at x
  const int id = pos.id;
  const double r1 = pos.bary[0], r2 = pos.bary[1], r3 = pos.bary[2], r4 = pos.bary[3];

  std::array<double, Tet::n_dofs_max> _Nu_;
  std::array<double, Tet::n_dofs_max> _Np_;
  std::array<double, Tet::n_dofs_max> _Nux_;
  std::array<double, Tet::n_dofs_max> _Nuy_;
  std::array<double, Tet::n_dofs_max> _Nuz_;

  tets_[id].linearbasis(r1, r2, r3, r4, _Nu_.data());

  if (include_pressure){
    tets_[id].linearbasis(r1, r2, r3, r4, _Np_.data());
  }

  std::array<double, Tet::n_dofs_max*3> u_prev_block;
  std::array<double, Tet::n_dofs_max*3> u_next_block;

  // Restrict solution to cell
  //const dolfin::GenericDofMap& u_dofmap = *u_space_->dofmap();
  //auto u_dofs = u_dofmap.cell_dofs(dolfin_cells_[id].index());

  const std::uint32_t* u_dofs = u_dofs_[id];
  for (std::size_t i = 0; i < u_dofs_.stride(); ++i){
      u_prev_block[i] = u_prev_data_[u_dofs[i]];
      u_next_block[i] = u_next_data_[u_dofs[i]];
  }

  // Evaluate
  Vector3d U_prev = {std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, u_prev_block.begin(), 0.0),
                     std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, &u_prev_block[ncoeffs_u], 0.0),
                     std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, &u_prev_block[2*ncoeffs_u], 0.0)};
  Vector3d U_next = {std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, u_next_block.begin(), 0.0),
                     std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, &u_next_block[ncoeffs_u], 0.0),
                     std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, &u_next_block[2*ncoeffs_u], 0.0)};

  // Gradient
  Matrix3d gradU_prev = Matrix3d::Zero(), gradU_next = Matrix3d::Zero();
  if (wants_gradient()){
    tets_[id].linearderiv(r1, r2, r3, r4, _Nux_.data(), _Nuy_.data(), _Nuz_.data());

    gradU_prev <<
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, u_prev_block.begin(), 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, u_prev_block.begin(), 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, u_prev_block.begin(), 0.0),
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, &u_prev_block[ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, &u_prev_block[ncoeffs_u], 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, &u_prev_block[ncoeffs_u], 0.0),
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, &u_prev_block[2*ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, &u_prev_block[2*ncoeffs_u], 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, &u_prev_block[2*ncoeffs_u], 0.0);
    gradU_next << 
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, u_next_block.begin(), 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, u_next_block.begin(), 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, u_next_block.begin(), 0.0),
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, &u_next_block[ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, &u_next_block[ncoeffs_u], 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, &u_next_block[ncoeffs_u], 0.0),
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, &u_next_block[2*ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, &u_next_block[2*ncoeffs_u], 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, &u_next_block[2*ncoeffs_u], 0.0);
  }


  // Update
  fields.U = _alpha_t * U_next + (1-_alpha_t) * U_prev;
  fields.A = stamp_rate(U_next, U_prev, t_prev, t_next);

  if (wants_gradient()){
    fields.gradU = _alpha_t * gradU_next + (1-_alpha_t) * gradU_prev;
    fields.gradA = stamp_rate(gradU_next, gradU_prev, t_prev, t_next);
  }

  if (include_pressure){
    std::array<double, Tet::n_dofs_max> p_prev_block;
    std::array<double, Tet::n_dofs_max> p_next_block;
    
    const std::uint32_t* p_dofs = p_dofs_[id];
    for (std::size_t i = 0; i < p_dofs_.stride(); ++i){
        p_prev_block[i] = p_prev_data_[p_dofs[i]];
        p_next_block[i] = p_next_data_[p_dofs[i]];
    }
    
    // Evaluate
    double P_prev = std::inner_product(_Np_.data(), _Np_.data()+ncoeffs_p, p_prev_block.begin(), 0.0);
    double P_next = std::inner_product(_Np_.data(), _Np_.data()+ncoeffs_p, p_next_block.begin(), 0.0);
    fields.P = _alpha_t * P_next + (1-_alpha_t) * P_prev;
  }

  if (include_phi){
    std::array<double, Tet::n_dofs_max> phi_prev_block;
    std::array<double, Tet::n_dofs_max> phi_next_block;
    
    const std::uint32_t* p_dofs = p_dofs_[id];
    for (std::size_t i = 0; i < p_dofs_.stride(); ++i){
        phi_prev_block[i] = phi_prev_data_[p_dofs[i]];
        phi_next_block[i] = phi_next_data_[p_dofs[i]];
    }

    double Phi_prev = std::inner_product(_Np_.data(), _Np_.data()+ncoeffs_p, phi_prev_block.begin(), 0.0);
    double Phi_next = std::inner_product(_Np_.data(), _Np_.data()+ncoeffs_p, phi_next_block.begin(), 0.0);
    fields.Phi = _alpha_t * Phi_next + (1-_alpha_t) * Phi_prev;
  }

  // cell_type
  fields.cell_type = cell_type_[id];
}

#endif
