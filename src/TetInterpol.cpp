#ifdef USE_DOLFIN
#include "TetInterpol.hpp"
#include <array>
#include <numeric>
#include "Timestamps.hpp"
#include "H5Cpp.h"
#include <boost/algorithm/string.hpp>
#include <cassert>
#include "dolfin_elements/P1_3.h"
#include "dolfin_elements/P2_3.h"
#include "dolfin_elements/vP1_3.h"
#include "dolfin_elements/vP2_3.h"
#include "PeriodicBC.hpp"
#include "dolfin_helpers.hpp"

using namespace H5;

TetInterpol::TetInterpol(const std::string& infilename)
  : Interpol(infilename)
{
  std::ifstream input(infilename);
  if (!input){
    std::cout << "File " << infilename <<" doesn't exist." << std::endl;
    exit(1);
  }

  // Default parameters
  dolfin_params["velocity_field"] = "u";
  dolfin_params["pressure_field"] = "p";
  dolfin_params["ignore_pressure"] = "false";

  // Overload from file
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

  // Precompute all tets Taylor-Hood P2-P1
  // FIXME compute on the fly and save
  tets_.resize(mesh->num_cells());
  dolfin_cells_.resize(mesh->num_cells());
  cell2cells_.resize(mesh->num_cells());

  auto constrained_domain = std::make_shared<PeriodicBC>(periodic, x_min, x_max, dim);

  // u_space_ = std::make_shared<vP2_3::FunctionSpace>(mesh, constrained_domain);
  // p_space_ = std::make_shared<P1_3::FunctionSpace>(mesh, constrained_domain);

  std::string u_el = dolfin_params["velocity_space"];
  std::string p_el = dolfin_params["pressure_space"];

  // Velocity
  if (u_el == "P1"){
    u_space_ = std::make_shared<vP1_3::FunctionSpace>(mesh, constrained_domain);
    ncoeffs_u = 4;
  }
  else if (u_el == "P2"){
    u_space_ = std::make_shared<vP2_3::FunctionSpace>(mesh, constrained_domain);
    ncoeffs_u = 10;
  }
  else {
    std::cout << "Unrecognized velocity element: " << u_el << std::endl;
    exit(1);
  }
  // Pressure
  if (include_pressure){
    if (p_el == "P1"){
      p_space_ = std::make_shared<P1_3::FunctionSpace>(mesh, constrained_domain);
      ncoeffs_p = 4;
    }
    else if (p_el == "P2"){
      p_space_ = std::make_shared<P2_3::FunctionSpace>(mesh, constrained_domain);
      ncoeffs_p = 10;
    }
    else {
      std::cout << "Unrecognized pressure element: " << p_el << std::endl;
      exit(1);
    }
  }
  else {
    std::cout << "Note: Ignoring pressure." << std::endl;
  }


  const std::size_t ncells = mesh->num_cells();
  const std::vector<std::uint32_t> order =
    cell_order(*u_space_->dofmap(), ncells, dolfin_params["renumber_cells"], dolfin2local_);
  for (std::size_t l = 0; l < ncells; ++l)
  {
    dolfin::Cell dolfin_cell(*mesh, order[l]);
    tets_[l] = Tet(dolfin_cell);
    dolfin_cells_[l] = dolfin_cell;
  }
  // Build cell neighbour list for lookup speed
  build_neighbor_list(cell2cells_, mesh, dolfin_cells_,
                      dolfin2local_.empty() ? nullptr : &dolfin2local_);

  double tol = 1e-12; // heuristic
  apply_periodic_boundaries(cell2cells_, periodic, x_min, x_max, mesh, dolfin_cells_, dim, tol);

  cell_type_.resize(mesh->num_cells());
  label_cell_type(cell_type_, cell2cells_, dim);

  cell_facets_.resize(mesh->num_cells());

  for ( Uint i=0; i < mesh->num_cells(); ++i)
  {
    if (cell_type_[i] == 1)
    {
      auto num_facets = dolfin_cells_[i].num_entities(dim-1);
      auto facets = dolfin_cells_[i].entities(dim-1);
      for ( std::size_t j = 0; j < num_facets; ++j ){
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
            std::vector<Vector3d> facet_loc;
            for (dolfin::VertexIterator v(dolfin_facet); !v.end(); ++v){
              Vector3d vloc(v->point().coordinates());
              facet_loc.push_back(vloc);
            }
            cell_facets_[i].push_back(facets_.size());
            facets_.push_back(facet_loc);
          }
        }
      }
    }
  }

  hmin = mesh->hmin();

  u_prev_ = std::make_shared<dolfin::Function>(u_space_);
  u_next_ = std::make_shared<dolfin::Function>(u_space_);

  // u_prev_coefficients_.resize(3*ncoeffs_u);
  // u_next_coefficients_.resize(3*ncoeffs_u);

  //Nu_.resize(ncoeffs_u);
  //Nux_.resize(ncoeffs_u);
  //Nuy_.resize(ncoeffs_u);
  //Nuz_.resize(ncoeffs_u);

  if (include_pressure){
    p_prev_ = std::make_shared<dolfin::Function>(p_space_);
    p_next_ = std::make_shared<dolfin::Function>(p_space_);

    //p_prev_coefficients_.resize(ncoeffs_p);
    //p_next_coefficients_.resize(ncoeffs_p);

    //Np_.resize(ncoeffs_p);
  }

  check_dofs_fit(ncoeffs_u, ncoeffs_p, Tet::n_dofs_max, "TetInterpol");

  std::cout << "Setting max threads: " << omp_get_max_threads() << std::endl;

  found_.resize(omp_get_max_threads());

  // Precomputing dofs
  const dolfin::GenericDofMap& u_dofmap = *u_space_->dofmap();

  u_dofs_.build(u_dofmap, dolfin_cells_, "TetInterpol");
  u_dofs_.check_stride(3*ncoeffs_u, "TetInterpol");

  if (include_pressure){
    const dolfin::GenericDofMap& p_dofmap = *p_space_->dofmap();

    p_dofs_.build(p_dofmap, dolfin_cells_, "TetInterpol");
    p_dofs_.check_stride(ncoeffs_p, "TetInterpol");
  }
  
  can_reflect = true;
}

void TetInterpol::update(const double t)
{

  StampPair sp = ts.get(t);
  // std::cout << sp.prev.filename << " " << sp.next.filename << std::endl;

  // Always load once; keep last bracket past t_max
  if ( !is_initialized || ((t_prev != sp.prev.t || t_next != sp.next.t) && t < ts.get_t_max()) ){

    if (is_initialized && t_next == sp.prev.t)
    {
      std::cout << "Prev: Timestep = " << sp.prev.t << ", swapping... "<< std::endl;
      u_prev_vec.swap(u_next_vec);
    }
    else {
      std::cout << "Prev: Timestep = " << sp.prev.t << ", filename = " << sp.prev.filename << std::endl;
      dolfin::HDF5File prevfile(MPI_COMM_WORLD, get_folder() + "/" + sp.prev.filename, "r");

      prevfile.read(*u_prev_, dolfin_params["velocity_field"]);
      u_prev_->vector()->get_local(u_prev_vec);
      if (include_pressure){
        prevfile.read(*p_prev_, dolfin_params["pressure_field"]);
        p_prev_->vector()->get_local(p_prev_vec);
      }
    }

    std::cout << "Next: Timestep = " << sp.next.t << ", filename = " << sp.next.filename << std::endl;
    // Single stamp: copy prev
    if (sp.next.filename == sp.prev.filename){
      u_next_vec = u_prev_vec;
      if (include_pressure)
        p_next_vec = p_prev_vec;
    }
    else {
      dolfin::HDF5File nextfile(MPI_COMM_WORLD, get_folder() + "/" + sp.next.filename, "r");
      nextfile.read(*u_next_, dolfin_params["velocity_field"]);
      u_next_->vector()->get_local(u_next_vec);
      if (include_pressure){
        nextfile.read(*p_next_, dolfin_params["pressure_field"]);
        p_next_->vector()->get_local(p_next_vec);
      }
    }

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  // alpha_t = sp.weight_next(t);
  t_update = t;
}



Vector3d TetInterpol::_modx(const Vector3d &x){
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


bool TetInterpol::locate(const Vector3d &x, const double t, CellPos& pos)
{
  const Vector3d xx = _modx(x);
  return locate_in_cells(tets_, cell2cells_, *mesh, dim, xx, pos,
                         found_, dolfin2local_.empty() ? nullptr : &dolfin2local_);
}

void TetInterpol::evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  // FIXME: better interpolation than using "restrict" (see Triangle)

  // Assuming inside fluid
  assert(t <= t_next && t >= t_prev);
  const double alpha_t = stamp_weight(t, t_prev, t_next);

  std::array<double, Tet::n_dofs_max> _Nu_;
  std::array<double, Tet::n_dofs_max> _Np_;
  std::array<double, Tet::n_dofs_max> _Nux_;
  std::array<double, Tet::n_dofs_max> _Nuy_;
  std::array<double, Tet::n_dofs_max> _Nuz_;

  // Compute P2-P1 basis at x
  const int id = pos.id;
  const double r1 = pos.bary[0], r2 = pos.bary[1], r3 = pos.bary[2], r4 = pos.bary[3];
  if (ncoeffs_u == 4){
    tets_[id].linearbasis(r1, r2, r3, r4, _Nu_.data());
  }
  else if (ncoeffs_u == 10){
    tets_[id].quadbasis(r1, r2, r3, r4, _Nu_.data());
  }
  else {
    std::cout << "Unrecognized ncoeffs_u = " << ncoeffs_u << std::endl;
    exit(1);
  }
  if (include_pressure){
    if (ncoeffs_p == 4){
      tets_[id].linearbasis(r1, r2, r3, r4, _Np_.data());
    }
    else if (ncoeffs_p == 10){
      tets_[id].quadbasis(r1, r2, r3, r4, _Np_.data());
    }
    else {
      std::cout << "Unrecognized ncoeffs_p = " << ncoeffs_p << std::endl;
      exit(1);
    }
  }

  std::array<double, Tet::n_dofs_max*3> u_prev_coefficients_;
  std::array<double, Tet::n_dofs_max*3> u_next_coefficients_;

  // Gathered: restrict() is not thread-safe

  const std::uint32_t* u_dofs = u_dofs_[id];
  for (std::size_t i=0; i < u_dofs_.stride(); ++i){
    u_prev_coefficients_[i] = u_prev_vec[u_dofs[i]];
    u_next_coefficients_[i] = u_next_vec[u_dofs[i]];
  }

  //std::cout << "VV " << vvec[0] << std::endl;

  //u_prev_->vector()->get_local(u_prev_coefficients_.data(), u_dofs_[id].size(), u_dofs_[id].data());
  //u_next_->vector()->get_local(u_next_coefficients_.data(), u_dofs_[id].size(), u_dofs_[id].data());

  //for (std::size_t i = 0; i < u_dofs_[id].size(); ++i){
  //  u_prev_block[i] = u_prev_data_[u_dofs_[id][i]];
  //  u_next_block[i] = u_next_data_[u_dofs_[id][i]];
  //}

  // Evaluate
  Vector3d U_prev = {std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, u_prev_coefficients_.begin(), 0.0),
    std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, &u_prev_coefficients_[1*ncoeffs_u], 0.0),
    std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, &u_prev_coefficients_[2*ncoeffs_u], 0.0)};
  Vector3d U_next = {std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, u_next_coefficients_.begin(), 0.0),
    std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, &u_next_coefficients_[1*ncoeffs_u], 0.0),
    std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, &u_next_coefficients_[2*ncoeffs_u], 0.0)};
  // else unrecognized element

  // Update
  fields.U = alpha_t * U_next + (1-alpha_t) * U_prev;
  fields.A = stamp_rate(U_next, U_prev, t_prev, t_next);

  if (include_pressure){
    std::array<double, Tet::n_dofs_max> p_prev_coefficients_;
    std::array<double, Tet::n_dofs_max> p_next_coefficients_;

    const std::uint32_t* p_dofs = p_dofs_[id];
    for (std::size_t i=0; i < p_dofs_.stride(); ++i){
      p_prev_coefficients_[i] = p_prev_vec[p_dofs[i]];
      p_next_coefficients_[i] = p_next_vec[p_dofs[i]];
    }

    // Evaluate
    double P_prev = std::inner_product(_Np_.data(), _Np_.data()+ncoeffs_p, p_prev_coefficients_.begin(), 0.0);
    double P_next = std::inner_product(_Np_.data(), _Np_.data()+ncoeffs_p, p_next_coefficients_.begin(), 0.0);
    // else unrecognized element

    // Update
    fields.P = alpha_t * P_next + (1-alpha_t) * P_prev;
  }
  else { // Unnecessary
    fields.P = 0.;
  }

  if (wants_gradient()){
    if (ncoeffs_u == 4){
      tets_[id].linearderiv(r1, r2, r3, r4, _Nux_.data(), _Nuy_.data(), _Nuz_.data());
    }
    else if (ncoeffs_u == 10){
      tets_[id].quadderiv(r1, r2, r3, r4, _Nux_.data(), _Nuy_.data(), _Nuz_.data());
    }
    else {
      std::cout << "Unrecognized ncoeffs_u = " << ncoeffs_u << std::endl;
      exit(1);
    }
    Matrix3d gradU_prev;
    gradU_prev << std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, u_prev_coefficients_.begin(), 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, u_prev_coefficients_.begin(), 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, u_prev_coefficients_.begin(), 0.0),
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, &u_prev_coefficients_[1*ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, &u_prev_coefficients_[1*ncoeffs_u], 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, &u_prev_coefficients_[1*ncoeffs_u], 0.0),
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, &u_prev_coefficients_[2*ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, &u_prev_coefficients_[2*ncoeffs_u], 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, &u_prev_coefficients_[2*ncoeffs_u], 0.0);
    Matrix3d gradU_next;
    gradU_next << std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, u_next_coefficients_.begin(), 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, u_next_coefficients_.begin(), 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, u_next_coefficients_.begin(), 0.0),
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, &u_next_coefficients_[1*ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, &u_next_coefficients_[1*ncoeffs_u], 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, &u_next_coefficients_[1*ncoeffs_u], 0.0),
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, &u_next_coefficients_[2*ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, &u_next_coefficients_[2*ncoeffs_u], 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, &u_next_coefficients_[2*ncoeffs_u], 0.0);

    // Update
    fields.gradU = alpha_t * gradU_next + (1-alpha_t) * gradU_prev;
    fields.gradA = stamp_rate(gradU_next, gradU_prev, t_prev, t_next);
  }
}

void TetInterpol::reflect(Vector3d &x, Vector3d &dx_new, const double t, const double dt, int& cell_id){
  if (cell_id == -1) locate(x, t, cell_id);
  //std::cout << "Reflecting: " << cell_id ; // << std::endl;
  //std::cout << " " << cell_type_[cell_id] << " ";
  //for (auto & facet : cell_facets_[cell_id]) {
  //  for (auto & pt : facet ){
  //    std::cout << pt << " ";
  //  }
    //std::cout << facet.size() << " ";
  //}
  //std::cout << std::endl;

  if (cell_type_[cell_id] != 1 || true){
    int cell_id_tmp = cell_id;
    Vector3d x_tmp = x;
    bool is_inside_tmp = true;
    int count = 0;
    int splits = 1;
    double ddx = 1e-3 * hmin;

    int num_splits = log(dx_new.norm()/ddx)/log(2);

    while (cell_type_[cell_id_tmp] != 1 || splits < num_splits)
    {
      //std::cout << count << " " << splits << " " << cell_type_[cell_id_tmp] << std::endl;
      Vector3d ddx_tmp = pow(0.5, splits) * dx_new;

      if ( (x_tmp+ddx_tmp-x).squaredNorm() > dx_new.squaredNorm() + 1e-7){
        //std::cout << "Failed..." << std::endl;
        exit(1);
        break;
      }

      is_inside_tmp = locate(x_tmp + ddx_tmp, t + pow(0.5, splits)*dt, cell_id_tmp);
      if (is_inside_tmp){
        x_tmp = x_tmp + ddx_tmp;
      }
      else {
        ++splits;
        //ddx_tmp /= 2;
      }
      ++count;
    }
    //std::cout << std::endl;
    cell_id = cell_id_tmp;
  }

  bool crossing = false;
  double beta;
  Vector3d N;

  for (auto & facet_id : cell_facets_[cell_id]) {
    auto facet = facets_[facet_id];
    crossing = _cross_facet(beta, N, x, dx_new, facet);
    if (crossing) break;
  }
  if (!crossing){
    for (auto & facet : facets_){
      crossing = _cross_facet(beta, N, x, dx_new, facet);
      if (crossing) break;
    }
  }
  if (crossing){
    //Vector3d xc = x + beta * dx_new;
    //Vector3d xpp = x + dx_new - 2*N*N.dot(dx_new) * (1-beta);
    //Vector3d xmm = x + 2*N*N.dot(dx_new) * beta;
    
    Vector3d dx_out = dx_new - 2*N*N.dot(dx_new) * (1-beta);

    dx_new = dx_out;
    // TODO: multiple reflections

  }

  //return {0., 0., 0.};
}

bool TetInterpol::_cross_facet(double& beta, Vector3d& N, const Vector3d& x, const Vector3d &dx_new, std::vector<Vector3d> &facet)
{
  Vector3d p0 = (x_max + x_min)/2; // To avoid periodic effects
  Vector3d p10 = facet[1]-facet[0];
  Vector3d p20 = facet[2]-facet[0];
  Vector3d s0 = _modx(x - facet[0] + p0);
  Vector3d s10 = dx_new;
  Vector3d sp0 = s0-p0;

  Matrix3d A;
  A.col(0) = -s10;
  A.col(1) = p10;
  A.col(2) = p20;

  Vector3d tuv = A.inverse() * sp0;
  bool crossing = (tuv[0] >= 0 && tuv[0] <= 1) && (tuv[1] >= 0 && tuv[1] <= 1) && (tuv[2] >= 0 && tuv[2] <= 1);

  beta = tuv[0];

  if (crossing){
    N = p10.cross(p20);
    N /= N.norm();
  }

  return crossing;
}

#endif
