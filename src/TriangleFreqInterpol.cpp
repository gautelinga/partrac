#ifdef USE_DOLFIN
#include "TriangleFreqInterpol.hpp"
#include <array>
#include <numeric>
#include "FreqStamps.hpp"
//#include "H5Cpp.h"
#include <boost/algorithm/string.hpp>
#include <cassert>
#include "dolfin_elements/P1_2.h"
#include "dolfin_elements/P2_2.h"
#include "dolfin_elements/vP1_2.h"
#include "dolfin_elements/vP2_2.h"
#include "PeriodicBC.hpp"
#include "dolfin_helpers.hpp"

TriangleFreqInterpol::TriangleFreqInterpol(const std::string& infilename)
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

  // base frequency
  double tau = stod(dolfin_params["tau"]);
  if (tau > 0)
    omega0 = 2 * M_PI / tau;

  std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  // Using the FreqStamps class to hold frequency data
  fs.initialize(get_folder() + "/" + dolfin_params["freqstamps"]);

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

  // Precompute all triangles Taylor-Hood P2-P1
  // FIXME compute on the fly and save
  triangles_.resize(mesh->num_cells());
  dolfin_cells_.resize(mesh->num_cells());
  cell2cells_.resize(mesh->num_cells());

  const std::vector<std::uint32_t> order =
    cell_order(*u_space_->dofmap(), mesh->num_cells(), dolfin_params["renumber_cells"], dolfin2local_);
  for (std::size_t l = 0; l < mesh->num_cells(); ++l)
  {
    dolfin::Cell dolfin_cell(*mesh, order[l]);
    triangles_[l] = Triangle(dolfin_cell);
    dolfin_cells_[l] = dolfin_cell;
  }
  // Build cell neighbour list for lookup speed
  build_neighbor_list(cell2cells_, mesh, dolfin_cells_,
                      dolfin2local_.empty() ? nullptr : &dolfin2local_);

  std::cout << "Built neighbour list" << std::endl;

  double tol = 1e-12; // heuristic
  apply_periodic_boundaries(cell2cells_, periodic, x_min, x_max, mesh, dolfin_cells_, dim, tol);

  // make structures
  u_ = std::make_shared<dolfin::Function>(u_space_);
  if (include_pressure)
    p_ = std::make_shared<dolfin::Function>(p_space_);

  u_coefficients_.resize(fs.size());
  if (include_pressure)
    p_coefficients_.resize(fs.size());
  
  
  for (std::size_t iFreq=0; iFreq < static_cast<std::size_t>(fs.size()); ++iFreq)
  {
    u_coefficients_[iFreq].resize(mesh->num_cells());
    if (include_pressure)
      p_coefficients_[iFreq].resize(mesh->num_cells());

    FreqStamp& f = fs.get(iFreq);
    dolfin::HDF5File file_i(MPI_COMM_WORLD, get_folder() + "/" + f.filename, "r");
    file_i.read(*u_, "u");
    if (include_pressure)
      file_i.read(*p_, "p");

    // Reused buffers
    std::vector<double> coordinate_dofs;
    ufc::cell ufc_cell;
    for (std::size_t id = 0; id < mesh->num_cells(); ++id)
    {
      dolfin_cells_[id].get_coordinate_dofs(coordinate_dofs);
      dolfin_cells_[id].get_cell_data(ufc_cell);
      u_coefficients_[iFreq][id].resize(dim*ncoeffs_u);
      u_->restrict(u_coefficients_[iFreq][id].data(), *u_space_->element(), dolfin_cells_[id],
                   coordinate_dofs.data(), ufc_cell);
      if (include_pressure){
        p_coefficients_[iFreq][id].resize(ncoeffs_p);
        p_->restrict(p_coefficients_[iFreq][id].data(), *p_space_->element(), dolfin_cells_[id],
                     coordinate_dofs.data(), ufc_cell);
      }
    }
  }

  check_dofs_fit(ncoeffs_u, ncoeffs_p, Triangle::n_dofs_max, "TriangleFreqInterpol");

  std::cout << "Setting max threads: " << omp_get_max_threads() << std::endl;

  found_.resize(omp_get_max_threads());

  // todo: remove below
}

void TriangleFreqInterpol::update(const double t)
{
  //std::cout << "Update t = " << t << std::endl;
  if (!is_initialized){
    is_initialized = true;
  }
  t_update = t;
}




Vector3d TriangleFreqInterpol::_modx(const Vector3d &x){
  Vector3d x_loc = x;
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

bool TriangleFreqInterpol::locate(const Vector3d &x, const double t, CellPos& pos)
{
  const Vector3d xx = _modx(x);
  return locate_in_cells(triangles_, cell2cells_, *mesh, dim, xx, pos,
                         found_, dolfin2local_.empty() ? nullptr : &dolfin2local_);
}

void TriangleFreqInterpol::evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  const int id = pos.id;
  // Assume found in fluid domain
  // std::cout << "probing..." << std::endl;

  // update frequency weights
  static thread_local std::vector<double> w_f_; w_f_.resize(fs.size());
  static thread_local std::vector<double> wt_f_; wt_f_.resize(fs.size());
  for (std::size_t iFreq=0; iFreq < static_cast<std::size_t>(fs.size()); ++iFreq){
    FreqStamp& f = fs.get(iFreq);
    double a = f.a; //fs.get_a(iFreq);
    double t_shift = f.t; // fs.get_t(iFreq);
    w_f_[iFreq] = a * cos(omega0 * (iFreq * t + t_shift));
    wt_f_[iFreq] = - a * omega0 * iFreq * sin(omega0 * (iFreq * t + t_shift));
  }

  // Compute Pk-Pl basis at x
  const double r1 = pos.bary[0], r2 = pos.bary[1], r3 = pos.bary[2];

  std::array<double, Triangle::n_dofs_max> Nu_;
  std::array<double, Triangle::n_dofs_max> Np_;
  std::array<double, Triangle::n_dofs_max> Nux_;
  std::array<double, Triangle::n_dofs_max> Nuy_;

  if (ncoeffs_u == 3){
    triangles_[id].linearbasis(r1, r2, r3, Nu_.data());
  }
  else if (ncoeffs_u == 6){
    triangles_[id].quadbasis(r1, r2, r3, Nu_.data());
  }
  else {
    std::cout << "Unrecognized ncoeffs_u = " << ncoeffs_u << std::endl;
    exit(1);
  }
  if (include_pressure){
    if (ncoeffs_p == 3){
      triangles_[id].linearbasis(r1, r2, r3, Np_.data());
    }
    else if (ncoeffs_p == 6){
      triangles_[id].quadbasis(r1, r2, r3, Np_.data());
    }
    else {
      std::cout << "Unrecognized ncoeffs_p = " << ncoeffs_p << std::endl;
      exit(1);
    }
  }

  static thread_local std::vector<double> ux_f_; ux_f_.resize(fs.size());
  static thread_local std::vector<double> uy_f_; uy_f_.resize(fs.size());
  static thread_local std::vector<double> p_f_; p_f_.resize(fs.size());

  for (std::size_t iFreq=0; iFreq < static_cast<std::size_t>(fs.size()); ++iFreq){
    ux_f_[iFreq] = std::inner_product(Nu_.data(), Nu_.data()+ncoeffs_u, u_coefficients_[iFreq][id].begin(), 0.0);
    uy_f_[iFreq] = std::inner_product(Nu_.data(), Nu_.data()+ncoeffs_u, &u_coefficients_[iFreq][id][ncoeffs_u], 0.0);
    if (include_pressure)
      p_f_[iFreq] = std::inner_product(Np_.data(), Np_.data()+ncoeffs_p, p_coefficients_[iFreq][id].begin(), 0.0);
  }

  // Update
  fields.U = { std::inner_product(w_f_.begin(), w_f_.end(), ux_f_.begin(), 0.0),
               std::inner_product(w_f_.begin(), w_f_.end(), uy_f_.begin(), 0.0),
               0.};
  fields.A = { std::inner_product(wt_f_.begin(), wt_f_.end(), ux_f_.begin(), 0.0),
               std::inner_product(wt_f_.begin(), wt_f_.end(), uy_f_.begin(), 0.0),
               0.};
  if (include_pressure){
    fields.P = std::inner_product(w_f_.begin(), w_f_.end(), p_f_.begin(), 0.0);
  }

  if (wants_gradient()){
    if (ncoeffs_u == 3){
      triangles_[id].linearderiv(r1, r2, r3, Nux_.data(), Nuy_.data());
    }
    else if (ncoeffs_u == 6){
      triangles_[id].quadderiv(r1, r2, r3, Nux_.data(), Nuy_.data());
    }
    else {
      std::cout << "Unrecognized ncoeffs_u = " << ncoeffs_u << std::endl;
      exit(1);
    }

    static thread_local std::vector<double> uxx_f_; uxx_f_.resize(fs.size());
    static thread_local std::vector<double> uxy_f_; uxy_f_.resize(fs.size());
    static thread_local std::vector<double> uyx_f_; uyx_f_.resize(fs.size());
    static thread_local std::vector<double> uyy_f_; uyy_f_.resize(fs.size());
  
    for (std::size_t iFreq=0; iFreq < static_cast<std::size_t>(fs.size()); ++iFreq){
      uxx_f_[iFreq] = std::inner_product(Nux_.data(), Nux_.data()+ncoeffs_u, u_coefficients_[iFreq][id].begin(), 0.0);
      uxy_f_[iFreq] = std::inner_product(Nuy_.data(), Nuy_.data()+ncoeffs_u, u_coefficients_[iFreq][id].begin(), 0.0);
      uyx_f_[iFreq] = std::inner_product(Nux_.data(), Nux_.data()+ncoeffs_u, &u_coefficients_[iFreq][id][ncoeffs_u], 0.0);
      uyy_f_[iFreq] = std::inner_product(Nuy_.data(), Nuy_.data()+ncoeffs_u, &u_coefficients_[iFreq][id][ncoeffs_u], 0.0);
    }

    fields.gradU(0, 0) = std::inner_product(w_f_.begin(), w_f_.end(), uxx_f_.begin(), 0.0);
    fields.gradU(0, 1) = std::inner_product(w_f_.begin(), w_f_.end(), uxy_f_.begin(), 0.0);
    fields.gradU(1, 0) = std::inner_product(w_f_.begin(), w_f_.end(), uyx_f_.begin(), 0.0);
    fields.gradU(1, 1) = std::inner_product(w_f_.begin(), w_f_.end(), uyy_f_.begin(), 0.0);

    fields.gradA(0, 0) = std::inner_product(wt_f_.begin(), wt_f_.end(), uxx_f_.begin(), 0.0);
    fields.gradA(0, 1) = std::inner_product(wt_f_.begin(), wt_f_.end(), uxy_f_.begin(), 0.0);
    fields.gradA(1, 0) = std::inner_product(wt_f_.begin(), wt_f_.end(), uyx_f_.begin(), 0.0);
    fields.gradA(1, 1) = std::inner_product(wt_f_.begin(), wt_f_.end(), uyy_f_.begin(), 0.0);
  }
}

void TriangleFreqInterpol::enable_reflection()
{
  build_facet_neighbours(facet_neigh_, mesh, dolfin_cells_,
                         dolfin2local_.empty() ? nullptr : &dolfin2local_,
                         periodic, x_min, x_max, dim, 1e-12);
  period_ = periodic_lengths(periodic, x_min, x_max, dim);
  can_reflect = true;
}

bool TriangleFreqInterpol::reflect(const Vector3d& x, Vector3d& dx, CellPos& pos)
{
  return reflect_in_cells(triangles_, facet_neigh_, 3, period_, x, dx, pos,
                          [this](const Vector3d& p){ return _modx(p); });
}

#endif
