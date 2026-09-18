#ifdef USE_DOLFIN
#include "geometry.hpp"
#include "loader_params.hpp"
#include "dolfin_spaces.hpp"
#include "p12_eval.hpp"
#include "TriangleFreqInterpol.hpp"
#include <array>
#include <numeric>
#include "FreqStamps.hpp"
//#include "H5Cpp.h"
#include <cassert>
#include "PeriodicBC.hpp"
#include "dolfin_helpers.hpp"

TriangleFreqInterpol::TriangleFreqInterpol(const std::string& infilename)
  : MeshInterpol<Triangle>(infilename)
{
  dolfin_params = partrac::parse_file_or_exit(triangle_freq_schema(), infilename);

  // base frequency
  double tau = dolfin_params.get<double>("tau");
  if (tau > 0)
    omega0 = 2 * M_PI / tau;

  std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  // Using the FreqStamps class to hold frequency data
  fs.initialize(get_folder() + "/" + dolfin_params.get<std::string>("freqstamps"));

  read_mesh_params();
  std::string meshfilename = get_folder() + "/" + dolfin_params.get<std::string>("mesh");
  dolfin::HDF5File meshfile(MPI_COMM_WORLD, meshfilename, "r");

  dolfin::Mesh mesh_in;
  meshfile.read(mesh_in, "mesh", false);

  mesh = std::make_shared<dolfin::Mesh>(mesh_in);
  init_mesh_geometry();

  auto constrained_domain = std::make_shared<PeriodicBC>(periodic, x_min, x_max, dim);
  std::cout << "Made periodic domain." << std::endl;

  std::string u_el = dolfin_params.get<std::string>("velocity_space");
  std::string p_el = dolfin_params.get<std::string>("pressure_space");
  
  taylor_hood_spaces<Triangle>(u_el, p_el, include_pressure, mesh, constrained_domain,
                         u_space_, p_space_, ncoeffs_u, ncoeffs_p);

  // Precompute all triangles Taylor-Hood P2-P1
  // FIXME compute on the fly and save

  build_cells(*u_space_->dofmap(), true);

  std::cout << "Built neighbour list" << std::endl;

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
    file_i.read(*u_, dolfin_params.get<std::string>("velocity_field"));
    if (include_pressure)
      file_i.read(*p_, dolfin_params.get<std::string>("pressure_field"));

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

  std::array<double, Triangle::n_dofs_max> Nu_;
  std::array<double, Triangle::n_dofs_max> Np_;
  std::array<double, Triangle::n_dofs_max> Nux_;
  std::array<double, Triangle::n_dofs_max> Nuy_;

  cell_basis(cells_[id], pos.bary, ncoeffs_u, Nu_.data(), "u");
  if (include_pressure)
    cell_basis(cells_[id], pos.bary, ncoeffs_p, Np_.data(), "p");

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
    cell_deriv(cells_[id], pos.bary, ncoeffs_u, Nux_.data(), Nuy_.data(), nullptr, "u");

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

#endif
