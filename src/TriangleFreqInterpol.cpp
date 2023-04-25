#ifdef USE_DOLFIN
#include "TriangleFreqInterpol.hpp"
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

  // base frequency
  omega0 = 2 * M_PI / stod(dolfin_params["tau"]);

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

  Nu_.resize(ncoeffs_u);
  Np_.resize(ncoeffs_p);

  Nux_.resize(ncoeffs_u);
  Nuy_.resize(ncoeffs_u);

  // make structures
  for (std::size_t iFreq=0; iFreq < fs.size(); ++iFreq){
    u__.push_back(std::make_shared<dolfin::Function>(u_space_));

    std::vector<double> u_coeff_i(dim*ncoeffs_u);
    u_coefficients__.push_back(u_coeff_i);

    if (include_pressure){
      p__.push_back(std::make_shared<dolfin::Function>(p_space_));
      
      std::vector<double> p_coeff_i(ncoeffs_p);
      p_coefficients__.push_back(p_coeff_i);  
    }
  }

  w_f_.resize(fs.size());
  wt_f_.resize(fs.size());

  ux_f_.resize(fs.size());
  uy_f_.resize(fs.size());
  p_f_.resize(fs.size());

  uxx_f_.resize(fs.size());
  uxy_f_.resize(fs.size());
  uyx_f_.resize(fs.size());
  uyy_f_.resize(fs.size());

  // read files
  for (std::size_t iFreq=0; iFreq < fs.size(); ++iFreq){
    FreqStamp& f = fs.get(iFreq);
    std::string filename = f.filename; // fs.get_filename(iFreq);
    double a = f.a; // fs.get_a(iFreq);
    double t_shift = f.t; // fs.get_t(iFreq);

    std::cout << "FreqStamp: iFreq = " << iFreq << ", t_shift = " << t_shift << ", amplitude = " << a << ", filename = " << filename << std::endl;

    dolfin::HDF5File file_i(MPI_COMM_WORLD, get_folder() + "/" + filename, "r");
    file_i.read(*(u__[iFreq]), "u");
    if (include_pressure)
      file_i.read(*(p__[iFreq]), "p");
  }
}

void TriangleFreqInterpol::update(const double t)
{
  //std::cout << "Update t = " << t << std::endl;
  if (!is_initialized){
    is_initialized = true;
  }
  t_update = t;
}

void TriangleFreqInterpol::probe(const Vector3d &x, const double t)
{
  int id_prev = -1;
  probe(x, t, id_prev);
}

void TriangleFreqInterpol::probe(const Vector3d &x, const double t, int& id_prev)
{
  // std::cout << "probing..." << std::endl;

  // update frequency weights
  for (std::size_t iFreq=0; iFreq < fs.size(); ++iFreq){
    FreqStamp& f = fs.get(iFreq);
    double a = f.a; //fs.get_a(iFreq);
    double t_shift = f.t; // fs.get_t(iFreq);
    w_f_[iFreq] = a * cos(omega0 * (iFreq * t + t_shift));
    wt_f_[iFreq] = - a * omega0 * iFreq * sin(omega0 * (iFreq * t + t_shift));
  }

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
    double r1, r2, r3;
    triangles_[id].xy2bary(x_loc[0], x_loc[1], r1, r2, r3);

    if (ncoeffs_u == 3){
      triangles_[id].linearbasis(r1, r2, r3, Nu_);
    }
    else if (ncoeffs_u == 6){
      triangles_[id].quadbasis(r1, r2, r3, Nu_);
    }
    else {
      std::cout << "Unrecognized ncoeffs_u = " << ncoeffs_u << std::endl;
      exit(0);
    }
    if (include_pressure){
      if (ncoeffs_p == 3){
        triangles_[id].linearbasis(r1, r2, r3, Np_);
      }
      else if (ncoeffs_p == 6){
        triangles_[id].quadbasis(r1, r2, r3, Np_);
      }
      else {
        std::cout << "Unrecognized ncoeffs_p = " << ncoeffs_p << std::endl;
        exit(0);
      }
    }

    // Restricting
    //std::cout << "Restricting..." << std::endl;

    for (std::size_t iFreq=0; iFreq < fs.size(); ++iFreq){
      // Restrict solution to cell
      u__[iFreq]->restrict(u_coefficients__[iFreq].data(), *u_space_->element(), dolfin_cells_[id],
                        coordinate_dofs_[id].data(), ufc_cells_[id]);
      if (include_pressure){
        p__[iFreq]->restrict(p_coefficients__[iFreq].data(), *p_space_->element(), dolfin_cells_[id],
                          coordinate_dofs_[id].data(), ufc_cells_[id]);
      }
    }

    // Evaluate
    //std::cout << "Evaluating U, A, P..." << std::endl;

    for (std::size_t iFreq=0; iFreq < fs.size(); ++iFreq){
      ux_f_[iFreq] = std::inner_product(Nu_.begin(), Nu_.end(), u_coefficients__[iFreq].begin(), 0.0);
      uy_f_[iFreq] = std::inner_product(Nu_.begin(), Nu_.end(), &u_coefficients__[iFreq][ncoeffs_u], 0.0);
      if (include_pressure)
        p_f_[iFreq] = std::inner_product(Np_.begin(), Np_.end(), p_coefficients__[iFreq].begin(), 0.0);
    }

    // Update
    //std::cout << "Updating U, A, P..." << std::endl;

    U = {std::inner_product(w_f_.begin(), w_f_.end(), ux_f_.begin(), 0.0),
         std::inner_product(w_f_.begin(), w_f_.end(), uy_f_.begin(), 0.0),
         0.};
    A = {std::inner_product(wt_f_.begin(), wt_f_.end(), ux_f_.begin(), 0.0),
         std::inner_product(wt_f_.begin(), wt_f_.end(), uy_f_.begin(), 0.0),
         0.};
    if (include_pressure){
      P = std::inner_product(w_f_.begin(), w_f_.end(), p_f_.begin(), 0.0);
    }

    if (this->int_order > 1){
      if (ncoeffs_u == 3){
        triangles_[id].linearderiv(r1, r2, r3, Nux_, Nuy_);
      }
      else if (ncoeffs_u == 6){
        triangles_[id].quadderiv(r1, r2, r3, Nux_, Nuy_);
      }

      for (std::size_t iFreq=0; iFreq < fs.size(); ++iFreq){
        uxx_f_[iFreq] = std::inner_product(Nux_.begin(), Nux_.end(), u_coefficients__[iFreq].begin(), 0.0);
        uxy_f_[iFreq] = std::inner_product(Nuy_.begin(), Nuy_.end(), u_coefficients__[iFreq].begin(), 0.0);
        uyx_f_[iFreq] = std::inner_product(Nux_.begin(), Nux_.end(), &u_coefficients__[iFreq][ncoeffs_u], 0.0);
        uyy_f_[iFreq] = std::inner_product(Nuy_.begin(), Nuy_.end(), &u_coefficients__[iFreq][ncoeffs_u], 0.0);
      }

      //std::cout << "Updating gradU, gradA..." << std::endl;
      //Matrix3d gradU_, gradA_;
      /*gradU = {
        std::inner_product(w_f_.begin(), w_f_.end(), uxx_f_.begin(), 0.0),
        std::inner_product(w_f_.begin(), w_f_.end(), uxy_f_.begin(), 0.0),
        0.0,
        std::inner_product(w_f_.begin(), w_f_.end(), uyx_f_.begin(), 0.0),
        std::inner_product(w_f_.begin(), w_f_.end(), uyy_f_.begin(), 0.0),
        0.0,
        0.0,
        0.0,
        0.0};

      gradA = {
        std::inner_product(wt_f_.begin(), wt_f_.end(), uxx_f_.begin(), 0.0),
        std::inner_product(wt_f_.begin(), wt_f_.end(), uxy_f_.begin(), 0.0),
        0.0,
        std::inner_product(wt_f_.begin(), wt_f_.end(), uyx_f_.begin(), 0.0),
        std::inner_product(wt_f_.begin(), wt_f_.end(), uyy_f_.begin(), 0.0),
        0.0,
        0.0,
        0.0,
        0.0 };
      */
      gradU(0, 0) = std::inner_product(w_f_.begin(), w_f_.end(), uxx_f_.begin(), 0.0);
      gradU(0, 1) = std::inner_product(w_f_.begin(), w_f_.end(), uxy_f_.begin(), 0.0);
      gradU(1, 0) = std::inner_product(w_f_.begin(), w_f_.end(), uyx_f_.begin(), 0.0);
      gradU(1, 1) = std::inner_product(w_f_.begin(), w_f_.end(), uyy_f_.begin(), 0.0);

      gradA(0, 0) = std::inner_product(wt_f_.begin(), wt_f_.end(), uxx_f_.begin(), 0.0);
      gradA(0, 1) = std::inner_product(wt_f_.begin(), wt_f_.end(), uxy_f_.begin(), 0.0);
      gradA(1, 0) = std::inner_product(wt_f_.begin(), wt_f_.end(), uyx_f_.begin(), 0.0);
      gradA(1, 1) = std::inner_product(wt_f_.begin(), wt_f_.end(), uyy_f_.begin(), 0.0);
      //gradU = gradU_;
      //gradA = gradA_;

    }
  }
}
#endif
