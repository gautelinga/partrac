#include "Interpol.hpp"
#include "Timestamps.hpp"
#include "H5Cpp.h"

#ifndef __STRUCTUREDINTERPOL_HPP
#define __STRUCTUREDINTERPOL_HPP

static void fix_boundary(const std::map<std::array<Uint, 3>, std::vector<std::array<Uint, 3>>>& solid_ids_neigh,
                         double*** ux_next,
                         double*** uy_next,
                         double*** uz_next,
                         double*** rho_next,
                         double*** p_next) {
  for (auto &el : solid_ids_neigh){
    auto i = el.first;
    auto neigh = el.second;
    Uint num_neigh = neigh.size();
    if (num_neigh > 0){
      double ux_next_sum = 0.;
      double uy_next_sum = 0.;
      double uz_next_sum = 0.;
      double rho_next_sum = 0.;
      double p_next_sum = 0.;
      for (auto &i_ : neigh){
        ux_next_sum += ux_next[i_[0]][i_[1]][i_[2]];
        uy_next_sum += uy_next[i_[0]][i_[1]][i_[2]];
        uz_next_sum += uz_next[i_[0]][i_[1]][i_[2]];
        rho_next_sum += rho_next[i_[0]][i_[1]][i_[2]];
        p_next_sum += p_next[i_[0]][i_[1]][i_[2]];
      }
      ux_next[i[0]][i[1]][i[2]] = ux_next_sum/num_neigh;
      uy_next[i[0]][i[1]][i[2]] = uy_next_sum/num_neigh;
      uz_next[i[0]][i[1]][i[2]] = uz_next_sum/num_neigh;
      rho_next[i[0]][i[1]][i[2]] = rho_next_sum/num_neigh;
      p_next[i[0]][i[1]][i[2]] = p_next_sum/num_neigh;
    }
  }
}

static void load_field(H5::H5File &h5file
                     , double*** u
                     , const std::string field
                     , const int nx, const int ny, const int nz
                     ){
  H5::DataSet dset = h5file.openDataSet(field);
  H5::DataSpace dspace = dset.getSpace();
  std::vector<double> Uv(nx*ny*nz);
  dset.read(Uv.data(), PredType::NATIVE_DOUBLE, dspace, dspace);
  for (int ix=0; ix<nx; ++ix){
    for (int iy=0; iy<ny; ++iy){
      for (int iz=0; iz<nz; ++iz){
  	    u[ix][iy][iz] = Uv[nx*ny*iz+nx*iy+ix];
      }
    }
  }
}

static void load_int_field(H5::H5File &h5file
                         , int*** u
                         , const std::string field
                         , const int nx, const int ny, const int nz
                         )
{
  H5::DataSet dset = h5file.openDataSet(field);
  H5::DataSpace dspace = dset.getSpace();
  std::vector<int> Uv(nx*ny*nz);
  dset.read(Uv.data(), PredType::NATIVE_INT, dspace, dspace);
  for (int ix=0; ix<nx; ++ix){
    for (int iy=0; iy<ny; ++iy){
      for (int iz=0; iz<nz; ++iz){
  	u[ix][iy][iz] = Uv[nx*ny*iz+nx*iy+ix];
      }
    }
  }
}

static void load_h5(const std::string h5filename
                  , double*** ux
                  , double*** uy
                  , double*** uz
                  , double*** rho
                  , double*** p
                  , const int nx
                  , const int ny
                  , const int nz
                  , const bool verbose
                  , const bool ignore_density
                  , const bool ignore_pressure
                  , const bool ignore_uz
       )
{
  // Assert that h5 file exists
  verify_file_exists(h5filename);
  if (verbose)
    std::cout << "Opening " << h5filename << std::endl;
  H5::H5File h5file(h5filename, H5F_ACC_RDONLY);
  load_field(h5file, ux, "u_x", nx, ny, nz);
  load_field(h5file, uy, "u_y", nx, ny, nz);
  if (!ignore_uz)
    load_field(h5file, uz, "u_z", nx, ny, nz);
  if (!ignore_density)
    load_field(h5file, rho, "density", nx, ny, nz);
  if (!ignore_pressure)
    load_field(h5file, p, "pressure", nx, ny, nz);
  h5file.close();
}

static double weighted_sum(double*** C,
                           const Uint ind[3][2],
                           const double w[2][2][2]){  
  double f = 0.0;
  for (Uint q0=0; q0<2; ++q0){
    for (Uint q1=0; q1<2; ++q1){
      for (Uint q2=0; q2<2; ++q2){
        f += C[ind[0][q0]][ind[1][q1]][ind[2][q2]]*w[q0][q1][q2];
      }
    }
  }
  return f;
}

class StructuredInterpol : public Interpol {
public:
  StructuredInterpol(const std::string infilename);
  void update(const double t);
  void probe(const Vector3d &x, const double t);
  bool inside_domain();
  bool inside_domain(const Vector3d &x);
  Uint get_nx() { return n[0]; };
  Uint get_ny() { return n[1]; };
  Uint get_nz() { return n[2]; };
  double get_ux();
  double get_uy();
  double get_uz();
  double get_ax();
  double get_ay();
  double get_az();
  double get_t_min() { return ts.get_t_min(); };
  double get_t_max() { return ts.get_t_max(); };
  double get_rho();
  double get_p();
  double get_uxx();
  double get_uxy();
  double get_uxz();
  double get_uyx();
  double get_uyy();
  double get_uyz();
  double get_uzx();
  double get_uzy();
  double get_uzz();
  Matrix3d get_grada();
protected:
  void probe_space(const Vector3d &x);
  void probe_grad();
  Timestamps ts;
  double t_prev = 0.;
  double t_next = 0.;
  double alpha_t;

  Uint n[3] = {0, 0, 0};
  Vector3d dx;

  double wq[3][2];
  double dwq[3][2];

  int*** isSolid;
  double*** levelZ;
  double*** ux_prev;
  double*** uy_prev;
  double*** uz_prev;
  double*** ux_next;
  double*** uy_next;
  double*** uz_next;
  double*** rho_prev;
  double*** rho_next;
  double*** p_prev;
  double*** p_next;

  Uint ind[3][2] = {{0, 0}, {0, 0}, {0, 0}};  // trilinear intp
  // Uint ind_pc[3] = {0, 0, 0};  // piecewise constant intp

  double w[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dw_x[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dw_y[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dw_z[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};

  std::map<std::string, std::string> felbm_params;

  std::map<std::array<Uint, 3>, std::vector<std::array<Uint, 3>>> solid_ids_neigh;
  double inside_domain_factor = 0.;
  double inside_domain_factor_x = 0.;
  double inside_domain_factor_y = 0.;
  double inside_domain_factor_z = 0.;
  double Ux = 0.;
  double Uy = 0.;
  double Uz = 0.;
  double Ax = 0.;
  double Ay = 0.;
  double Az = 0.;
  int boundary_mode = 0; // 0: Rounded extrapolation, 1: Sharp extrapolation, 2: Sharp simple
  bool ignore_density = false;
  bool ignore_pressure = false;
  bool ignore_uz = false;
};

StructuredInterpol::StructuredInterpol(const std::string infilename) : Interpol(infilename) {
  std::ifstream input(infilename);
  if (!input){
    std::cout << "File " << infilename <<" doesn't exist." << std::endl;
    exit(0);
  }
  // Default params
  felbm_params["timestamps"] = "timestamps.dat";
  felbm_params["is_solid_file"] = "output_is_solid.h5";
  felbm_params["boundary_mode"] = "rounded";

  felbm_params["ignore_pressure"] = "false";
  felbm_params["ignore_density"] = "false";
  felbm_params["ignore_uz"] = "false";

  // Overload from file
  size_t found;
  std::string key, val;
  for (std::string line; getline(input, line); ){
    found = line.find('=');
    if (found != std::string::npos){
      key = line.substr(0, found);
      val = line.substr(found+1);
      boost::trim(key);
      boost::trim(val);
      felbm_params[key] = val;
      //std::cout << "--> " << key << ": " << val << std::endl;
    }
  }
  std::cout << "Chosen parameters:" << std::endl;
  for (auto & k : felbm_params) {
    std::cout << k.first << ": " << k.second << std::endl;
  }

  if (felbm_params["ignore_pressure"] == "true" || felbm_params["ignore_pressure"] == "True"){
    ignore_pressure = true;
  }
  if (felbm_params["ignore_density"] == "true" || felbm_params["ignore_density"] == "True"){
    ignore_density = true;
  }
  if (felbm_params["ignore_uz"] == "true" || felbm_params["ignore_uz"] == "True"){
    ignore_uz = true;
  }

  std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));
  ts.initialize(get_folder() + "/" + felbm_params["timestamps"]);

  std::string solid_filename = get_folder() + "/" + felbm_params["is_solid_file"];
  verify_file_exists(solid_filename);

  if (felbm_params["boundary_mode"] == "rounded"){
    boundary_mode = 0;
  }
  else if (felbm_params["boundary_mode"] == "sharp"){
    boundary_mode = 1;
  }
  else if (felbm_params["boundary_mode"] == "simple"){
    boundary_mode = 2;
  }
  else {
    std::cout << "StructuredInterpol: Unknown boundary mode." << std::endl;
    exit(0);
  }

  H5::H5File solid_file(solid_filename, H5F_ACC_RDONLY);
  H5::DataSet dset_solid = solid_file.openDataSet("is_solid");
  H5::DataSpace dspace_solid = dset_solid.getSpace();

  hsize_t dims[3];
  dspace_solid.getSimpleExtentDims(dims, NULL);
  for (Uint i=0; i<3; ++i)
    n[i] = dims[i];
  x_min << 0., 0., 0.;
  x_max << n[0], n[1], n[2];

  dx << this->get_Lx()/n[0], this->get_Ly()/n[1], this->get_Lz()/n[2];
  for (Uint i=0; i<3; ++i){
    dwq[i][0] = -1./dx[i];
    dwq[i][1] =  1./dx[i];
  }

  // Create arrays
  isSolid = new int**[n[0]];
  levelZ = new double**[n[0]];
  ux_prev = new double**[n[0]];
  uy_prev = new double**[n[0]];
  uz_prev = new double**[n[0]];
  ux_next = new double**[n[0]];
  uy_next = new double**[n[0]];
  uz_next = new double**[n[0]];
  rho_prev = new double**[n[0]];
  rho_next = new double**[n[0]];
  p_prev = new double**[n[0]];
  p_next = new double**[n[0]];
  for (Uint ix=0; ix<n[0]; ++ix){
    isSolid[ix] = new int*[n[1]];
    levelZ[ix] = new double*[n[1]];
    ux_prev[ix] = new double*[n[1]];
    uy_prev[ix] = new double*[n[1]];
    uz_prev[ix] = new double*[n[1]];
    ux_next[ix] = new double*[n[1]];
    uy_next[ix] = new double*[n[1]];
    uz_next[ix] = new double*[n[1]];
    rho_prev[ix] = new double*[n[1]];
    rho_next[ix] = new double*[n[1]];
    p_prev[ix] = new double*[n[1]];
    p_next[ix] = new double*[n[1]];
    for (Uint iy=0; iy<n[1]; ++iy){
      isSolid[ix][iy] = new int[n[2]];
      levelZ[ix][iy] = new double[n[2]];
      ux_prev[ix][iy] = new double[n[2]];
      uy_prev[ix][iy] = new double[n[2]];
      uz_prev[ix][iy] = new double[n[2]];
      ux_next[ix][iy] = new double[n[2]];
      uy_next[ix][iy] = new double[n[2]];
      uz_next[ix][iy] = new double[n[2]];
      rho_prev[ix][iy] = new double[n[2]];
      rho_next[ix][iy] = new double[n[2]];
      p_prev[ix][iy] = new double[n[2]];
      p_next[ix][iy] = new double[n[2]];
    }
  }



  Uint nnlist[6][3];
  Uint nnum=0;

  Uint nnnlist[12][3];
  Uint nnnum=0;

  Uint nnnnlist[8][3];
  Uint nnnnum=0;

  for (int i=-1; i<=1; ++i){
    for (int j=-1; j<=1; ++j){
      for (int k=-1; k<=1; ++k){
        Uint sumabsijk = abs(i)+abs(j)+abs(k);
        if (sumabsijk == 1){
          nnlist[nnum][0] = i;
          nnlist[nnum][1] = j;
          nnlist[nnum][2] = k;
          ++nnum;
        }
        else if (sumabsijk == 2){
          nnnlist[nnnum][0] = i;
          nnnlist[nnnum][1] = j;
          nnnlist[nnnum][2] = k;
          ++nnnum;
        }
        else if (sumabsijk == 3){
          nnnnlist[nnnnum][0] = i;
          nnnnlist[nnnnum][1] = j;
          nnnnlist[nnnnum][2] = k;
          ++nnnnum;
        }
      }
    }
  }

  load_int_field(solid_file, isSolid, "is_solid", n[0], n[1], n[2]);
  // Smooth solid-liquid interface
  if (boundary_mode == 0){
    for (Uint ix=0; ix<n[0]; ++ix){
      for (Uint iy=0; iy<n[1]; ++iy){
        for (Uint iz=0; iz<n[2]; ++iz){
          levelZ[ix][iy][iz] = -2*isSolid[ix][iy][iz]+1;
          if (isSolid[ix][iy][iz]){
            std::vector<std::array<Uint, 3>> neigh;
            for (Uint inn=0; inn<6; ++inn){
              Uint ix_ = imodulo(ix+nnlist[inn][0], n[0]);
              Uint iy_ = imodulo(iy+nnlist[inn][1], n[1]);
              Uint iz_ = imodulo(iz+nnlist[inn][2], n[2]);
              if (!isSolid[ix_][iy_][iz_]){
                neigh.push_back({ix_, iy_, iz_});
              }
            }
            if (neigh.size() == 0){
              for (Uint innn=0; innn<12; ++innn){
                Uint ix_ = imodulo(ix+nnnlist[innn][0], n[0]);
                Uint iy_ = imodulo(iy+nnnlist[innn][1], n[1]);
                Uint iz_ = imodulo(iz+nnnlist[innn][2], n[2]);
                if (!isSolid[ix_][iy_][iz_]){
                  neigh.push_back({ix_, iy_, iz_});
                }
              }
            }
            if (neigh.size() == 0){
              for (Uint innnn=0; innnn<8; ++innnn){
                Uint ix_ = imodulo(ix+nnnnlist[innnn][0], n[0]);
                Uint iy_ = imodulo(iy+nnnnlist[innnn][1], n[1]);
                Uint iz_ = imodulo(iz+nnnnlist[innnn][2], n[2]);
                if (!isSolid[ix_][iy_][iz_]){
                  neigh.push_back({ix_, iy_, iz_});
                }
              }
            }
            std::pair<std::array<Uint, 3>, std::vector<std::array<Uint, 3>>> pp({ix, iy, iz}, neigh);
            solid_ids_neigh.insert(pp);
          }
        }
      }
    }
  }
  else { // if (boundary_mode == 1 || boundary_mode == 2){
    for (Uint ix=0; ix<n[0]; ++ix){
      for (Uint iy=0; iy<n[1]; ++iy){
        for (Uint iz=0; iz<n[2]; ++iz){
          levelZ[ix][iy][iz] = 1.0-isSolid[ix][iy][iz];
        }
      }
    }
  }
}

void StructuredInterpol::update(const double t){
  StampPair sp = ts.get(t);

  if (!is_initialized || t_prev != sp.prev.t || t_next != sp.next.t){
    if (is_initialized && t_next == sp.prev.t){
      std::swap(ux_prev, ux_next);
      std::swap(uy_prev, uy_next);
      if (!ignore_uz)
        std::swap(uz_prev, uz_prev);
      if (!ignore_density)
        std::swap(rho_prev, rho_next);
      if (!ignore_pressure)
        std::swap(p_prev, p_next);
    }
    else {
      std::cout << "Previous: Timestep = " << sp.prev.t << ", filename = " << sp.prev.filename << std::endl;
      load_h5(folder + "/" + sp.prev.filename,
              ux_prev, uy_prev, uz_prev, rho_prev, p_prev,
              n[0], n[1], n[2],
              verbose, ignore_density, ignore_pressure, ignore_uz);
      // Extrapolate fields into solid
      if (boundary_mode == 0){
        fix_boundary(solid_ids_neigh, ux_prev, uy_prev, uz_prev, rho_prev, p_prev);
      }
    }

    std::cout << "Next: Timestep = " << sp.next.t << ", filename = " << sp.next.filename << std::endl;
    load_h5(folder + "/" + sp.next.filename,
            ux_next, uy_next, uz_next, rho_next, p_next,
            n[0], n[1], n[2], 
            verbose, ignore_density, ignore_pressure, ignore_uz);

    // Extrapolate fields into solid
    if (boundary_mode == 0){
      fix_boundary(solid_ids_neigh, ux_next, uy_next, uz_next, rho_next, p_next);
    }

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  // alpha_t = sp.weight_next(t);
  t_update = t;
}

void StructuredInterpol::probe(const Vector3d &x, const double t){
  assert(t <= t_next && t >= t_prev);
  alpha_t = (t-t_prev)/(t_next-t_prev);

  probe_space(x);
  probe_grad();

  // Not in use
  //for (Uint i=0; i<3; ++i){
  //  ind_pc[i] = imodulo(round(x[i]/dx[i]), n[i]);
  //}

  // Precompute velocities
  double Ux_prev = weighted_sum(ux_prev, ind, w);
  double Ux_next = weighted_sum(ux_next, ind, w);
  Ux = alpha_t * Ux_next + (1-alpha_t) * Ux_prev;
  Ax = (Ux_next-Ux_prev)/(t_next-t_prev);

  double Uy_prev = weighted_sum(uy_prev, ind, w);
  double Uy_next = weighted_sum(uy_next, ind, w);
  Uy = alpha_t * Uy_next + (1-alpha_t) * Uy_prev;
  Ay = (Uy_next-Uy_prev)/(t_next-t_prev);

  if (!ignore_uz){
    double Uz_prev = weighted_sum(uz_prev, ind, w);
    double Uz_next = weighted_sum(uz_next, ind, w);
    Uz = alpha_t * Uz_next + (1-alpha_t) * Uz_prev;
    Az = (Uz_next-Uz_prev)/(t_next-t_prev);
  }
}

bool StructuredInterpol::inside_domain(){
  return inside_domain_factor > 0.0;
}

void StructuredInterpol::probe_space(const Vector3d &x){
  // Computes ind, wq, w and inside_domain_factor
  int ix_fl[3];
  for (Uint i=0; i<3; ++i){
    ix_fl[i] = floor(x[i]/dx[i]);
  }
  for (Uint i=0; i<3; ++i){
    ind[i][0] = imodulo(ix_fl[i], n[i]);
    ind[i][1] = imodulo(ind[i][0] + 1, n[i]);
  }
  for (Uint i=0; i<3; ++i){
    double wxi = (x[i]-dx[i]*ix_fl[i])/dx[i];
    wq[i][0] = 1-wxi;
    wq[i][1] =   wxi;
  }
  for (Uint i=0; i<2; ++i){
    for (Uint j=0; j<2; ++j){
      for (Uint k=0; k<2; ++k){
        w[i][j][k] = wq[0][i]*wq[1][j]*wq[2][k];
      }
    }
  }
  if (boundary_mode == 0){
    inside_domain_factor = std::max(weighted_sum(levelZ, ind, w), 0.0);
  }
  else if (boundary_mode == 1) {
    inside_domain_factor = weighted_sum(levelZ, ind, w) > 0.0 ? 1.0 : 0.0;
  }
  else { // boundary_mode == 2
    inside_domain_factor = weighted_sum(levelZ, ind, w) < 1.0 ? 0.0 : 1.0;
  }
}

void StructuredInterpol::probe_grad(){
  for (Uint i=0; i<2; ++i){
    for (Uint j=0; j<2; ++j){
      for (Uint k=0; k<2; ++k){
        dw_x[i][j][k] = dwq[0][i]*wq[1][j]*wq[2][k];
        dw_y[i][j][k] = wq[0][i]*dwq[1][j]*wq[2][k];
        dw_z[i][j][k] = wq[0][i]*wq[1][j]*dwq[2][k];
      }
    }
  }

  // Factor used to determine whether point is in the domain or not
  // and to strictly enforce no-slip boundary conditions
  if (boundary_mode == 0){
    //inside_domain_factor = std::max(weighted_sum(levelZ, ind, w), 0.0);
    inside_domain_factor_x = weighted_sum(levelZ, ind, dw_x);
    inside_domain_factor_y = weighted_sum(levelZ, ind, dw_y);
    inside_domain_factor_z = weighted_sum(levelZ, ind, dw_z);
  }
}

// Interpolate in space and time and enforce BCs
double StructuredInterpol::get_ux(){
  return Ux*inside_domain_factor;
}

double StructuredInterpol::get_uy(){
  return Uy*inside_domain_factor;
}

double StructuredInterpol::get_uz(){
  return Uz*inside_domain_factor;
}

double StructuredInterpol::get_ax(){
  return Ax*inside_domain_factor;
}

double StructuredInterpol::get_ay(){
  return Ay*inside_domain_factor;
}

double StructuredInterpol::get_az(){
  return Az*inside_domain_factor;
}

double StructuredInterpol::get_rho(){
  double Rho_prev = weighted_sum(rho_prev, ind, w);
  double Rho_next = weighted_sum(rho_next, ind, w);
  //double Rho_prev = rho_prev[ind_pc[0]][ind_pc[1]][ind_pc[2]];
  //double Rho_next = rho_next[ind_pc[0]][ind_pc[1]][ind_pc[2]];
  return alpha_t * Rho_next + (1-alpha_t) * Rho_prev;
}

double StructuredInterpol::get_p(){
  double P_prev = weighted_sum(p_prev, ind, w);
  double P_next = weighted_sum(p_next, ind, w);
  return alpha_t * P_next + (1-alpha_t) * P_prev;
}

double StructuredInterpol::get_uxx(){
  double Uxx_prev = weighted_sum(ux_prev, ind, dw_x);
  double Uxx_next = weighted_sum(ux_next, ind, dw_x);
  double Uxx = alpha_t * Uxx_next + (1-alpha_t) * Uxx_prev;
  return Ux * inside_domain_factor_x + Uxx * inside_domain_factor;
}

double StructuredInterpol::get_uxy(){
  double Uxy_prev = weighted_sum(ux_prev, ind, dw_y);
  double Uxy_next = weighted_sum(ux_next, ind, dw_y);
  double Uxy = alpha_t * Uxy_next + (1-alpha_t) * Uxy_prev;
  return Ux * inside_domain_factor_y + Uxy * inside_domain_factor;
}

double StructuredInterpol::get_uxz(){
  double Uxz_prev = weighted_sum(ux_prev, ind, dw_z);
  double Uxz_next = weighted_sum(ux_next, ind, dw_z);
  double Uxz = alpha_t * Uxz_next + (1-alpha_t) * Uxz_prev;
  return Ux * inside_domain_factor_z + Uxz * inside_domain_factor;
}

double StructuredInterpol::get_uyx(){
  double Uyx_prev = weighted_sum(uy_prev, ind, dw_x);
  double Uyx_next = weighted_sum(uy_next, ind, dw_x);
  double Uyx = alpha_t * Uyx_next + (1-alpha_t) * Uyx_prev;
  return Uy * inside_domain_factor_x + Uyx * inside_domain_factor;
}

double StructuredInterpol::get_uyy(){
  double Uyy_prev = weighted_sum(uy_prev, ind, dw_y);
  double Uyy_next = weighted_sum(uy_next, ind, dw_y);
  double Uyy = alpha_t * Uyy_next + (1-alpha_t) * Uyy_prev;
  return Uy * inside_domain_factor_y + Uyy * inside_domain_factor;
}

double StructuredInterpol::get_uyz(){
  double Uyz_prev = weighted_sum(uy_prev, ind, dw_z);
  double Uyz_next = weighted_sum(uy_next, ind, dw_z);
  double Uyz = alpha_t * Uyz_next + (1-alpha_t) * Uyz_prev;
  return Uy * inside_domain_factor_z + Uyz * inside_domain_factor;
}

double StructuredInterpol::get_uzx(){
  double Uzx_prev = weighted_sum(uz_prev, ind, dw_x);
  double Uzx_next = weighted_sum(uz_next, ind, dw_x);
  double Uzx = alpha_t * Uzx_next + (1-alpha_t) * Uzx_prev;
  return Uz * inside_domain_factor_x + Uzx * inside_domain_factor;
}

double StructuredInterpol::get_uzy(){
  double Uzy_prev = weighted_sum(uz_prev, ind, dw_y);
  double Uzy_next = weighted_sum(uz_next, ind, dw_y);
  double Uzy = alpha_t * Uzy_next + (1-alpha_t) * Uzy_prev;
  return Uz * inside_domain_factor_y + Uzy * inside_domain_factor;
}

double StructuredInterpol::get_uzz(){
  double Uzz_prev = weighted_sum(uz_prev, ind, dw_z);
  double Uzz_next = weighted_sum(uz_next, ind, dw_z);
  double Uzz = alpha_t * Uzz_next + (1-alpha_t) * Uzz_prev;
  return Uz * inside_domain_factor_z + Uzz * inside_domain_factor;
}

Matrix3d StructuredInterpol::get_grada(){
  double dt_inv = 1./(t_next-t_prev);
  double Axx = dt_inv*(weighted_sum(ux_next, ind, dw_x)-weighted_sum(ux_prev, ind, dw_x));
  double Axy = dt_inv*(weighted_sum(ux_next, ind, dw_y)-weighted_sum(ux_prev, ind, dw_y));
  double Axz = dt_inv*(weighted_sum(ux_next, ind, dw_z)-weighted_sum(ux_prev, ind, dw_z));
  double Ayx = dt_inv*(weighted_sum(uy_next, ind, dw_x)-weighted_sum(uy_prev, ind, dw_x));
  double Ayy = dt_inv*(weighted_sum(uy_next, ind, dw_y)-weighted_sum(uy_prev, ind, dw_y));
  double Ayz = dt_inv*(weighted_sum(uy_next, ind, dw_z)-weighted_sum(uy_prev, ind, dw_z));
  double Azx = dt_inv*(weighted_sum(uz_next, ind, dw_x)-weighted_sum(uz_prev, ind, dw_x));
  double Azy = dt_inv*(weighted_sum(uz_next, ind, dw_y)-weighted_sum(uz_prev, ind, dw_y));
  double Azz = dt_inv*(weighted_sum(uz_next, ind, dw_z)-weighted_sum(uz_prev, ind, dw_z));
  Matrix3d M;
  M << Axx, Axy, Axz,
    Ayx, Ayy, Ayz,
    Azx, Azy, Azz;
  return M;
}

#endif
