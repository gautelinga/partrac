#ifndef __STRUCTUREDINTERPOL_HPP
#define __STRUCTUREDINTERPOL_HPP

#include "Interpol.hpp"
#include "Timestamps.hpp"
#include "H5Cpp.h"

void compute_ind_pc(Uint* ind_pc, const Vector3d &x, const Vector3d& dx, const Uint n[3]){
  // Constant
  for (Uint i=0; i<3; ++i){
    ind_pc[i] = imodulo(round(x[i]/dx[i]), n[i]);
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
  dset.read(Uv.data(), H5::PredType::NATIVE_DOUBLE, dspace, dspace);
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
  dset.read(Uv.data(), H5::PredType::NATIVE_INT, dspace, dspace);
  for (int ix=0; ix<nx; ++ix){
    for (int iy=0; iy<ny; ++iy){
      for (int iz=0; iz<nz; ++iz){
  	    u[ix][iy][iz] = Uv[nx*ny*iz+nx*iy+ix];
      }
    }
  }
}

static void load_int_field_as_bool(H5::H5File &h5file
                                 , bool*** u
                                 , const std::string field
                                 , const int nx
                                 , const int ny
                                 , const int nz
)
{
  H5::DataSet dset = h5file.openDataSet(field);
  H5::DataSpace dspace = dset.getSpace();
  std::vector<int> Uv(nx*ny*nz);
  dset.read(Uv.data(), H5::PredType::NATIVE_INT, dspace, dspace);
  for (int ix=0; ix<nx; ++ix){
    for (int iy=0; iy<ny; ++iy){
      for (int iz=0; iz<nz; ++iz){
  	    u[ix][iy][iz] = Uv[nx*ny*iz+nx*iy+ix] != 0;
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

static double inner_product(const double ww[2][2][2], const double bb[2][2][2]){
  double sum = 0.;
  for (Uint i=0; i<2; ++i){
    for (Uint j=0; j<2; ++j){
      for (Uint k=0; k<2; ++k){
        sum += ww[i][j][k] * bb[i][j][k];
      }
    }
  }
  return sum;
}
static void matrix_product(double v[3][3][3], double*** u, const Uint ind[3][2], const double W[3][3][3][2][2][2]){
  for (Uint i=0; i<3; ++i){
    for (Uint j=0; j<3; ++j){
      for (Uint k=0; k<3; ++k){
        v[i][j][k] = 0.;
        for (Uint l=0; l<2; ++l){
          for (Uint m=0; m<2; ++m){
            for (Uint n=0; n<2; ++n){
              v[i][j][k] += W[i][j][k][l][m][n] * u[ind[0][l]][ind[1][m]][ind[2][n]];
            }
          }
        }
      }
    }
  }
}

static void enforce_noslip(double v[3][3][3], bool*** isSolid, const Uint ind[3][2]){
  for (Uint i=0; i<2; ++i){
    for (Uint j=0; j<2; ++j){
      for (Uint k=0; k<2; ++k){
        if (isSolid[ind[0][i]][ind[1][j]][ind[2][k]]){
          for (Uint l=0; l<2; ++l){
            for (Uint m=0; m<2; ++m){
              for (Uint n=0; n<2; ++n){
                v[i+l][j+m][k+n] = 0.;
              }
            }
          }
        }
      }
    }
  }
}

static void compute_solid_local(bool is_solid_3[3][3][3], bool*** isSolid, const Uint ind[3][2]){
  for (Uint i=0; i<3; ++i){
    for (Uint j=0; j<3; ++j){
      for (Uint k=0; k<3; ++k){
        is_solid_3[i][j][k] = false;
      }
    }
  }
  for (Uint i=0; i<2; ++i){
    for (Uint j=0; j<2; ++j){
      for (Uint k=0; k<2; ++k){
        if (isSolid[ind[0][i]][ind[1][j]][ind[2][k]]){
          for (Uint l=0; l<2; ++l){
            for (Uint m=0; m<2; ++m){
              for (Uint n=0; n<2; ++n){
                is_solid_3[i+l][j+m][k+n] = true;
              }
            }
          }
        }
      }
    }
  }
}

static void enforce_noslip(double V[2][2][2], const bool is_solid_2[2][2][2]){
  for (Uint i=0; i<2; ++i){
    for (Uint j=0; j<2; ++j){
      for (Uint k=0; k<2; ++k){
        if (is_solid_2[i][j][k])
          V[i][j][k] = 0.;
      }
    }
  }
}

template<typename T>
static void get_subcube(T V[2][2][2], const T v[3][3][3], const bool sub_x[3]){
  for (Uint i=0; i<2; ++i){
    for (Uint j=0; j<2; ++j){
      for (Uint k=0; k<2; ++k){
        V[i][j][k] = v[ sub_x[0] + i][ sub_x[1] + j ][ sub_x[2] + k ];
      }
    }
  }
}

//static void compute_velocity_subcube(double V[2][2][2], double*** u, bool*** isSolid, const bool sub_x[3], const Uint ind[3][2], const double W[3][3][3][2][2][2]){
static void compute_velocity_subcube(double V[2][2][2], double*** u, const bool is_solid_2[2][2][2], const bool sub_x[3], const Uint ind[3][2], const double W[3][3][3][2][2][2]){
  double v[3][3][3];
  matrix_product(v, u, ind, W);
  get_subcube(V, v, sub_x);
  enforce_noslip(V, is_solid_2);
}

class StructuredInterpol : public Interpol {
public:
  StructuredInterpol(const std::string& infilename);
  void update(const double t);
  void probe(const Vector3d &x, const double t);
  void probe(const Vector3d &x, const double t, int& cell_id) { probe(x, t); };
  bool probe_light(const Vector3d &x, const double t, int& cell_id);
  void probe_heavy(const Vector3d &x, const double t, const int cell_id, PointValues& fields);
  bool compute_ind(const Vector3d &x, Uint _ind[3][2], int _ix_fl[3]);
  void probe_space_bulk(const Vector3d &x, 
    const Uint _ind[3][2],
    const int _ix_fl[3],
    double _w[2][2][2],
    double _dw_x[2][2][2],
    double _dw_y[2][2][2],
    double _dw_z[2][2][2]);
  void probe_space_boundary(
    const Vector3d &x, 
    const Uint _ind[3][2],
    const int _ix_fl[3],
    bool _is_solid_2[2][2][2],
    bool _sub_x[3],
    double _wux[2][2][2],
    double _wuy[2][2][2],
    double _wuz[2][2][2],
    double _dwux_x[2][2][2],
    double _dwux_y[2][2][2],
    double _dwux_z[2][2][2],
    double _dwuy_x[2][2][2],
    double _dwuy_y[2][2][2],
    double _dwuy_z[2][2][2],
    double _dwuz_x[2][2][2],
    double _dwuz_y[2][2][2],
    double _dwuz_z[2][2][2]);
  bool inside_domain() const;
  //bool inside_domain(const Vector3d &x) const;
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
  void probe(const Vector3d &x){ probe(x, this->t_update); };
protected:
  void probe_space(const Vector3d &x);
  Timestamps ts;
  double t_prev = 0.;
  double t_next = 0.;
  double alpha_t;

  Uint n[3] = {0, 0, 0};
  Vector3d dx;

  double wq[3][2];
  double dwq[3][2];

  bool*** isSolid;
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
  Uint ind_pc[3] = {0, 0, 0};  // piecewise constant intp

  double w[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dw_x[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dw_y[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dw_z[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};

  double wux[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dwux_x[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dwux_y[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dwux_z[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};

  double wuy[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dwuy_x[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dwuy_y[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dwuy_z[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};

  double wuz[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dwuz_x[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dwuz_y[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dwuz_z[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};

  double W[3][3][3][2][2][2];

  double Vx_prev[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double Vy_prev[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double Vz_prev[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double Vx_next[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double Vy_next[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double Vz_next[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};

  std::map<std::string, std::string> felbm_params;

  bool is_bulk = false;
  bool is_inside_domain = false;

  double Ux = 0.;
  double Uy = 0.;
  double Uz = 0.;
  double Ax = 0.;
  double Ay = 0.;
  double Az = 0.;

  bool ignore_density = false;
  bool ignore_pressure = false;
  bool ignore_uz = false;
};

StructuredInterpol::StructuredInterpol(const std::string& infilename) : Interpol(infilename) {
  std::ifstream input(infilename);
  if (!input){
    std::cout << "File " << infilename <<" doesn't exist." << std::endl;
    exit(0);
  }
  // Default params
  felbm_params["timestamps"] = "timestamps.dat";
  felbm_params["is_solid_file"] = "output_is_solid.h5";

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
    //dwwq[i][0] = -2./dx[i];
    //dwwq[i][1] =  2./dx[i];
  }

  // Create arrays
  isSolid = new bool**[n[0]];
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
    isSolid[ix] = new bool*[n[1]];
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
      isSolid[ix][iy] = new bool[n[2]];
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

  load_int_field_as_bool(solid_file, isSolid, "is_solid", n[0], n[1], n[2]);

  double wwx[3][2];
  for (Uint i=0; i<3; ++i){
    for (Uint j=0; j<3; ++j){
      for (Uint k=0; k<3; ++k){            
        wwx[0][0] = 1. - (double(i))/2;
        wwx[1][0] = 1. - (double(j))/2;
        wwx[2][0] = 1. - (double(k))/2;
        for (Uint d=0; d<3; ++d){
          wwx[d][1] = 1. - wwx[d][0];
        }
        for (Uint l=0; l<2; ++l){
          for (Uint m=0; m<2; ++m){
            for (Uint n=0; n<2; ++n){
              W[i][j][k][l][m][n] = wwx[0][l] * wwx[1][m] * wwx[2][n];
            }
          }
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
    }

    std::cout << "Next: Timestep = " << sp.next.t << ", filename = " << sp.next.filename << std::endl;
    load_h5(folder + "/" + sp.next.filename,
            ux_next, uy_next, uz_next, rho_next, p_next,
            n[0], n[1], n[2], 
            verbose, ignore_density, ignore_pressure, ignore_uz);

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

  double Ux_prev, Uy_prev, Uz_prev;
  double Ux_next, Uy_next, Uz_next;

  // Precompute velocities
  if (is_bulk)
  {
    Ux_prev = weighted_sum(ux_prev, ind, w);
    Ux_next = weighted_sum(ux_next, ind, w);

    Uy_prev = weighted_sum(uy_prev, ind, w);
    Uy_next = weighted_sum(uy_next, ind, w);

    Uz_prev = weighted_sum(uz_prev, ind, w);
    Uz_next = weighted_sum(uz_next, ind, w);
  }
  else
  {
    Ux_prev = inner_product(wux, Vx_prev);
    Ux_next = inner_product(wux, Vx_next);

    Uy_prev = inner_product(wuy, Vy_prev);
    Uy_next = inner_product(wuy, Vy_next);

    Uz_prev = inner_product(wuz, Vz_prev);
    Uz_next = inner_product(wuz, Vz_next);
  }

  Ux = alpha_t * Ux_next + (1-alpha_t) * Ux_prev;
  Ax = (Ux_next-Ux_prev)/(t_next-t_prev);

  Uy = alpha_t * Uy_next + (1-alpha_t) * Uy_prev;
  Ay = (Uy_next-Uy_prev)/(t_next-t_prev);

  Uz = alpha_t * Uz_next + (1-alpha_t) * Uz_prev;
  Az = (Uz_next-Uz_prev)/(t_next-t_prev);
}

bool StructuredInterpol::probe_light(const Vector3d &x, const double t, int& cell_id){
  Uint _ind_pc[3];
  compute_ind_pc(_ind_pc, x, dx, n);
  return !isSolid[_ind_pc[0]][_ind_pc[1]][_ind_pc[2]];
}

void StructuredInterpol::probe_heavy(const Vector3d &x, const double t, const int cell_id, PointValues& fields){
  // Assuming probe_light has already been called and found that the cell is not in solid
  double Ux_prev, Uy_prev, Uz_prev;
  double Ux_next, Uy_next, Uz_next;
  double Rho_prev, Rho_next;
  double P_prev, P_next;

  double Uxx_prev, Uxx_next, Uxy_prev, Uxy_next, Uxz_prev, Uxz_next;
  double Uyx_prev, Uyx_next, Uyy_prev, Uyy_next, Uyz_prev, Uyz_next;
  double Uzx_prev, Uzx_next, Uzy_prev, Uzy_next, Uzz_prev, Uzz_next;

  Uint _ind[3][2] = {{0, 0}, {0, 0}, {0, 0}};
  int _ix_fl[3];
  bool _is_bulk = compute_ind(x, _ind, _ix_fl);

  // Precompute velocities
  if (_is_bulk) // Bulk cell
  {
    //_w, _dw_x, _dw_y, _dw_z
    double _w[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _dw_x[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _dw_y[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _dw_z[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};

    probe_space_bulk(x, _ind, _ix_fl, _w, _dw_x, _dw_y, _dw_z);

    Ux_prev = weighted_sum(ux_prev, _ind, _w);
    Ux_next = weighted_sum(ux_next, _ind, _w);

    Uy_prev = weighted_sum(uy_prev, _ind, _w);
    Uy_next = weighted_sum(uy_next, _ind, _w);

    Uz_prev = weighted_sum(uz_prev, _ind, _w);
    Uz_next = weighted_sum(uz_next, _ind, _w);

    Rho_prev = weighted_sum(rho_prev, _ind, _w);
    Rho_next = weighted_sum(rho_next, _ind, _w);

    P_prev = weighted_sum(p_prev, _ind, _w);
    P_next = weighted_sum(p_next, _ind, _w);

    Uxx_prev = weighted_sum(ux_prev, _ind, _dw_x);
    Uxx_next = weighted_sum(ux_next, _ind, _dw_x);
    Uxy_prev = weighted_sum(ux_prev, _ind, _dw_y);
    Uxy_next = weighted_sum(ux_next, _ind, _dw_y);
    Uxz_prev = weighted_sum(ux_prev, _ind, _dw_z);
    Uxz_next = weighted_sum(ux_next, _ind, _dw_z);
    Uyx_prev = weighted_sum(uy_prev, _ind, _dw_x);
    Uyx_next = weighted_sum(uy_next, _ind, _dw_x);
    Uyy_prev = weighted_sum(uy_prev, _ind, _dw_y);
    Uyy_next = weighted_sum(uy_next, _ind, _dw_y);
    Uyz_prev = weighted_sum(uy_prev, _ind, _dw_z);
    Uyz_next = weighted_sum(uy_next, _ind, _dw_z);
    Uzx_prev = weighted_sum(uz_prev, _ind, _dw_x);
    Uzx_next = weighted_sum(uz_next, _ind, _dw_x);
    Uzy_prev = weighted_sum(uz_prev, _ind, _dw_y);
    Uzy_next = weighted_sum(uz_next, _ind, _dw_y);
    Uzz_prev = weighted_sum(uz_prev, _ind, _dw_z);
    Uzz_next = weighted_sum(uz_next, _ind, _dw_z);
  }
  else // Close to boundary
  {
    bool _is_solid_2[2][2][2];
    bool _sub_x[3];

    double _wux[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _dwux_x[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _dwux_y[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _dwux_z[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};

    double _wuy[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _dwuy_x[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _dwuy_y[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _dwuy_z[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};

    double _wuz[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _dwuz_x[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _dwuz_y[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _dwuz_z[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};

    probe_space_boundary(x, _ind, _ix_fl, _is_solid_2, _sub_x,
                         _wux, _wuy, _wuz,
                         _dwux_x, _dwux_y, _dwux_z, 
                         _dwuy_x, _dwuy_y, _dwuy_z,
                         _dwuz_x, _dwuz_y, _dwuz_z);


    double _Vx_prev[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _Vy_prev[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _Vz_prev[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _Vx_next[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _Vy_next[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
    double _Vz_next[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};

    compute_velocity_subcube(_Vx_prev, ux_prev, _is_solid_2, _sub_x, _ind, W);
    compute_velocity_subcube(_Vy_prev, uy_prev, _is_solid_2, _sub_x, _ind, W);
    compute_velocity_subcube(_Vz_prev, uz_prev, _is_solid_2, _sub_x, _ind, W);

    compute_velocity_subcube(_Vx_next, ux_next, _is_solid_2, _sub_x, _ind, W);
    compute_velocity_subcube(_Vy_next, uy_next, _is_solid_2, _sub_x, _ind, W);
    compute_velocity_subcube(_Vz_next, uz_next, _is_solid_2, _sub_x, _ind, W);

    Ux_prev = inner_product(_wux, _Vx_prev);
    Ux_next = inner_product(_wux, _Vx_next);

    Uy_prev = inner_product(_wuy, _Vy_prev);
    Uy_next = inner_product(_wuy, _Vy_next);

    Uz_prev = inner_product(_wuz, _Vz_prev);
    Uz_next = inner_product(_wuz, _Vz_next);

    Uint _ind_pc[3];
    compute_ind_pc(_ind_pc, x, dx, n);
    Rho_prev = rho_prev[_ind_pc[0]][_ind_pc[1]][_ind_pc[2]];
    Rho_next = rho_next[_ind_pc[0]][_ind_pc[1]][_ind_pc[2]];

    P_prev = p_prev[_ind_pc[0]][_ind_pc[1]][_ind_pc[2]];
    P_next = p_next[_ind_pc[0]][_ind_pc[1]][_ind_pc[2]];

    Uxx_prev = inner_product(_dwux_x, _Vx_prev);
    Uxx_next = inner_product(_dwux_x, _Vx_prev);
    Uxy_prev = inner_product(_dwux_y, _Vx_prev);
    Uxy_next = inner_product(_dwux_y, _Vx_next);
    Uxz_prev = inner_product(_dwux_z, _Vx_prev);
    Uxz_next = inner_product(_dwux_z, _Vx_next);
    Uyx_prev = inner_product(_dwux_x, _Vy_prev);
    Uyx_next = inner_product(_dwux_x, _Vy_prev);
    Uyy_prev = inner_product(_dwux_y, _Vy_prev);
    Uyy_next = inner_product(_dwux_y, _Vy_next);
    Uyz_prev = inner_product(_dwux_z, _Vy_prev);
    Uyz_next = inner_product(_dwux_z, _Vy_next);
    Uzx_prev = inner_product(_dwux_x, _Vz_prev);
    Uzx_next = inner_product(_dwux_x, _Vz_prev);
    Uzy_prev = inner_product(_dwux_y, _Vz_prev);
    Uzy_next = inner_product(_dwux_y, _Vz_next);
    Uzz_prev = inner_product(_dwux_z, _Vz_prev);
    Uzz_next = inner_product(_dwux_z, _Vz_next);
  }

  fields.U = { alpha_t * Ux_next + (1-alpha_t) * Ux_prev,
               alpha_t * Uy_next + (1-alpha_t) * Uy_prev,
               alpha_t * Uz_next + (1-alpha_t) * Uz_prev };
  fields.A = { (Ux_next-Ux_prev)/(t_next-t_prev),
               (Uy_next-Uy_prev)/(t_next-t_prev),
               (Uz_next-Uz_prev)/(t_next-t_prev) };

  fields.P = alpha_t * P_next + (1-alpha_t) * P_prev;
  fields.Rho = alpha_t * Rho_next + (1-alpha_t) * Rho_prev;

  Matrix3d gradU_prev;
  gradU_prev << Uxx_prev, Uxy_prev, Uxz_prev,
                Uyx_prev, Uyy_prev, Uyz_prev,
                Uzx_prev, Uzy_prev, Uzz_prev;

  Matrix3d gradU_next;
  gradU_next << Uxx_next, Uxy_next, Uxz_next,
                Uyx_next, Uyy_next, Uyz_next,
                Uzx_next, Uzy_next, Uzz_next;

  fields.gradU = alpha_t * gradU_next + (1-alpha_t) * gradU_prev;

}

bool StructuredInterpol::inside_domain() const {
  return is_inside_domain;
}

bool StructuredInterpol::compute_ind(const Vector3d &x, Uint _ind[3][2], int _ix_fl[3]){
  // Assuming this cell is not inside the solid phase
  for (Uint i=0; i<3; ++i){
    _ix_fl[i] = floor(x[i]/dx[i]);
  }

  //Uint _ind[3][2] = {{0, 0}, {0, 0}, {0, 0}};  // trilinear intp
  for (Uint i=0; i<3; ++i){
    _ind[i][0] = imodulo(_ix_fl[i], n[i]);
    _ind[i][1] = imodulo(_ind[i][0] + 1, n[i]);
  }

  for (Uint i=0; i<2; ++i){
    for (Uint j=0; j<2; ++j){
      for (Uint k=0; k<2; ++k){
        if (isSolid[_ind[0][i]][_ind[1][j]][_ind[2][k]]) {
          return false;
        }
      }
    }
  }
  return true;
}

void StructuredInterpol::probe_space_bulk(const Vector3d &x, 
    const Uint _ind[3][2],
    const int _ix_fl[3],
    double _w[2][2][2],
    double _dw_x[2][2][2],
    double _dw_y[2][2][2],
    double _dw_z[2][2][2]
  ){
  // Computes ind, w, dw...
  double _wq[3][2];

  for (Uint i=0; i<3; ++i){
    double wxi = (x[i]-dx[i]*_ix_fl[i])/dx[i];
    _wq[i][0] = 1 - wxi;
    _wq[i][1] =     wxi;
  }

  for (Uint i=0; i<2; ++i){
    for (Uint j=0; j<2; ++j){
      for (Uint k=0; k<2; ++k){
        _w[i][j][k] = _wq[0][i] * _wq[1][j] * _wq[2][k];
      }
    }
  }
  for (Uint i=0; i<2; ++i){
    for (Uint j=0; j<2; ++j){
      for (Uint k=0; k<2; ++k){
        _dw_x[i][j][k] = dwq[0][i] * _wq[1][j] * _wq[2][k];
        _dw_y[i][j][k] = _wq[0][i] * dwq[1][j] * _wq[2][k];
        _dw_z[i][j][k] = _wq[0][i] * _wq[1][j] * dwq[2][k];
      }
    }
  }
}

void StructuredInterpol::probe_space_boundary(
  const Vector3d &x, 
  const Uint _ind[3][2],
  const int _ix_fl[3],
  bool _is_solid_2[2][2][2],
  bool _sub_x[3],
  double _wux[2][2][2],
  double _wuy[2][2][2],
  double _wuz[2][2][2],
  double _dwux_x[2][2][2],
  double _dwux_y[2][2][2],
  double _dwux_z[2][2][2],
  double _dwuy_x[2][2][2],
  double _dwuy_y[2][2][2],
  double _dwuy_z[2][2][2],
  double _dwuz_x[2][2][2],
  double _dwuz_y[2][2][2],
  double _dwuz_z[2][2][2]
  )
{
  double xd[3];
  for (Uint i=0; i<3; ++i){
    xd[i] = x[i]/dx[i] - _ix_fl[i];
  }

  //bool sub_x[3];
  for (Uint i=0; i<3; ++i){
    _sub_x[i] = xd[i] >= 0.5;
  }

  bool is_solid_3[3][3][3];
  //bool is_solid_2[2][2][2];
  compute_solid_local(is_solid_3, isSolid, _ind);
  get_subcube(_is_solid_2, is_solid_3, _sub_x);

  double _wq[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};;
  for (Uint i=0; i<3; ++i){
    double wxi = _sub_x[i] ? 2 * xd[i] - 1.0: 2 * xd[i];
    _wq[i][0] = 1 - wxi;
    _wq[i][1] =     wxi;
  }

  double gamma = 2;
  for (Uint i=0; i<2; ++i){
    for (Uint j=0; j<2; ++j){
      for (Uint k=0; k<2; ++k){
        double wqux = wq[0][i];
        double dwqux = dwq[0][i];
        double wquy = wq[1][j];
        double dwquy = dwq[1][j];
        double wquz = wq[2][k];
        double dwquz = dwq[2][k];

        if (_is_solid_2[i == 0 ? 1 : 0][j][k]){
          wqux = pow(_wq[0][i], gamma);
          dwqux = gamma * pow(_wq[0][i], gamma-1) * dwq[0][i];
        }

        if (_is_solid_2[i][j == 0 ? 1 : 0][k]){
          wquy = pow(_wq[1][j], gamma);
          dwquy = gamma * pow(_wq[1][j], gamma-1) * dwq[1][j];
        }

        if (_is_solid_2[i][j][k == 0 ? 1 : 0]){
          wquz = pow(_wq[2][k], gamma);
          dwquz = gamma * pow(_wq[2][k], gamma-1) * dwq[2][k];
        }

        _wux[i][j][k] = wqux     * _wq[1][j] * _wq[2][k];
        _wuy[i][j][k] = _wq[0][i] * wquy     * _wq[2][k];
        _wuz[i][j][k] = _wq[0][i] * _wq[1][j] * wquz;

        _dwux_x[i][j][k] = 2 * dwqux * _wq[1][j] * _wq[2][k];
        _dwux_y[i][j][k] = 2 *  wqux * dwq[1][j] * _wq[2][k];
        _dwux_z[i][j][k] = 2 *  wqux * _wq[1][j] * dwq[2][k];

        _dwuy_x[i][j][k] = 2 * dwq[0][i] *  wquy * _wq[2][k];
        _dwuy_y[i][j][k] = 2 * _wq[0][i] * dwquy * _wq[2][k];
        _dwuy_z[i][j][k] = 2 * _wq[0][i] *  wquy * dwq[2][k];

        _dwuz_x[i][j][k] = 2 * dwq[0][i] * _wq[1][j] *  wquz;
        _dwuz_y[i][j][k] = 2 * _wq[0][i] * dwq[1][j] *  wquz;
        _dwuz_z[i][j][k] = 2 * _wq[0][i] * _wq[1][j] * dwquz;
      }
    }
  }
}

void StructuredInterpol::probe_space(const Vector3d &x){
  // Computes ind, wq, w and inside_domain_factor

  // Constant
  for (Uint i=0; i<3; ++i){
    ind_pc[i] = imodulo(round(x[i]/dx[i]), n[i]);
  }
  is_inside_domain = !isSolid[ind_pc[0]][ind_pc[1]][ind_pc[2]];

  if (!is_inside_domain){
    // All weights to zero
    for (Uint i=0; i<2; ++i){
      for (Uint j=0; j<2; ++j){
        for (Uint k=0; k<2; ++k){
          w[i][j][k] = 0.;
          dw_x[i][j][k] = 0.;
          dw_y[i][j][k] = 0.;
          dw_z[i][j][k] = 0.;
        }
      }
    }
  }
  else {
    int ix_fl[3];
    for (Uint i=0; i<3; ++i){
      ix_fl[i] = floor(x[i]/dx[i]);
    }

    for (Uint i=0; i<3; ++i){
      ind[i][0] = imodulo(ix_fl[i], n[i]);
      ind[i][1] = imodulo(ind[i][0] + 1, n[i]);
    }

    is_bulk = true;
    for (Uint i=0; i<2; ++i){
      for (Uint j=0; j<2; ++j){
        for (Uint k=0; k<2; ++k){
          is_bulk = is_bulk && !isSolid[ind[0][i]][ind[1][j]][ind[2][k]];
        }
      }
    }

    // Testing:
    //is_bulk = false;

    if (is_bulk){
      for (Uint i=0; i<3; ++i){
        double wxi = (x[i]-dx[i]*ix_fl[i])/dx[i];
        wq[i][0] = 1 - wxi;
        wq[i][1] =     wxi;
      }

      for (Uint i=0; i<2; ++i){
        for (Uint j=0; j<2; ++j){
          for (Uint k=0; k<2; ++k){
            w[i][j][k] = wq[0][i] * wq[1][j] * wq[2][k];
          }
        }
      }
      for (Uint i=0; i<2; ++i){
        for (Uint j=0; j<2; ++j){
          for (Uint k=0; k<2; ++k){
            dw_x[i][j][k] = dwq[0][i] *  wq[1][j] *  wq[2][k];
            dw_y[i][j][k] =  wq[0][i] * dwq[1][j] *  wq[2][k];
            dw_z[i][j][k] =  wq[0][i] *  wq[1][j] * dwq[2][k];
          }
        }
      }
    }
    else {
      double xd[3];
      for (Uint i=0; i<3; ++i){
        xd[i] = x[i]/dx[i] - ix_fl[i];
      }

      bool sub_x[3];
      for (Uint i=0; i<3; ++i){
        sub_x[i] = xd[i] >= 0.5;
      }

      bool is_solid_3[3][3][3];
      bool is_solid_2[2][2][2];
      compute_solid_local(is_solid_3, isSolid, ind);
      get_subcube(is_solid_2, is_solid_3, sub_x);

      compute_velocity_subcube(Vx_prev, ux_prev, is_solid_2, sub_x, ind, W);
      compute_velocity_subcube(Vy_prev, uy_prev, is_solid_2, sub_x, ind, W);
      compute_velocity_subcube(Vz_prev, uz_prev, is_solid_2, sub_x, ind, W);

      compute_velocity_subcube(Vx_next, ux_next, is_solid_2, sub_x, ind, W);
      compute_velocity_subcube(Vy_next, uy_next, is_solid_2, sub_x, ind, W);
      compute_velocity_subcube(Vz_next, uz_next, is_solid_2, sub_x, ind, W);

      for (Uint i=0; i<3; ++i){
        double wxi = sub_x[i] ? 2 * xd[i] - 1.0: 2 * xd[i];
        wq[i][0] = 1 - wxi;
        wq[i][1] =     wxi;
      }

      double gamma = 2;
      for (Uint i=0; i<2; ++i){
        for (Uint j=0; j<2; ++j){
          for (Uint k=0; k<2; ++k){
            double wqux = wq[0][i];
            double dwqux = dwq[0][i];
            if (is_solid_2[i == 0 ? 1 : 0][j][k]){
              wqux = pow(wq[0][i], gamma);
              dwqux = gamma * pow(wq[0][i], gamma-1) * dwq[0][i];
            }
            double wquy = wq[1][j];
            double dwquy = dwq[1][j];
            if (is_solid_2[i][j == 0 ? 1 : 0][k]){
              wquy = pow(wq[1][j], gamma);
              dwquy = gamma * pow(wq[1][j], gamma-1) * dwq[1][j];
            }
            double wquz = wq[2][k];
            double dwquz = dwq[2][k];
            if (is_solid_2[i][j][k == 0 ? 1 : 0]){
              wquz = pow(wq[2][k], gamma);
              dwquz = gamma * pow(wq[2][k], gamma-1) * dwq[2][k];
            }

            wux[i][j][k] = wqux     * wq[1][j] * wq[2][k];
            wuy[i][j][k] = wq[0][i] * wquy     * wq[2][k];
            wuz[i][j][k] = wq[0][i] * wq[1][j] * wquz;

            dwux_x[i][j][k] = 2 * dwqux *  wq[1][j] *  wq[2][k];
            dwux_y[i][j][k] = 2 *  wqux * dwq[1][j] *  wq[2][k];
            dwux_z[i][j][k] = 2 *  wqux *  wq[1][j] * dwq[2][k];

            dwuy_x[i][j][k] = 2 * dwq[0][i] *  wquy *  wq[2][k];
            dwuy_y[i][j][k] = 2 *  wq[0][i] * dwquy *  wq[2][k];
            dwuy_z[i][j][k] = 2 *  wq[0][i] *  wquy * dwq[2][k];

            dwuz_x[i][j][k] = 2 * dwq[0][i] *  wq[1][j] *  wquz;
            dwuz_y[i][j][k] = 2 *  wq[0][i] * dwq[1][j] *  wquz;
            dwuz_z[i][j][k] = 2 *  wq[0][i] *  wq[1][j] * dwquz;
          }
        }
      }
    }
  }
}

// Interpolate in space and time and enforce BCs
double StructuredInterpol::get_ux(){
  return Ux;
}

double StructuredInterpol::get_uy(){
  return Uy;
}

double StructuredInterpol::get_uz(){
  return Uz;
}

double StructuredInterpol::get_ax(){
  return Ax;
}

double StructuredInterpol::get_ay(){
  return Ay;
}

double StructuredInterpol::get_az(){
  return Az;
}

double StructuredInterpol::get_rho(){
  double Rho_prev, Rho_next;
  if (is_bulk){
    Rho_prev = weighted_sum(rho_prev, ind, w);
    Rho_next = weighted_sum(rho_next, ind, w);
  }
  else {
    Rho_prev = rho_prev[ind_pc[0]][ind_pc[1]][ind_pc[2]];
    Rho_next = rho_next[ind_pc[0]][ind_pc[1]][ind_pc[2]];
  }
  return alpha_t * Rho_next + (1-alpha_t) * Rho_prev;
}

double StructuredInterpol::get_p(){
  double P_prev, P_next;
  if (is_bulk){
    P_prev = weighted_sum(p_prev, ind, w);
    P_next = weighted_sum(p_next, ind, w);
  }
  else {
    P_prev = p_prev[ind_pc[0]][ind_pc[1]][ind_pc[2]];
    P_next = p_next[ind_pc[0]][ind_pc[1]][ind_pc[2]];
  }
  return alpha_t * P_next + (1-alpha_t) * P_prev;
}

double StructuredInterpol::get_uxx(){
  double Uxx_prev, Uxx_next;
  if (is_bulk){
    Uxx_prev = weighted_sum(ux_prev, ind, dw_x);
    Uxx_next = weighted_sum(ux_next, ind, dw_x);
  }
  else {
    Uxx_prev = inner_product(dwux_x, Vx_prev);
    Uxx_next = inner_product(dwux_x, Vx_prev);
  }
  double Uxx = alpha_t * Uxx_next + (1-alpha_t) * Uxx_prev;
  return Uxx;
}

double StructuredInterpol::get_uxy(){
  double Uxy_prev, Uxy_next;
  if (is_bulk){
    Uxy_prev = weighted_sum(ux_prev, ind, dw_y);
    Uxy_next = weighted_sum(ux_next, ind, dw_y);
  }
  else {
    Uxy_prev = inner_product(dwux_y, Vx_prev);
    Uxy_next = inner_product(dwux_y, Vx_next);
  }
  double Uxy = alpha_t * Uxy_next + (1-alpha_t) * Uxy_prev;
  return Uxy;
}

double StructuredInterpol::get_uxz(){
  double Uxz_prev, Uxz_next;
  if (is_bulk){
    Uxz_prev = weighted_sum(ux_prev, ind, dw_z);
    Uxz_next = weighted_sum(ux_next, ind, dw_z);
  }
  else {
    Uxz_prev = inner_product(dwux_z, Vx_prev);
    Uxz_next = inner_product(dwux_z, Vx_next);    
  }
  double Uxz = alpha_t * Uxz_next + (1-alpha_t) * Uxz_prev;
  return Uxz;
}

double StructuredInterpol::get_uyx(){
  double Uyx_prev, Uyx_next;
  if (is_bulk){
    Uyx_prev = weighted_sum(uy_prev, ind, dw_x);
    Uyx_next = weighted_sum(uy_next, ind, dw_x);
  }
  else {
    Uyx_prev = inner_product(dwuy_x, Vy_prev);
    Uyx_next = inner_product(dwuy_x, Vy_next);  
  }
  double Uyx = alpha_t * Uyx_next + (1-alpha_t) * Uyx_prev;
  return Uyx;
}

double StructuredInterpol::get_uyy(){
  double Uyy_prev, Uyy_next;
  if (is_bulk){
    Uyy_prev = weighted_sum(uy_prev, ind, dw_y);
    Uyy_next = weighted_sum(uy_next, ind, dw_y);
  }
  else {
    Uyy_prev = inner_product(dwuy_y, Vy_prev);
    Uyy_next = inner_product(dwuy_y, Vy_next);
  }
  double Uyy = alpha_t * Uyy_next + (1-alpha_t) * Uyy_prev;
  return Uyy;
}

double StructuredInterpol::get_uyz(){
  double Uyz_prev, Uyz_next;
  if (is_bulk){
    Uyz_prev = weighted_sum(uy_prev, ind, dw_z);
    Uyz_next = weighted_sum(uy_next, ind, dw_z);
  }
  else {
    Uyz_prev = inner_product(dwuy_z, Vy_prev);
    Uyz_next = inner_product(dwuy_z, Vy_next);
  }
  double Uyz = alpha_t * Uyz_next + (1-alpha_t) * Uyz_prev;
  return Uyz;
}

double StructuredInterpol::get_uzx(){
  double Uzx_prev, Uzx_next;
  if (is_bulk){
    Uzx_prev = weighted_sum(uz_prev, ind, dw_x);
    Uzx_next = weighted_sum(uz_next, ind, dw_x);
  }
  else {
    Uzx_prev = inner_product(dwuz_x, Vz_prev);
    Uzx_next = inner_product(dwuz_x, Vz_next);
  }
  double Uzx = alpha_t * Uzx_next + (1-alpha_t) * Uzx_prev;
  return Uzx;
}

double StructuredInterpol::get_uzy(){
  double Uzy_prev, Uzy_next;
  if (is_bulk){
    Uzy_prev = weighted_sum(uz_prev, ind, dw_y);
    Uzy_next = weighted_sum(uz_next, ind, dw_y);
  }
  else {
    Uzy_prev = inner_product(dwuz_y, Vz_prev);
    Uzy_next = inner_product(dwuz_y, Vz_next);
  }
  double Uzy = alpha_t * Uzy_next + (1-alpha_t) * Uzy_prev;
  return Uzy;
}

double StructuredInterpol::get_uzz(){
  double Uzz_prev, Uzz_next;
  if (is_bulk){
    Uzz_prev = weighted_sum(uz_prev, ind, dw_z);
    Uzz_next = weighted_sum(uz_next, ind, dw_z);
  }
  else {
    Uzz_prev = inner_product(dwuz_z, Vz_prev);
    Uzz_next = inner_product(dwuz_z, Vz_next);
  }
  double Uzz = alpha_t * Uzz_next + (1-alpha_t) * Uzz_prev;
  return Uzz;
}

Matrix3d StructuredInterpol::get_grada(){
  double dt_inv = 1./(t_next-t_prev);
  double Axx, Axy, Axz, Ayx, Ayy, Ayz, Azx, Azy, Azz;
  if (is_bulk){
    Axx = dt_inv*(weighted_sum(ux_next, ind, dw_x)-weighted_sum(ux_prev, ind, dw_x));
    Axy = dt_inv*(weighted_sum(ux_next, ind, dw_y)-weighted_sum(ux_prev, ind, dw_y));
    Axz = dt_inv*(weighted_sum(ux_next, ind, dw_z)-weighted_sum(ux_prev, ind, dw_z));
    Ayx = dt_inv*(weighted_sum(uy_next, ind, dw_x)-weighted_sum(uy_prev, ind, dw_x));
    Ayy = dt_inv*(weighted_sum(uy_next, ind, dw_y)-weighted_sum(uy_prev, ind, dw_y));
    Ayz = dt_inv*(weighted_sum(uy_next, ind, dw_z)-weighted_sum(uy_prev, ind, dw_z));
    Azx = dt_inv*(weighted_sum(uz_next, ind, dw_x)-weighted_sum(uz_prev, ind, dw_x));
    Azy = dt_inv*(weighted_sum(uz_next, ind, dw_y)-weighted_sum(uz_prev, ind, dw_y));
    Azz = dt_inv*(weighted_sum(uz_next, ind, dw_z)-weighted_sum(uz_prev, ind, dw_z));
  }
  else {
    Axx = dt_inv*(inner_product(dwux_x, Vx_next)-inner_product(dwux_x, Vx_prev));
    Axy = dt_inv*(inner_product(dwux_y, Vx_next)-inner_product(dwux_y, Vx_prev));
    Axz = dt_inv*(inner_product(dwux_z, Vx_next)-inner_product(dwux_z, Vx_prev));
    Ayx = dt_inv*(inner_product(dwuy_x, Vy_next)-inner_product(dwuy_x, Vy_prev));
    Ayy = dt_inv*(inner_product(dwuy_y, Vy_next)-inner_product(dwuy_y, Vy_prev));
    Ayz = dt_inv*(inner_product(dwuy_z, Vy_next)-inner_product(dwuy_z, Vy_prev));
    Azx = dt_inv*(inner_product(dwuz_x, Vz_next)-inner_product(dwuz_x, Vz_prev));
    Azy = dt_inv*(inner_product(dwuz_y, Vz_next)-inner_product(dwuz_y, Vz_prev));
    Azz = dt_inv*(inner_product(dwuz_z, Vz_next)-inner_product(dwuz_z, Vz_prev));
  }
  Matrix3d M;
  M << Axx, Axy, Axz,
       Ayx, Ayy, Ayz,
       Azx, Azy, Azz;
  return M;
}

#endif
