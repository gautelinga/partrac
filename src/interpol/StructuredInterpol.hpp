#ifndef __STRUCTUREDINTERPOL_HPP
#define __STRUCTUREDINTERPOL_HPP

#include <algorithm>
#include <array>
#include <limits>
#include "Error.hpp"
#include "Interpol.hpp"
#include "Params.hpp"
#include "loader_params.hpp"
#include "H5Cpp.h"
#include "files.hpp"
#include "geometry.hpp"
#include "strings.hpp"
#include "Timestamps.hpp"

// std::round to int for |v| < 2^31, without the libm call
inline int round_to_int(const double v){
  const int r = static_cast<int>(v);
  const double f = v - r;
  return r + (f >= 0.5) - (f <= -0.5);
}

inline void compute_ind_pc(Uint* ind_pc, const Vector3d &x, const Vector3d& dx, const Uint n[3]){
  // Constant
  for (Uint i=0; i<3; ++i){
    ind_pc[i] = imodulo(round_to_int(x[i]/dx[i]), n[i]);
  }
}

// Unit normal toward the solid axis neighbours of x's nearest node; zero if none
template<typename Solid>
inline Vector3d lattice_wall_normal(const Vector3d& x, const Vector3d& dx, const Uint n[3], const Solid& solid){
  int idx[3];
  for (Uint i=0; i<3; ++i)
    idx[i] = round_to_int(x[i]/dx[i]);
  Vector3d nrm = Vector3d::Zero();
  for (Uint a=0; a<3; ++a)
    for (const int s : {-1, 1}){
      int j[3] = {idx[0], idx[1], idx[2]};
      j[a] += s;
      if (solid(imodulo(j[0], n[0]), imodulo(j[1], n[1]), imodulo(j[2], n[2])))
        nrm[a] += s;
    }
  const double len = nrm.norm();
  return len > 0. ? Vector3d(nrm/len) : nrm;
}

// View of one GridBlock field; valid while the block is unresized
template<typename T>
struct Grid3 {
  T* p = nullptr;
  std::size_t sx = 0, sy = 0, sz = 0;
  T& operator()(const Uint i, const Uint j, const Uint k){ return p[i*sx + j*sy + k*sz]; }
  const T& operator()(const Uint i, const Uint j, const Uint k) const { return p[i*sx + j*sy + k*sz]; }
};

// All fields in one block, z fastest, a node's fields adjacent
template<typename T>
struct GridBlock {
  Uint ny = 0, nz = 0, nf = 0;   // nx not stored
  std::vector<T> v;
  void resize(const Uint nx, const Uint ny_, const Uint nz_, const Uint nf_){
    ny = ny_; nz = nz_; nf = nf_;
    v.assign(std::size_t(nx)*ny*nz*nf, T());
  }
  Grid3<T> field(const Uint f){
    return {v.data() + f, std::size_t(ny)*nz*nf, std::size_t(nz)*nf, nf};
  }
};

inline void load_field(H5::H5File &h5file
                     , Grid3<double>& u
                     , const std::string field
                     , const int nx, const int ny, const int nz
                     ){
  std::vector<double> Uv(nx*ny*nz);
  try {
    H5::DataSet dset = h5file.openDataSet(field);
    H5::DataSpace dspace = dset.getSpace();
    dset.read(Uv.data(), H5::PredType::NATIVE_DOUBLE, dspace, dspace);
  } catch (const H5::Exception&){
    partrac::fail(h5file.getFileName(), ": cannot read the dataset '", field, "'");
  }
  for (int ix=0; ix<nx; ++ix){
    for (int iy=0; iy<ny; ++iy){
      for (int iz=0; iz<nz; ++iz){
  	    u(ix, iy, iz) = Uv[nx*ny*iz+nx*iy+ix];
      }
    }
  }
}

inline void load_int_field_as_bool(H5::H5File &h5file
                                 , Grid3<unsigned char>& u
                                 , const std::string field
                                 , const int nx
                                 , const int ny
                                 , const int nz
)
{
  std::vector<int> Uv(nx*ny*nz);
  try {
    H5::DataSet dset = h5file.openDataSet(field);
    H5::DataSpace dspace = dset.getSpace();
    dset.read(Uv.data(), H5::PredType::NATIVE_INT, dspace, dspace);
  } catch (const H5::Exception&){
    partrac::fail(h5file.getFileName(), ": cannot read the dataset '", field, "'");
  }
  for (int ix=0; ix<nx; ++ix){
    for (int iy=0; iy<ny; ++iy){
      for (int iz=0; iz<nz; ++iz){
  	    u(ix, iy, iz) = Uv[nx*ny*iz+nx*iy+ix] != 0;
      }
    }
  }
}

inline void load_h5(const std::string h5filename
                  , Grid3<double>& ux
                  , Grid3<double>& uy
                  , Grid3<double>& uz
                  , Grid3<double>& rho
                  , Grid3<double>& p
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
  H5::H5File h5file;
  try {
    h5file.openFile(h5filename, H5F_ACC_RDONLY);
  } catch (const H5::Exception&){
    partrac::fail("cannot open ", h5filename);
  }
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

inline double weighted_sum(const Grid3<double>& C,
                           const Uint ind[3][2],
                           const double w[2][2][2]){  
  double f = 0.0;
  for (Uint q0=0; q0<2; ++q0){
    for (Uint q1=0; q1<2; ++q1){
      for (Uint q2=0; q2<2; ++q2){
        f += C(ind[0][q0], ind[1][q1], ind[2][q2])*w[q0][q1][q2];
      }
    }
  }
  return f;
}

inline double inner_product(const double ww[2][2][2], const double bb[2][2][2]){
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
inline void matrix_product(double v[3][3][3], const Grid3<double>& u, const Uint ind[3][2], const double W[3][3][3][2][2][2]){
  for (Uint i=0; i<3; ++i){
    for (Uint j=0; j<3; ++j){
      for (Uint k=0; k<3; ++k){
        v[i][j][k] = 0.;
        for (Uint l=0; l<2; ++l){
          for (Uint m=0; m<2; ++m){
            for (Uint n=0; n<2; ++n){
              v[i][j][k] += W[i][j][k][l][m][n] * u(ind[0][l], ind[1][m], ind[2][n]);
            }
          }
        }
      }
    }
  }
}

inline void compute_solid_local(bool is_solid_3[3][3][3], const Grid3<unsigned char>& isSolid, const Uint ind[3][2]){
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
        if (isSolid(ind[0][i], ind[1][j], ind[2][k])){
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

inline void enforce_noslip(double V[2][2][2], const bool is_solid_2[2][2][2]){
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
inline void get_subcube(T V[2][2][2], const T v[3][3][3], const bool sub_x[3]){
  for (Uint i=0; i<2; ++i){
    for (Uint j=0; j<2; ++j){
      for (Uint k=0; k<2; ++k){
        V[i][j][k] = v[ sub_x[0] + i][ sub_x[1] + j ][ sub_x[2] + k ];
      }
    }
  }
}

inline void compute_velocity_subcube(double V[2][2][2], const Grid3<double>& u, const bool is_solid_2[2][2][2], const bool sub_x[3], const Uint ind[3][2], const double W[3][3][3][2][2][2]){
  double v[3][3][3];
  matrix_product(v, u, ind, W);
  get_subcube(V, v, sub_x);
  enforce_noslip(V, is_solid_2);
}

// A region of the lattice, x never wrapped: a floor cell whose nodes are all
// fluid, or one of the eight sub-cubes of a floor cell next to a solid node,
// whose nearest node is its corner on that side
struct LatticeRegion {
  int id = -1;                      // the wrapped node: the floor cell's, or the sub-cube's nearest
  std::array<int, 3> f = {0, 0, 0}; // the floor cell
  int sub = -1;                     // -1: the floor cell; else bit a set for the upper half along axis a
};

// A floor cell's or a sub-cube's velocity nodes at both stamps, [stamp][component][i][j][k]:
// a sub-cube's from the refined nodes with no slip, beside its solid corners
struct LatticeNodes {
  double v[2][3][2][2][2];
  bool solid[2][2][2];
  bool bulk;
};

// The felbm lattice: loading, stamps, locate and reflect; the two below evaluate on it
class StructuredLattice
  : public Interpol {
public:
  StructuredLattice(const std::string& infilename, const std::string& interpolation);
  void update(const double t);
  // The two stamps blended at t, held for every later time
  void freeze(const double t) override;
  bool locate(const Vector3d &x, const double t, CellPos& pos);
  bool reflect(const Vector3d& x, Vector3d& dx, CellPos& pos);
  void enable_reflection() { can_reflect = true; };
  // Outward unit normal toward the solid nodes next to x's node; zero away from walls
  Vector3d get_boundary_normal(const Vector3d &x, int& cell_id);
  double hmin() const { return dx.minCoeff(); };
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
  // The weights of a floor cell from its 1D weights
  void bulk_weights(const double _wq[3][2], double _w[2][2][2], double _dw_x[2][2][2], double _dw_y[2][2][2],
                    double _dw_z[2][2][2]) const;
  // probe_space_boundary's weights of a sub-cube from its 1D weights, for the held evaluation
  void sub_weights(const double _wq[3][2], const bool _is_solid_2[2][2][2], double _wux[2][2][2],
                   double _wuy[2][2][2], double _wuz[2][2][2], double _dwux_x[2][2][2], double _dwux_y[2][2][2],
                   double _dwux_z[2][2][2], double _dwuy_x[2][2][2], double _dwuy_y[2][2][2],
                   double _dwuy_z[2][2][2], double _dwuz_x[2][2][2], double _dwuz_y[2][2][2],
                   double _dwuz_z[2][2][2]) const;
  double get_t_min() { return ts.get_t_min(); };
  double get_t_max() { return ts.get_t_max(); };
  double next_stamp_after(const double t) const override { return ts.next_after(t); }
  using Interpol::locate;
  using Interpol::evaluate;
protected:
  Timestamps ts;
  double t_prev = 0.;
  double t_next = 0.;

  Uint n[3] = {0, 0, 0};
  Vector3d dx;

  double dwq[3][2];

  GridBlock<unsigned char> solid_;   // not bool (vector<bool>)
  GridBlock<double> fields_;
  Grid3<unsigned char> isSolid;
  Grid3<double> ux_prev, uy_prev, uz_prev;
  Grid3<double> ux_next, uy_next, uz_next;
  Grid3<double> rho_prev, rho_next;
  Grid3<double> p_prev, p_next;

  double W[3][3][3][2][2][2];

  partrac::Params felbm_params;

  bool ignore_density = false;
  bool ignore_pressure = false;
  bool ignore_uz = false;
};

// Trilinear in space, no slip at the walls
class StructuredInterpol final
  : public StructuredLattice {
public:
  StructuredInterpol(const std::string& infilename) : StructuredLattice(infilename, "linear") {}
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  using Interpol::evaluate;
  // Regions (regions.hpp's interface): levels the region's 1D weights, the
  // lower face's then the upper face's along each axis; a point on an upper
  // face is the next region's
  using region_type = LatticeRegion;
  using levels_type = std::array<double, 6>;
  static constexpr int n_levels = 6;
  static constexpr bool half_open = true;
  bool region_of(const int id, const Vector3d& x, const double band, LatticeRegion& R, levels_type& lev) const;
  LatticeRegion region_in(const int id, const Vector3d& x) const;
  Vector3d region_point(const LatticeRegion& R, const Vector3d& x, levels_type& lev) const;
  int across(const LatticeRegion& R, const int k, const Vector3d& x, LatticeRegion& next) const;
  void level_rates(const LatticeRegion& R, const Vector3d& v, levels_type& rate) const;
  Vector3d level_grad(const LatticeRegion& R, const int k) const;
  bool is_wall(const LatticeRegion& R, const int k) const;
  // The sub-cubes'
  double region_size(const int) const { return 0.5*hmin(); }
  // Ids a region takes: its node's
  Uint cell_count() const { return n[0]*n[1]*n[2]; }
  using Held = LatticeNodes;
  void hold(const LatticeRegion& R, Held& h) const;
  // evaluate's velocity, its rate and gradient at the levels of the held region
  void held_motion(const int id, const levels_type& lev, const double t, const Held& h, PointValues& fields) const;
  // The same without the gradient
  void held_velocity(const int id, const levels_type& lev, const double t, const Held& h, PointValues& fields) const;
private:
  // The 1D weights of a floor cell, and of the region R at x
  void bulk_wq(const Vector3d& x, const int ix_fl[3], double wq[3][2]) const;
  void region_wq(const LatticeRegion& R, const Vector3d& x, double wq[3][2]) const;
  // The wrapped node indices of floor cell f
  void cell_ind(const std::array<int, 3>& f, Uint ind[3][2]) const;
  bool fluid_cell(const std::array<int, 3>& f) const;
  int node_id(const std::array<int, 3>& f, const int sub) const;
  bool solid_node(const std::array<int, 3>& f, const int sub) const;
  // The region holding x: the floor cell, else the sub-cube of x's halves
  LatticeRegion region_at(const Vector3d& x) const;
  // At a face shared by sub-cubes: the one whose node is id, if its levels at x are down to -tol
  bool region_named(const int id, const Vector3d& x, const double tol, LatticeRegion& R) const;
  void gather(const Uint ind[3][2], const bool bulk, const bool sub_x[3], LatticeNodes& h) const;
  // A weighted sum of stamp s's component c of the gathered nodes
  struct NodeVals {
    const double (&v)[2][3][2][2][2];
    template<int s, int c>
    double sum(const double w[2][2][2]) const { return inner_product(w, v[s][c]); }
  };
  // Velocity, its rate and gradient from the 1D weights and the gathered nodes; Velocity: no gradient
  template<bool Velocity>
  void motion_at(const double wq[3][2], const LatticeNodes& h, const double alpha_t, PointValues& fields) const;
};

// The nearest node in space; no gradient
class StructuredConstInterpol final
  : public StructuredLattice {
public:
  StructuredConstInterpol(const std::string& infilename) : StructuredLattice(infilename, "constant") {}
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  void check_gradient() const;
  using Interpol::evaluate;
};

inline StructuredLattice::StructuredLattice(const std::string& infilename, const std::string& interpolation)
  : Interpol(infilename) {
  felbm_params = partrac::parse_file_or_exit(felbm_schema(), infilename);
  std::cout << "Chosen parameters:" << std::endl;
  felbm_params.print();

  if (felbm_params.get<bool>("ignore_pressure")){
    ignore_pressure = true;
  }
  if (felbm_params.get<bool>("ignore_density")){
    ignore_density = true;
  }
  if (felbm_params.get<bool>("ignore_uz")){
    ignore_uz = true;
  }
  if (felbm_params.get<std::string>("interpolation") != interpolation){
    partrac::fail("felbm_params.dat: interpolation=", felbm_params.get<std::string>("interpolation"),
                  ", built as ", interpolation);
  }

  std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));
  ts.initialize(get_folder() + "/" + felbm_params.get<std::string>("timestamps"));

  std::string solid_filename = get_folder() + "/" + felbm_params.get<std::string>("is_solid_file");
  verify_file_exists(solid_filename);

  H5::H5File solid_file;
  hsize_t dims[3];
  try {
    solid_file.openFile(solid_filename, H5F_ACC_RDONLY);
    H5::DataSpace dspace_solid = solid_file.openDataSet("is_solid").getSpace();
    if (dspace_solid.getSimpleExtentNdims() != 3)
      partrac::fail(solid_filename, ": is_solid is not a three-dimensional array");
    dspace_solid.getSimpleExtentDims(dims, NULL);
  } catch (const H5::Exception&){
    partrac::fail(solid_filename, ": cannot read the dataset 'is_solid'");
  }
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
  solid_.resize(n[0], n[1], n[2], 1);
  isSolid = solid_.field(0);
  fields_.resize(n[0], n[1], n[2], 10);
  Uint f = 0;
  for (Grid3<double>* g : {&ux_prev, &ux_next, &uy_prev, &uy_next, &uz_prev, &uz_next,
                           &rho_prev, &rho_next, &p_prev, &p_next})
    *g = fields_.field(f++);

  load_int_field_as_bool(solid_file, isSolid, "is_solid", n[0], n[1], n[2]);
  solid_file.close();

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

inline Vector3d StructuredLattice::get_boundary_normal(const Vector3d &x, int&){
  return lattice_wall_normal(x, dx, n, [this](const Uint i, const Uint j, const Uint k){
    return isSolid(i, j, k) != 0;
  });
}

inline void StructuredLattice::freeze(const double t){
  update(t);
  const double a = stamp_weight(t, t_prev, t_next);
  const auto hold = [&](Grid3<double>& prev, Grid3<double>& next){
    for (Uint i=0; i<n[0]; ++i)
      for (Uint j=0; j<n[1]; ++j)
        for (Uint k=0; k<n[2]; ++k){
          if (a != 0.)
            prev(i, j, k) = a*next(i, j, k) + (1 - a)*prev(i, j, k);
          next(i, j, k) = prev(i, j, k);
        }
  };
  hold(ux_prev, ux_next);
  hold(uy_prev, uy_next);
  if (!ignore_uz)
    hold(uz_prev, uz_next);
  if (!ignore_density)
    hold(rho_prev, rho_next);
  if (!ignore_pressure)
    hold(p_prev, p_next);
  t_prev = t;
  t_next = t;
  std::cout << "Fields frozen at t = " << t << std::endl;
}

inline void StructuredLattice::update(const double t){
  StampPair sp = ts.get(t);

  // Always load once; past the last stamp both ends are the last stamp, held
  if (!is_initialized || t_prev != sp.prev.t || t_next != sp.next.t){
    if (is_initialized && t_next == sp.prev.t){
      std::swap(ux_prev, ux_next);
      std::swap(uy_prev, uy_next);
      if (!ignore_uz)
        std::swap(uz_prev, uz_next);
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
  t_update = t;
}


inline bool StructuredLattice::locate(const Vector3d &x, const double t, CellPos& pos){
  Uint _ind_pc[3];
  compute_ind_pc(_ind_pc, x, dx, n);
  return !isSolid(_ind_pc[0], _ind_pc[1], _ind_pc[2]);
}

// Walk the node lattice; a wall is the mid-plane between a fluid and a solid node
__attribute__((noinline))
inline bool StructuredLattice::reflect(const Vector3d& x, Vector3d& dx_move, CellPos& pos){
  constexpr int max_bounces = 8;
  constexpr int max_crossings = 4096;
  // Lattice units, unwrapped node indices
  Vector3d p = x.cwiseQuotient(dx);
  Vector3d d = dx_move.cwiseQuotient(dx);
  int idx[3];
  for (Uint i=0; i<3; ++i)
    idx[i] = round_to_int(p[i]);
  const auto solid = [&](){
    return isSolid(imodulo(idx[0], n[0]), imodulo(idx[1], n[1]), imodulo(idx[2], n[2]));
  };
  if (solid())
    return false;
  Vector3d walked = Vector3d::Zero();
  int bounces = 0;
  for (int crossing = 0; crossing < max_crossings; ++crossing){
    // Nearest mid-plane ahead
    int axis = -1;
    double s = 1.;
    for (int i=0; i<3; ++i){
      if (d[i] == 0.) continue;
      const double plane = idx[i] + (d[i] > 0. ? 0.5 : -0.5);
      const double si = std::max((plane - p[i])/d[i], 0.);
      if (si < s){ s = si; axis = i; }
    }
    if (axis < 0){
      dx_move = (walked + d).cwiseProduct(dx);
      return StructuredLattice::locate(x + dx_move, t_update, pos);
    }
    const Vector3d part = s*d;
    p += part;
    walked += part;
    d -= part;
    const int step = d[axis] > 0. ? 1 : -1;
    idx[axis] += step;
    if (solid()){
      idx[axis] -= step;
      if (++bounces > max_bounces)
        return false;
      d[axis] = -d[axis];
    }
  }
  return false;
}

// The nearest node's values, blended in time
inline void StructuredConstInterpol::evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields){
  const double alpha_t = stamp_weight(t, t_prev, t_next);
  Uint i[3];
  compute_ind_pc(i, x, dx, n);
  const Vector3d U_prev(ux_prev(i[0], i[1], i[2]), uy_prev(i[0], i[1], i[2]), uz_prev(i[0], i[1], i[2]));
  const Vector3d U_next(ux_next(i[0], i[1], i[2]), uy_next(i[0], i[1], i[2]), uz_next(i[0], i[1], i[2]));
  fields.U = alpha_t * U_next + (1-alpha_t) * U_prev;
  fields.A = stamp_rate(U_next, U_prev, t_prev, t_next);
  fields.P = alpha_t * p_next(i[0], i[1], i[2]) + (1-alpha_t) * p_prev(i[0], i[1], i[2]);
  fields.Rho = alpha_t * rho_next(i[0], i[1], i[2]) + (1-alpha_t) * rho_prev(i[0], i[1], i[2]);
  fields.gradU = Matrix3d::Zero();
}

inline void StructuredConstInterpol::check_gradient() const {
  if (wants_gradient()){
    partrac::fail("felbm_params.dat: interpolation=constant has no velocity gradient, ",
                  "which int_order=2, vectors and tensors need");
  }
}

inline void StructuredInterpol::bulk_wq(const Vector3d& x, const int ix_fl[3], double wq[3][2]) const {
  for (Uint i=0; i<3; ++i){
    double wxi = (x[i]-dx[i]*ix_fl[i])/dx[i];
    wq[i][0] = 1 - wxi;
    wq[i][1] =     wxi;
  }
}

inline void StructuredInterpol::region_wq(const LatticeRegion& R, const Vector3d& x, double wq[3][2]) const {
  if (R.sub < 0){
    bulk_wq(x, R.f.data(), wq);
    return;
  }
  for (Uint i=0; i<3; ++i){
    const double xd = x[i]/dx[i] - R.f[i];
    double wxi = (R.sub >> i) & 1 ? 2 * xd - 1.0: 2 * xd;
    wq[i][0] = 1 - wxi;
    wq[i][1] =     wxi;
  }
}

inline void StructuredInterpol::cell_ind(const std::array<int, 3>& f, Uint ind[3][2]) const {
  for (Uint i=0; i<3; ++i){
    ind[i][0] = imodulo(f[i], n[i]);
    ind[i][1] = imodulo(ind[i][0] + 1, n[i]);
  }
}

inline bool StructuredInterpol::fluid_cell(const std::array<int, 3>& f) const {
  Uint ind[3][2];
  cell_ind(f, ind);
  for (Uint i=0; i<2; ++i)
    for (Uint j=0; j<2; ++j)
      for (Uint k=0; k<2; ++k)
        if (isSolid(ind[0][i], ind[1][j], ind[2][k]))
          return false;
  return true;
}

inline int StructuredInterpol::node_id(const std::array<int, 3>& f, const int sub) const {
  const int s = sub < 0 ? 0 : sub;
  return int((imodulo(f[0] + (s & 1), n[0])*n[1] + imodulo(f[1] + ((s >> 1) & 1), n[1]))*n[2]
             + imodulo(f[2] + ((s >> 2) & 1), n[2]));
}

inline bool StructuredInterpol::solid_node(const std::array<int, 3>& f, const int sub) const {
  return isSolid(imodulo(f[0] + (sub & 1), n[0]), imodulo(f[1] + ((sub >> 1) & 1), n[1]),
                 imodulo(f[2] + ((sub >> 2) & 1), n[2]));
}

inline LatticeRegion StructuredInterpol::region_at(const Vector3d& x) const {
  LatticeRegion R;
  for (Uint i=0; i<3; ++i)
    R.f[i] = floor(x[i]/dx[i]);
  if (!fluid_cell(R.f)){
    R.sub = 0;
    for (Uint i=0; i<3; ++i)
      if (x[i]/dx[i] - R.f[i] >= 0.5) R.sub |= 1 << i;
  }
  R.id = node_id(R.f, R.sub);
  return R;
}

inline bool StructuredInterpol::region_named(const int id, const Vector3d& x, const double tol, LatticeRegion& R) const {
  if (R.sub < 0 || id < 0)
    return false;
  for (int s = 0; s < 8; ++s){
    LatticeRegion S = R;
    S.sub = s;
    S.id = node_id(S.f, s);
    if (S.id != id) continue;
    levels_type lev;
    region_point(S, x, lev);
    if (*std::min_element(lev.begin(), lev.end()) >= -tol){
      R = S;
      return true;
    }
  }
  return false;
}

inline bool StructuredInterpol::region_of(const int id, const Vector3d& x, const double band, LatticeRegion& R,
                                          levels_type& lev) const {
  R = region_at(x);
  if (R.sub >= 0 && solid_node(R.f, R.sub) && !region_named(id, x, band, R))
    return false;
  region_point(R, x, lev);
  return true;
}

inline LatticeRegion StructuredInterpol::region_in(const int id, const Vector3d& x) const {
  LatticeRegion R = region_at(x);
  // A tie on a face, or a point just across a wall: the stored node's sub-cube
  if (R.id != id)
    region_named(id, x, R.sub >= 0 && solid_node(R.f, R.sub) ? std::numeric_limits<double>::infinity() : 0., R);
  return R;
}

inline Vector3d StructuredInterpol::region_point(const LatticeRegion& R, const Vector3d& x, levels_type& lev) const {
  double wq[3][2];
  region_wq(R, x, wq);
  for (int a = 0; a < 3; ++a){
    lev[2*a] = wq[a][1];
    lev[2*a + 1] = wq[a][0];
  }
  return x;
}

inline int StructuredInterpol::across(const LatticeRegion& R, const int k, const Vector3d& x, LatticeRegion& next) const {
  const int a = k >> 1, bit = 1 << a;
  const bool up = k & 1;
  const int entered = 2*a + !up;
  next = R;
  // A sub-cube's mid-plane: the next sub-cube of the cell, a wall if its node is solid
  if (R.sub >= 0 && up != bool(R.sub & bit)){
    next.sub = R.sub ^ bit;
    if (solid_node(next.f, next.sub))
      return -1;
    next.id = node_id(next.f, next.sub);
    return entered;
  }
  next.f[a] += up ? 1 : -1;
  if (fluid_cell(next.f))
    next.sub = -1;
  // Next to a solid: the sub-cube on the face entered, its node on that face
  else if (R.sub >= 0)
    next.sub = R.sub ^ bit;
  else {
    next.sub = up ? 0 : bit;
    for (int b = 0; b < 3; ++b)
      if (b != a && x[b]/dx[b] - next.f[b] >= 0.5) next.sub |= 1 << b;
  }
  next.id = node_id(next.f, next.sub);
  return entered;
}

inline void StructuredInterpol::level_rates(const LatticeRegion& R, const Vector3d& v, levels_type& rate) const {
  const double scale = R.sub < 0 ? 1. : 2.;
  for (int a = 0; a < 3; ++a){
    rate[2*a] = scale*v[a]/dx[a];
    rate[2*a + 1] = -rate[2*a];
  }
}

inline Vector3d StructuredInterpol::level_grad(const LatticeRegion& R, const int k) const {
  Vector3d g = Vector3d::Zero();
  g[k >> 1] = (k & 1 ? -1. : 1.)*(R.sub < 0 ? 1. : 2.)/dx[k >> 1];
  return g;
}

inline bool StructuredInterpol::is_wall(const LatticeRegion& R, const int k) const {
  const int bit = 1 << (k >> 1);
  return R.sub >= 0 && bool(k & 1) != bool(R.sub & bit) && solid_node(R.f, R.sub ^ bit);
}

inline void StructuredInterpol::gather(const Uint ind[3][2], const bool bulk, const bool sub_x[3], LatticeNodes& h) const {
  const Grid3<double>* g[2][3] = {{&ux_prev, &uy_prev, &uz_prev}, {&ux_next, &uy_next, &uz_next}};
  h.bulk = bulk;
  if (bulk){
    for (Uint s=0; s<2; ++s)
      for (Uint c=0; c<3; ++c)
        for (Uint i=0; i<2; ++i)
          for (Uint j=0; j<2; ++j)
            for (Uint k=0; k<2; ++k)
              h.v[s][c][i][j][k] = (*g[s][c])(ind[0][i], ind[1][j], ind[2][k]);
    return;
  }
  bool is_solid_3[3][3][3];
  compute_solid_local(is_solid_3, isSolid, ind);
  get_subcube(h.solid, is_solid_3, sub_x);
  for (Uint s=0; s<2; ++s)
    for (Uint c=0; c<3; ++c)
      compute_velocity_subcube(h.v[s][c], *g[s][c], h.solid, sub_x, ind, W);
}

template<bool Velocity>
__attribute__((always_inline))
inline void StructuredInterpol::motion_at(const double _wq[3][2], const LatticeNodes& h, const double alpha_t,
                                          PointValues& fields) const {
  double Ux_prev, Uy_prev, Uz_prev;
  double Ux_next, Uy_next, Uz_next;

  double Uxx_prev, Uxx_next, Uxy_prev, Uxy_next, Uxz_prev, Uxz_next;
  double Uyx_prev, Uyx_next, Uyy_prev, Uyy_next, Uyz_prev, Uyz_next;
  double Uzx_prev, Uzx_next, Uzy_prev, Uzy_next, Uzz_prev, Uzz_next;

  const NodeVals v{h.v};
  if (h.bulk) // Bulk cell
  {
    double _w[2][2][2], _dw_x[2][2][2], _dw_y[2][2][2], _dw_z[2][2][2];
    bulk_weights(_wq, _w, _dw_x, _dw_y, _dw_z);

    Ux_prev = v.template sum<0, 0>(_w);
    Ux_next = v.template sum<1, 0>(_w);
    Uy_prev = v.template sum<0, 1>(_w);
    Uy_next = v.template sum<1, 1>(_w);
    Uz_prev = v.template sum<0, 2>(_w);
    Uz_next = v.template sum<1, 2>(_w);

    if constexpr (!Velocity){
      Uxx_prev = v.template sum<0, 0>(_dw_x);
      Uxx_next = v.template sum<1, 0>(_dw_x);
      Uxy_prev = v.template sum<0, 0>(_dw_y);
      Uxy_next = v.template sum<1, 0>(_dw_y);
      Uxz_prev = v.template sum<0, 0>(_dw_z);
      Uxz_next = v.template sum<1, 0>(_dw_z);
      Uyx_prev = v.template sum<0, 1>(_dw_x);
      Uyx_next = v.template sum<1, 1>(_dw_x);
      Uyy_prev = v.template sum<0, 1>(_dw_y);
      Uyy_next = v.template sum<1, 1>(_dw_y);
      Uyz_prev = v.template sum<0, 1>(_dw_z);
      Uyz_next = v.template sum<1, 1>(_dw_z);
      Uzx_prev = v.template sum<0, 2>(_dw_x);
      Uzx_next = v.template sum<1, 2>(_dw_x);
      Uzy_prev = v.template sum<0, 2>(_dw_y);
      Uzy_next = v.template sum<1, 2>(_dw_y);
      Uzz_prev = v.template sum<0, 2>(_dw_z);
      Uzz_next = v.template sum<1, 2>(_dw_z);
    }
  }
  else // Close to boundary
  {
    double _wux[2][2][2], _dwux_x[2][2][2], _dwux_y[2][2][2], _dwux_z[2][2][2];
    double _wuy[2][2][2], _dwuy_x[2][2][2], _dwuy_y[2][2][2], _dwuy_z[2][2][2];
    double _wuz[2][2][2], _dwuz_x[2][2][2], _dwuz_y[2][2][2], _dwuz_z[2][2][2];
    sub_weights(_wq, h.solid, _wux, _wuy, _wuz, _dwux_x, _dwux_y, _dwux_z,
                _dwuy_x, _dwuy_y, _dwuy_z, _dwuz_x, _dwuz_y, _dwuz_z);

    Ux_prev = v.template sum<0, 0>(_wux);
    Ux_next = v.template sum<1, 0>(_wux);
    Uy_prev = v.template sum<0, 1>(_wuy);
    Uy_next = v.template sum<1, 1>(_wuy);
    Uz_prev = v.template sum<0, 2>(_wuz);
    Uz_next = v.template sum<1, 2>(_wuz);

    if constexpr (!Velocity){
      Uxx_prev = v.template sum<0, 0>(_dwux_x);
      Uxx_next = v.template sum<1, 0>(_dwux_x);
      Uxy_prev = v.template sum<0, 0>(_dwux_y);
      Uxy_next = v.template sum<1, 0>(_dwux_y);
      Uxz_prev = v.template sum<0, 0>(_dwux_z);
      Uxz_next = v.template sum<1, 0>(_dwux_z);
      Uyx_prev = v.template sum<0, 1>(_dwuy_x);
      Uyx_next = v.template sum<1, 1>(_dwuy_x);
      Uyy_prev = v.template sum<0, 1>(_dwuy_y);
      Uyy_next = v.template sum<1, 1>(_dwuy_y);
      Uyz_prev = v.template sum<0, 1>(_dwuy_z);
      Uyz_next = v.template sum<1, 1>(_dwuy_z);
      Uzx_prev = v.template sum<0, 2>(_dwuz_x);
      Uzx_next = v.template sum<1, 2>(_dwuz_x);
      Uzy_prev = v.template sum<0, 2>(_dwuz_y);
      Uzy_next = v.template sum<1, 2>(_dwuz_y);
      Uzz_prev = v.template sum<0, 2>(_dwuz_z);
      Uzz_next = v.template sum<1, 2>(_dwuz_z);
    }
  }

  fields.U = { alpha_t * Ux_next + (1-alpha_t) * Ux_prev,
               alpha_t * Uy_next + (1-alpha_t) * Uy_prev,
               alpha_t * Uz_next + (1-alpha_t) * Uz_prev };
  fields.A = { stamp_rate(Ux_next, Ux_prev, t_prev, t_next),
               stamp_rate(Uy_next, Uy_prev, t_prev, t_next),
               stamp_rate(Uz_next, Uz_prev, t_prev, t_next) };

  if constexpr (!Velocity){
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
}

inline void StructuredInterpol::evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields){
  // Time weight between stamps
  const double alpha_t = stamp_weight(t, t_prev, t_next);
  // Assuming locate has already been called and found that the cell is not in solid
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
    Rho_prev = rho_prev(_ind_pc[0], _ind_pc[1], _ind_pc[2]);
    Rho_next = rho_next(_ind_pc[0], _ind_pc[1], _ind_pc[2]);

    P_prev = p_prev(_ind_pc[0], _ind_pc[1], _ind_pc[2]);
    P_next = p_next(_ind_pc[0], _ind_pc[1], _ind_pc[2]);

    Uxx_prev = inner_product(_dwux_x, _Vx_prev);
    Uxx_next = inner_product(_dwux_x, _Vx_next);
    Uxy_prev = inner_product(_dwux_y, _Vx_prev);
    Uxy_next = inner_product(_dwux_y, _Vx_next);
    Uxz_prev = inner_product(_dwux_z, _Vx_prev);
    Uxz_next = inner_product(_dwux_z, _Vx_next);
    Uyx_prev = inner_product(_dwuy_x, _Vy_prev);
    Uyx_next = inner_product(_dwuy_x, _Vy_next);
    Uyy_prev = inner_product(_dwuy_y, _Vy_prev);
    Uyy_next = inner_product(_dwuy_y, _Vy_next);
    Uyz_prev = inner_product(_dwuy_z, _Vy_prev);
    Uyz_next = inner_product(_dwuy_z, _Vy_next);
    Uzx_prev = inner_product(_dwuz_x, _Vz_prev);
    Uzx_next = inner_product(_dwuz_x, _Vz_next);
    Uzy_prev = inner_product(_dwuz_y, _Vz_prev);
    Uzy_next = inner_product(_dwuz_y, _Vz_next);
    Uzz_prev = inner_product(_dwuz_z, _Vz_prev);
    Uzz_next = inner_product(_dwuz_z, _Vz_next);
  }

  fields.U = { alpha_t * Ux_next + (1-alpha_t) * Ux_prev,
               alpha_t * Uy_next + (1-alpha_t) * Uy_prev,
               alpha_t * Uz_next + (1-alpha_t) * Uz_prev };
  fields.A = { stamp_rate(Ux_next, Ux_prev, t_prev, t_next),
               stamp_rate(Uy_next, Uy_prev, t_prev, t_next),
               stamp_rate(Uz_next, Uz_prev, t_prev, t_next) };

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

inline void StructuredInterpol::hold(const LatticeRegion& R, Held& h) const {
  Uint ind[3][2];
  cell_ind(R.f, ind);
  const bool sub_x[3] = {bool(R.sub & 1), bool(R.sub & 2), bool(R.sub & 4)};
  gather(ind, R.sub < 0, sub_x, h);
}

inline void StructuredInterpol::held_motion(const int, const levels_type& lev, const double t, const Held& h,
                                            PointValues& fields) const {
  const double wq[3][2] = {{lev[1], lev[0]}, {lev[3], lev[2]}, {lev[5], lev[4]}};
  motion_at<false>(wq, h, stamp_weight(t, t_prev, t_next), fields);
}

inline void StructuredInterpol::held_velocity(const int, const levels_type& lev, const double t, const Held& h,
                                              PointValues& fields) const {
  const double wq[3][2] = {{lev[1], lev[0]}, {lev[3], lev[2]}, {lev[5], lev[4]}};
  motion_at<true>(wq, h, stamp_weight(t, t_prev, t_next), fields);
}

inline bool StructuredLattice::compute_ind(const Vector3d &x, Uint _ind[3][2], int _ix_fl[3]){
  // Assuming this cell is not inside the solid phase
  for (Uint i=0; i<3; ++i){
    _ix_fl[i] = floor(x[i]/dx[i]);
  }

  for (Uint i=0; i<3; ++i){
    _ind[i][0] = imodulo(_ix_fl[i], n[i]);
    _ind[i][1] = imodulo(_ind[i][0] + 1, n[i]);
  }

  for (Uint i=0; i<2; ++i){
    for (Uint j=0; j<2; ++j){
      for (Uint k=0; k<2; ++k){
        if (isSolid(_ind[0][i], _ind[1][j], _ind[2][k])) {
          return false;
        }
      }
    }
  }
  return true;
}

inline void StructuredLattice::probe_space_bulk(const Vector3d &x, 
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

  bulk_weights(_wq, _w, _dw_x, _dw_y, _dw_z);
}

__attribute__((always_inline))
inline void StructuredLattice::bulk_weights(const double _wq[3][2],
    double _w[2][2][2],
    double _dw_x[2][2][2],
    double _dw_y[2][2][2],
    double _dw_z[2][2][2]
  ) const {
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

inline void StructuredLattice::probe_space_boundary(
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

  for (Uint i=0; i<3; ++i){
    _sub_x[i] = xd[i] >= 0.5;
  }

  bool is_solid_3[3][3][3];
  compute_solid_local(is_solid_3, isSolid, _ind);
  get_subcube(_is_solid_2, is_solid_3, _sub_x);

  double _wq[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
  for (Uint i=0; i<3; ++i){
    double wxi = _sub_x[i] ? 2 * xd[i] - 1.0: 2 * xd[i];
    _wq[i][0] = 1 - wxi;
    _wq[i][1] =     wxi;
  }

  double gamma = 2;
  for (Uint i=0; i<2; ++i){
    for (Uint j=0; j<2; ++j){
      for (Uint k=0; k<2; ++k){
        double wqux = _wq[0][i];
        double dwqux = dwq[0][i];
        double wquy = _wq[1][j];
        double dwquy = dwq[1][j];
        double wquz = _wq[2][k];
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

__attribute__((always_inline))
inline void StructuredLattice::sub_weights(const double _wq[3][2],
  const bool _is_solid_2[2][2][2],
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
  ) const
{
  double gamma = 2;
  for (Uint i=0; i<2; ++i){
    for (Uint j=0; j<2; ++j){
      for (Uint k=0; k<2; ++k){
        double wqux = _wq[0][i];
        double dwqux = dwq[0][i];
        double wquy = _wq[1][j];
        double dwquy = dwq[1][j];
        double wquz = _wq[2][k];
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

#endif
