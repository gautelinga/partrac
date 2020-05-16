#include "io.hpp"
#include "Timestamps.hpp"
#include "H5Cpp.h"

#ifndef __INTERPOL_HPP
#define __INTERPOL_HPP

using namespace H5;
using namespace std;

class Interpol {
public:
  Interpol(string infilename);
  ~Interpol() {};
  void update(const double t);
  void set_folder(string folder){ this->folder=folder; };
  void probe(const double x, const double y, const double z);
  bool inside_domain();
  double get_ux();
  double get_uy();
  double get_uz();
  string get_folder(){ return folder; };
  double get_Lx() { return Lx; };
  double get_Ly() { return Ly; };
  double get_Lz() { return Lz; };
  Uint get_nx() { return nx; };
  Uint get_ny() { return ny; };
  Uint get_nz() { return nz; };
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
  double get_divu() { return get_uxx()+get_uyy()+get_uzz(); };
  double get_vortx() { return get_uzy()-get_uyz(); };
  double get_vorty() { return get_uxz()-get_uzx(); };
  double get_vortz() { return get_uyx()-get_uxy(); };
  double get_Jux() { return (get_uxx()*get_ux()
                             + get_uxy()*get_uy()
                             + get_uxz()*get_uz()); };
  double get_Juy() { return (get_uyx()*get_ux()
                             + get_uyy()*get_uy()
                             + get_uyz()*get_uz()); };
  double get_Juz() { return (get_uzx()*get_ux()
                             + get_uzy()*get_uy()
                             + get_uzz()*get_uz()); };
  bool get_nodal_inside(const Uint ix, const Uint iy, const Uint iz){
    return !isSolid[ix][iy][iz];
  }
  double get_nodal_ux(const Uint ix, const Uint iy, const Uint iz){
    return alpha_t * ux_next[ix][iy][iz] + (1-alpha_t) * ux_prev[ix][iy][iz];
  }
  double get_nodal_uy(const Uint ix, const Uint iy, const Uint iz){
    return alpha_t * uy_next[ix][iy][iz] + (1-alpha_t) * uy_prev[ix][iy][iz];
  }
  double get_nodal_uz(const Uint ix, const Uint iy, const Uint iz){
    return alpha_t * uz_next[ix][iy][iz] + (1-alpha_t) * uz_prev[ix][iy][iz];
  }
private:
  Timestamps ts;
  string folder;
  bool is_initialized = false;
  bool verbose = true;
  double t_prev = 0.;
  double t_next = 0.;
  double alpha_t;
  Uint nx;
  Uint ny;
  Uint nz;
  double Lx;
  double Ly;
  double Lz;
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
  Uint ind_pc[3] = {0, 0, 0};  // piecewise constant intp

  double w[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dw_x[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dw_y[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
  double dw_z[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};

  map<array<Uint, 3>, vector<array<Uint, 3>>> solid_ids_neigh;
  double inside_domain_factor = 0.;
  double inside_domain_factor_x = 0.;
  double inside_domain_factor_y = 0.;
  double inside_domain_factor_z = 0.;
  double Ux = 0.;
  double Uy = 0.;
  double Uz = 0.;
};

Interpol::Interpol(string infilename) : ts(infilename) {
  folder = ts.get_folder();
  string solid_filename = ts.get_is_solid();
  verify_file_exists(solid_filename);

  H5File solid_file(solid_filename, H5F_ACC_RDONLY);
  DataSet dset_solid = solid_file.openDataSet("is_solid");
  DataSpace dspace_solid = dset_solid.getSpace();

  hsize_t dims[3];
  dspace_solid.getSimpleExtentDims(dims, NULL);
  this->nx = dims[0];
  this->ny = dims[1];
  this->nz = dims[2];
  Lx = nx;
  Ly = ny;
  Lz = nz;

  // Create arrays
  isSolid = new int**[nx];
  levelZ = new double**[nx];
  ux_prev = new double**[nx];
  uy_prev = new double**[nx];
  uz_prev = new double**[nx];
  ux_next = new double**[nx];
  uy_next = new double**[nx];
  uz_next = new double**[nx];
  rho_prev = new double**[nx];
  rho_next = new double**[nx];
  p_prev = new double**[nx];
  p_next = new double**[nx];
  for (Uint ix=0; ix<nx; ++ix){
    isSolid[ix] = new int*[ny];
    levelZ[ix] = new double*[ny];
    ux_prev[ix] = new double*[ny];
    uy_prev[ix] = new double*[ny];
    uz_prev[ix] = new double*[ny];
    ux_next[ix] = new double*[ny];
    uy_next[ix] = new double*[ny];
    uz_next[ix] = new double*[ny];
    rho_prev[ix] = new double*[ny];
    rho_next[ix] = new double*[ny];
    p_prev[ix] = new double*[ny];
    p_next[ix] = new double*[ny];
    for (Uint iy=0; iy<ny; ++iy){
      isSolid[ix][iy] = new int[nz];
      levelZ[ix][iy] = new double[nz];
      ux_prev[ix][iy] = new double[nz];
      uy_prev[ix][iy] = new double[nz];
      uz_prev[ix][iy] = new double[nz];
      ux_next[ix][iy] = new double[nz];
      uy_next[ix][iy] = new double[nz];
      uz_next[ix][iy] = new double[nz];
      rho_prev[ix][iy] = new double[nz];
      rho_next[ix][iy] = new double[nz];
      p_prev[ix][iy] = new double[nz];
      p_next[ix][iy] = new double[nz];
    }
  }

  Uint nnlist[27][3];
  Uint nnum=0;
  for (int i=-1; i<=1; ++i){
    for (int j=-1; j<=1; ++j){
      for (int k=-1; k<=1; ++k){
        nnlist[nnum][0] = i;
        nnlist[nnum][1] = j;
        nnlist[nnum][2] = k;
        ++nnum;
      }
    }
  }

  // Smooth solid-liquid interface
  load_int_field(solid_file, isSolid, "is_solid", nx, ny, nz);
  for (Uint ix=0; ix<nx; ++ix){
    for (Uint iy=0; iy<ny; ++iy){
      for (Uint iz=0; iz<nz; ++iz){
        levelZ[ix][iy][iz] = -2*isSolid[ix][iy][iz]+1;
        if (isSolid[ix][iy][iz]){
          vector<array<Uint, 3>> neigh;
          for (Uint inn=0; inn<27; ++inn){
            Uint ix_ = imodulo(ix+nnlist[inn][0], nx);
            Uint iy_ = imodulo(iy+nnlist[inn][1], ny);
            Uint iz_ = imodulo(iz+nnlist[inn][2], nz);
            if (!isSolid[ix_][iy_][iz_]){
              neigh.push_back({ix_, iy_, iz_});
            }
          }
          pair<array<Uint, 3>, vector<array<Uint, 3>>> pp({ix, iy, iz}, neigh);
          solid_ids_neigh.insert(pp);
        }
      }
    }
  }
}

void Interpol::update(const double t){
  StampPair sp = ts.get(t);
  // cout << sp.prev.filename << " " << sp.next.filename << endl;

  if (!is_initialized || t_prev != sp.prev.t || t_next != sp.next.t){
    cout << "Prev: Timestep = " << sp.prev.t << ", filename = " << sp.prev.filename << endl;
    load_h5(folder + "/" + sp.prev.filename,
            ux_prev, uy_prev, uz_prev, rho_prev, p_prev,
            nx, ny, nz, verbose);

    cout << "Next: Timestep = " << sp.next.t << ", filename = " << sp.next.filename << endl;
    load_h5(folder + "/" + sp.next.filename,
            ux_next, uy_next, uz_next, rho_next, p_next,
            nx, ny, nz, verbose);

    // Extrapolate fields into solid
    for (map<array<Uint, 3>, vector<array<Uint, 3>>>::iterator it = solid_ids_neigh.begin();
         it != solid_ids_neigh.end(); ++it){
      array<Uint, 3> inds = it->first;
      vector<array<Uint, 3>> neigh = it->second;
      Uint ix = inds[0];
      Uint iy = inds[1];
      Uint iz = inds[2];
      Uint num_neigh = neigh.size();
      if (num_neigh > 0){
        double ux_prev_sum = 0.;
        double uy_prev_sum = 0.;
        double uz_prev_sum = 0.;
        double ux_next_sum = 0.;
        double uy_next_sum = 0.;
        double uz_next_sum = 0.;
        double rho_prev_sum = 0.;
        double rho_next_sum = 0.;
        double p_prev_sum = 0.;
        double p_next_sum = 0.;
        for (vector<array<Uint, 3>>::iterator vit=neigh.begin();
             vit != neigh.end(); ++vit){
          array<Uint, 3> other_inds = *vit;
          Uint ix_ = other_inds[0];
          Uint iy_ = other_inds[1];
          Uint iz_ = other_inds[2];
          ux_prev_sum += ux_prev[ix_][iy_][iz_];
          uy_prev_sum += uy_prev[ix_][iy_][iz_];
          uz_prev_sum += uz_prev[ix_][iy_][iz_];
          ux_next_sum += ux_next[ix_][iy_][iz_];
          uy_next_sum += uy_next[ix_][iy_][iz_];
          uz_next_sum += uz_next[ix_][iy_][iz_];
          rho_prev_sum += rho_prev[ix_][iy_][iz_];
          rho_next_sum += rho_next[ix_][iy_][iz_];
          p_prev_sum += p_prev[ix_][iy_][iz_];
          p_next_sum += p_next[ix_][iy_][iz_];
        }
        ux_prev[ix][iy][iz] = ux_prev_sum/num_neigh;
        uy_prev[ix][iy][iz] = uy_prev_sum/num_neigh;
        uz_prev[ix][iy][iz] = uz_prev_sum/num_neigh;
        ux_next[ix][iy][iz] = ux_next_sum/num_neigh;
        uy_next[ix][iy][iz] = uy_next_sum/num_neigh;
        uz_next[ix][iy][iz] = uz_next_sum/num_neigh;
        rho_prev[ix][iy][iz] = rho_prev_sum/num_neigh;
        rho_next[ix][iy][iz] = rho_next_sum/num_neigh;
        p_prev[ix][iy][iz] = p_prev_sum/num_neigh;
        p_next[ix][iy][iz] = p_next_sum/num_neigh;
      }
    }

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  alpha_t = sp.weight_next(t);
}

void Interpol::probe(const double x,
                     const double y,
                     const double z){
  double dx = Lx/nx;
  double dy = Ly/ny;
  double dz = Lz/nz;

  int ix_fl = floor(x/dx);
  int iy_fl = floor(y/dy);
  int iz_fl = floor(z/dz);

  ind[0][0] = imodulo(ix_fl, nx);
  ind[0][1] = imodulo(ind[0][0] + 1, nx);
  ind[1][0] = imodulo(iy_fl, ny);
  ind[1][1] = imodulo(ind[1][0] + 1, ny);
  ind[2][0] = imodulo(iz_fl, nz);
  ind[2][1] = imodulo(ind[2][0] + 1, nz);

  ind_pc[0] = imodulo(round(x/dx), nx);
  ind_pc[1] = imodulo(round(y/dy), ny);
  ind_pc[2] = imodulo(round(z/dz), nz);

  double wx = (x-dx*ix_fl)/dx;
  double wy = (y-dy*iy_fl)/dy;
  double wz = (z-dz*iz_fl)/dz;

  double wq[3][2] = {{1-wx, wx},
                     {1-wy, wy},
                     {1-wz, wz}};

  double dwq[3][2] = {{-1./dx, 1./dx},
                      {-1./dy, 1./dy},
                      {-1./dz, 1./dz}};

  for (Uint i=0; i<2; ++i){
    for (Uint j=0; j<2; ++j){
      for (Uint k=0; k<2; ++k){
        w[i][j][k] = wq[0][i]*wq[1][j]*wq[2][k];
        dw_x[i][j][k] = dwq[0][i]*wq[1][j]*wq[2][k];
        dw_y[i][j][k] = wq[0][i]*dwq[1][j]*wq[2][k];
        dw_z[i][j][k] = wq[0][i]*wq[1][j]*dwq[2][k];
      }
    }
  }

  // Factor used to determine whether point is in the domain or not
  // and to strictly enforce no-slip boundary conditions
  inside_domain_factor = max(weighted_sum(levelZ, ind, w), 0.0);
  inside_domain_factor_x = weighted_sum(levelZ, ind, dw_x);
  inside_domain_factor_y = weighted_sum(levelZ, ind, dw_y);
  inside_domain_factor_z = weighted_sum(levelZ, ind, dw_z);

  // Precompute velocities
  double Ux_prev = weighted_sum(ux_prev, ind, w);
  double Ux_next = weighted_sum(ux_next, ind, w);
  Ux = alpha_t * Ux_next + (1-alpha_t) * Ux_prev;

  double Uy_prev = weighted_sum(uy_prev, ind, w);
  double Uy_next = weighted_sum(uy_next, ind, w);
  Uy = alpha_t * Uy_next + (1-alpha_t) * Uy_prev;

  double Uz_prev = weighted_sum(uz_prev, ind, w);
  double Uz_next = weighted_sum(uz_next, ind, w);
  Uz = alpha_t * Uz_next + (1-alpha_t) * Uz_prev;
}

bool Interpol::inside_domain(){
  return inside_domain_factor > 0.0;
}

// Interpolate in space and time and enforce BCs
double Interpol::get_ux(){
  return Ux*inside_domain_factor;
}

double Interpol::get_uy(){
  return Uy*inside_domain_factor;
}

double Interpol::get_uz(){
  return Uz*inside_domain_factor;
}

double Interpol::get_rho(){
  double Rho_prev = weighted_sum(rho_prev, ind, w);
  double Rho_next = weighted_sum(rho_next, ind, w);
  //double Rho_prev = rho_prev[ind_pc[0]][ind_pc[1]][ind_pc[2]];
  //double Rho_next = rho_next[ind_pc[0]][ind_pc[1]][ind_pc[2]];
  return alpha_t * Rho_next + (1-alpha_t) * Rho_prev;
}

double Interpol::get_p(){
  double P_prev = weighted_sum(p_prev, ind, w);
  double P_next = weighted_sum(p_next, ind, w);
  return alpha_t * P_next + (1-alpha_t) * P_prev;
}

double Interpol::get_uxx(){
  double Uxx_prev = weighted_sum(ux_prev, ind, dw_x);
  double Uxx_next = weighted_sum(ux_next, ind, dw_x);
  double Uxx = alpha_t * Uxx_next + (1-alpha_t) * Uxx_prev;
  return Ux * inside_domain_factor_x + Uxx * inside_domain_factor;
}

double Interpol::get_uxy(){
  double Uxy_prev = weighted_sum(ux_prev, ind, dw_y);
  double Uxy_next = weighted_sum(ux_next, ind, dw_y);
  double Uxy = alpha_t * Uxy_next + (1-alpha_t) * Uxy_prev;
  return Ux * inside_domain_factor_y + Uxy * inside_domain_factor;
}

double Interpol::get_uxz(){
  double Uxz_prev = weighted_sum(ux_prev, ind, dw_z);
  double Uxz_next = weighted_sum(ux_next, ind, dw_z);
  double Uxz = alpha_t * Uxz_next + (1-alpha_t) * Uxz_prev;
  return Ux * inside_domain_factor_z + Uxz * inside_domain_factor;
}

double Interpol::get_uyx(){
  double Uyx_prev = weighted_sum(uy_prev, ind, dw_x);
  double Uyx_next = weighted_sum(uy_next, ind, dw_x);
  double Uyx = alpha_t * Uyx_next + (1-alpha_t) * Uyx_prev;
  return Uy * inside_domain_factor_x + Uyx * inside_domain_factor;
}

double Interpol::get_uyy(){
  double Uyy_prev = weighted_sum(uy_prev, ind, dw_y);
  double Uyy_next = weighted_sum(uy_next, ind, dw_y);
  double Uyy = alpha_t * Uyy_next + (1-alpha_t) * Uyy_prev;
  return Uy * inside_domain_factor_y + Uyy * inside_domain_factor;
}

double Interpol::get_uyz(){
  double Uyz_prev = weighted_sum(uy_prev, ind, dw_z);
  double Uyz_next = weighted_sum(uy_next, ind, dw_z);
  double Uyz = alpha_t * Uyz_next + (1-alpha_t) * Uyz_prev;
  return Uy * inside_domain_factor_z + Uyz * inside_domain_factor;
}

double Interpol::get_uzx(){
  double Uzx_prev = weighted_sum(uz_prev, ind, dw_x);
  double Uzx_next = weighted_sum(uz_next, ind, dw_x);
  double Uzx = alpha_t * Uzx_next + (1-alpha_t) * Uzx_prev;
  return Uz * inside_domain_factor_x + Uzx * inside_domain_factor;
}

double Interpol::get_uzy(){
  double Uzy_prev = weighted_sum(uz_prev, ind, dw_y);
  double Uzy_next = weighted_sum(uz_next, ind, dw_y);
  double Uzy = alpha_t * Uzy_next + (1-alpha_t) * Uzy_prev;
  return Uz * inside_domain_factor_y + Uzy * inside_domain_factor;
}

double Interpol::get_uzz(){
  double Uzz_prev = weighted_sum(uz_prev, ind, dw_z);
  double Uzz_next = weighted_sum(uz_next, ind, dw_z);
  double Uzz = alpha_t * Uzz_next + (1-alpha_t) * Uzz_prev;
  return Uz * inside_domain_factor_z + Uzz * inside_domain_factor;
}

#endif
