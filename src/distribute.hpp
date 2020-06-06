#include <iomanip>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include "Interpol.hpp"
#include "utils.hpp"

#ifndef __DISTRIBUTE_HPP
#define __DISTRIBUTE_HPP

using namespace std;

struct less_than_op {
  inline bool operator() (const Vector3d &a, const Vector3d &b){
    return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]) || (a[0] == b[0] && a[1] == b[1] && a[2] < b[2]);
  }
};

void load_positions(string input_file,
                    vector<Vector3d> &pos_init,
                    const Uint Nrw){
  ifstream infile(input_file);
  double x, y, z;
  while (infile >> x >> y >> z){
    pos_init.push_back({x, y, z});
  }
  infile.close();
  if (Nrw != pos_init.size()){
    cout << "Wrong dimensions..." << endl;
    exit(0);
  }
}

void dump_positions(string output_file,
                    Vector3d* x_rw,
                    const Uint Nrw){
  ofstream outfile(output_file);
  for (Uint irw=0; irw<Nrw; ++irw){
    outfile << std::setprecision(12)
            << x_rw[irw][0] << " "
            << x_rw[irw][1] << " "
            << x_rw[irw][2] << endl;
  }
  outfile.close();
}

void load_faces(string input_file,
                FacesType& faces){
  ifstream infile(input_file);
  Uint first, second, third;
  double fourth;
  while (infile >> first >> second >> third >> fourth){
    faces.push_back({{first, second, third}, fourth});
  }
  infile.close();
}

void load_edges(string input_file,
                EdgesType &edges){
  ifstream infile(input_file);
  Uint first, second;
  double third;
  while (infile >> first >> second >> third){
    edges.push_back({{first, second}, third});
  }
  infile.close();
}

void dump_faces(string output_file,
                FacesType &faces){
  ofstream outfile(output_file);
  for (FacesType::iterator faceit = faces.begin();
       faceit != faces.end(); ++faceit){
    outfile << faceit->first[0] << " " << faceit->first[1] << " " << faceit->first[2]
            << " " << faceit->second << endl;
  }
  outfile.close();
}

void dump_edges(string output_file,
                EdgesType &edges){
  ofstream outfile(output_file);
  for (EdgesType::iterator edgeit = edges.begin();
       edgeit != edges.end(); ++edgeit){
    outfile << edgeit->first[0] << " " << edgeit->first[1] << " " << edgeit->second << endl;
  }
  outfile.close();
}

void load_colors(string input_file,
                 double* c_rw, const Uint Nrw){
  ifstream infile(input_file);
  for (Uint irw=0; irw < Nrw; ++irw){
    infile >> c_rw[irw];
  }
  infile.close();
}

void dump_colors(string output_file,
                 double* c_rw, const Uint Nrw){
  ofstream outfile(output_file);
  for (Uint irw=0; irw < Nrw; ++irw){
    outfile << setprecision(12) << c_rw[irw] << endl;
  }
  outfile.close();
}

vector<Vector3d> initial_positions(const string init_mode,
                                   const string init_weight,
                                   Uint &Nrw,
                                   const Vector3d &x0,
                                   const double La,
                                   const double Lb,
                                   const double ds,
                                   Interpol &intp,
                                   mt19937 &gen,
                                   EdgesType &edges,
                                   FacesType &faces
                                   ){
  Vector3d x;
  double Lx = intp.get_Lx();
  double Ly = intp.get_Ly();
  double Lz = intp.get_Lz();

  Uint Nx = 1;
  Uint Ny = 1;
  Uint Nz = 1;

  if (init_mode == "sheet_xy" ||
      init_mode == "sheet_xz" ||
      init_mode == "sheet_yz"){
    vector<Vector3d> pos_init;

    Vector3d n(0., 0., 0.);
    Vector3d ta(0., 0., 0.);
    Vector3d tb(0., 0., 0.);
    if (init_mode == "sheet_xy"){
      n[2] = 1.0;
      ta[0] = 1.0;
      tb[1] = 1.0;
    }
    if (init_mode == "sheet_xz"){
      n[1] = 1.0;
      ta[0] = 1.0;
      tb[2] = 1.0;
    }
    if (init_mode == "sheet_yz"){
      n[0] = 1.0;
      ta[1] = 1.0;
      tb[2] = 1.0;
    }

    Vector3d x00 = x0;
    x00[0] += - La/2*ta[0] - Lb/2*tb[0];
    x00[1] += - La/2*ta[1] - Lb/2*tb[1];
    x00[2] += - La/2*ta[2] - Lb/2*tb[2];

    Vector3d x01 = x0;
    x01[0] += La/2*ta[0] + Lb/2*tb[0];
    x01[1] += La/2*ta[1] - Lb/2*tb[1];
    x01[2] += - La/2*ta[2] - Lb/2*tb[2];

    Vector3d x10 = x0;
    x10[0] += La/2*ta[0] + Lb/2*tb[0];
    x10[1] += La/2*ta[1] + Lb/2*tb[1];
    x10[2] += La/2*ta[2] + Lb/2*tb[2];

    Vector3d x11 = x0;
    x11[0] += - La/2*ta[0] - Lb/2*tb[0];
    x11[1] += - La/2*ta[1] + Lb/2*tb[1];
    x11[2] += - La/2*ta[2] + Lb/2*tb[2];

    intp.probe(x00);
    bool inside_00 = intp.inside_domain();
    intp.probe(x01);
    bool inside_01 = intp.inside_domain();
    intp.probe(x10);
    bool inside_10 = intp.inside_domain();
    intp.probe(x11);
    bool inside_11 = intp.inside_domain();
    if (inside_00 && inside_01 && inside_10 && inside_11){
      cout << "Sheet inside domain." << endl;
    }
    else {
      cout << "Sheet not inside domain" << endl;
      exit(0);
    }

    pos_init.push_back(x00);
    pos_init.push_back(x01);
    pos_init.push_back(x10);
    pos_init.push_back(x11);

    edges.push_back({{0, 1}, dist(pos_init[0], pos_init[1])});
    edges.push_back({{0, 2}, dist(pos_init[0], pos_init[2])});
    edges.push_back({{1, 2}, dist(pos_init[1], pos_init[2])});
    edges.push_back({{2, 3}, dist(pos_init[2], pos_init[3])});
    edges.push_back({{3, 0}, dist(pos_init[3], pos_init[4])});

    faces.push_back({{0, 2, 1}, La*Lb/2});
    faces.push_back({{1, 3, 4}, La*Lb/2});

    Nrw = pos_init.size();
    return pos_init;
  }
  else if (init_mode == "pair_xyz" ||
           init_mode == "pair_xy" ||
           init_mode == "pair_xz" ||
           init_mode == "pair_yz" ||
           init_mode == "pairs_xyz_x" ||
           init_mode == "pairs_xyz_y" ||
           init_mode == "pairs_xyz_z" ||
           init_mode == "pairs_xy_x" ||
           init_mode == "pairs_xy_y" ||
           init_mode == "pairs_xy_z" ||
           init_mode == "pairs_xz_x" ||
           init_mode == "pairs_xz_y" ||
           init_mode == "pairs_xz_z" ||
           init_mode == "pairs_yz_x" ||
           init_mode == "pairs_yz_y" ||
           init_mode == "pairs_yz_z"
           ){
    vector<Vector3d> pos_init;

    uniform_real_distribution<> uni_dist_x(0, Lx);
    uniform_real_distribution<> uni_dist_y(0, Ly);
    uniform_real_distribution<> uni_dist_z(0, Lz);
    normal_distribution<double> rnd_normal(0.0, 1.0);

    // uniform_real_distribution<> uni_dist_theta(0., 2*M_PI);
    Uint Npairs = (init_mode == "pair_xyz" ||
                   init_mode == "pair_xy" ||
                   init_mode == "pair_xz" ||
                   init_mode == "pair_yz") ? 1 : Nrw/2;

    cout << "Npairs = " << Npairs << endl;

    Vector3d x0_ = x0;
    Uint ipair=0;
    while (ipair < Npairs){
      if (init_mode == "pairs_xyz_x" ||
          init_mode == "pairs_xy_x" ||
          init_mode == "pairs_xz_x" ||
          init_mode == "pairs_yz_x"){
        x0_[0] = uni_dist_x(gen);
      }
      if (init_mode == "pairs_xyz_y" ||
          init_mode == "pairs_xy_y" ||
          init_mode == "pairs_xz_y" ||
          init_mode == "pairs_yz_y"){
        x0_[0] = uni_dist_y(gen);
      }
      if (init_mode == "pairs_xyz_z" ||
          init_mode == "pairs_xy_z" ||
          init_mode == "pairs_xz_z" ||
          init_mode == "pairs_yz_z"){
        x0_[0] = uni_dist_z(gen);
      }

      Vector3d dx(0., 0., 0.);
      if (init_mode == "pair_yz" ||
          init_mode == "pairs_yz_x" ||
          init_mode == "pairs_yz_y" ||
          init_mode == "pairs_yz_z"){
        dx[0] = 0.;
      }
      else {
        dx[0] = rnd_normal(gen);
      }
      if (init_mode == "pair_xz" ||
          init_mode == "pairs_xz_x" ||
          init_mode == "pairs_xz_y" ||
          init_mode == "pairs_xz_z"){
        dx[1] = 0.;
      }
      else {
        dx[1] = rnd_normal(gen);
      }
      if (init_mode == "pair_xy" ||
          init_mode == "pairs_xy_x" ||
          init_mode == "pairs_xy_y" ||
          init_mode == "pairs_xy_z"){
        dx[2] = 0.;
      }
      else {
        dx[2] = rnd_normal(gen);
      }
      dx *= 0.5*ds/dx.norm();

      Vector3d x_a = x0_ + dx;
      intp.probe(x_a);
      bool inside_a = intp.inside_domain();
      Vector3d x_b = x0_ - dx;
      intp.probe(x_b);
      bool inside_b = intp.inside_domain();
      if (inside_a && inside_b){
        cout << "INSIDE" << endl;
        pos_init.push_back(x_a);
        pos_init.push_back(x_b);
        double ds0 = (x_a-x_b).norm();
        edges.push_back({{2*ipair, 2*ipair+1}, ds0});
        ++ipair;
      }
      else if (init_mode == "pair_xyz" ||
               init_mode == "pair_xy" ||
               init_mode == "pair_xz" ||
               init_mode == "pair_yz"){
        cout << "Pair not inside domain" << endl;
        exit(0);
      }
    }
    Nrw = 2*Npairs;
    return pos_init;
  }

  bool init_rand_x = false;
  bool init_rand_y = false;
  bool init_rand_z = false;
  if (init_mode == "line_x" ||
      init_mode == "plane_xy" ||
      init_mode == "plane_xz" ||
      init_mode == "volume"){
    init_rand_x = true;
  }
  if (init_mode == "line_y" ||
      init_mode == "plane_xy" ||
      init_mode == "plane_yz" ||
      init_mode == "volume"){
    init_rand_y = true;

  }
  if (init_mode == "line_z" ||
      init_mode == "plane_xz" ||
      init_mode == "plane_yz" ||
      init_mode == "volume"){
    init_rand_z = true;
  }

  double dx_est = pow(Lx*Ly*Lz/10000000, 1./3);
  if (init_rand_x) Nx = Lx/dx_est;
  if (init_rand_y) Ny = Ly/dx_est;
  if (init_rand_z) Nz = Lz/dx_est;
  double dx = Lx/Nx;
  double dy = Ly/Ny;
  double dz = Lz/Nz;

  double ww;

  vector<double> wei;
  vector<Vector3d> pos;
  for (Uint ix=0; ix<Nx; ++ix){
    for (Uint iy=0; iy<Ny; ++iy){
      for (Uint iz=0; iz<Nz; ++iz){
        x = x0;
        if (init_rand_x) x[0] = (ix+0.5)*dx;
        if (init_rand_y) x[1] = (iy+0.5)*dy;
        if (init_rand_z) x[2] = (iz+0.5)*dz;
        intp.probe(x);
        if (init_weight == "ux"){
          ww = abs(intp.get_ux());
        }
        else if (init_weight == "uy"){
          ww = abs(intp.get_uy());
        }
        else if (init_weight == "uz"){
          ww = abs(intp.get_uz());
        }
        else if (init_weight == "u"){
          ww = sqrt(pow(intp.get_ux(), 2)
                    + pow(intp.get_uy(), 2)
                    + pow(intp.get_uz(), 2));
        }
        else {
          ww = 1.;
        }
        wei.push_back(ww);
        pos.push_back(x);
      }
    }
  }
  uniform_real_distribution<> uni_dist_dx(-0.5*dx, 0.5*dx);
  uniform_real_distribution<> uni_dist_dy(-0.5*dy, 0.5*dy);
  uniform_real_distribution<> uni_dist_dz(-0.5*dz, 0.5*dz);
  discrete_distribution<Uint> discrete_dist(wei.begin(), wei.end());

  vector<Vector3d> pos_init;
  for (Uint irw=0; irw<Nrw; ++irw){
    do {
      Uint ind = discrete_dist(gen);
      x = pos[ind];
      if (init_rand_x) x[0] += uni_dist_dx(gen);
      if (init_rand_y) x[1] += uni_dist_dy(gen);
      if (init_rand_z) x[2] += uni_dist_dz(gen);
      intp.probe(x);
    } while (!intp.inside_domain());
    pos_init.push_back(x);
  }

  sort(pos_init.begin(), pos_init.end(), less_than_op());

  for (Uint irw=1; irw < Nrw; ++irw){
    double ds0 = dist(pos_init[irw-1], pos_init[irw]);
    if (ds0 < 2)  // 2 lattice units
      edges.push_back({{irw-1, irw}, ds0});
    // Needs customization for 2D/3D applications
  }

  return pos_init;
}

#endif
