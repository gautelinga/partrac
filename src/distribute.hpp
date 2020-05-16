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

void load_positions(string input_file,
                    vector<array<double, 3>> &pos_init,
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
                    double* x_rw, double* y_rw, double* z_rw,
                    const Uint Nrw){
  ofstream outfile(output_file);
  for (Uint irw=0; irw<Nrw; ++irw){
    outfile << std::setprecision(12) << x_rw[irw] << " " << y_rw[irw] << " " << z_rw[irw] << endl;
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

vector<array<double, 3>> initial_positions(string init_mode,
                                           string init_weight,
                                           Uint &Nrw,
                                           double x0,
                                           double y0,
                                           double z0,
                                           double La,
                                           double Lb,
                                           const double ds,
                                           Interpol &intp,
                                           mt19937 &gen,
                                           EdgesType &edges,
                                           FacesType &faces
                                           ){
  double x, y, z;
  double Lx = intp.get_Lx();
  double Ly = intp.get_Ly();
  double Lz = intp.get_Lz();

  Uint Nx = 1;
  Uint Ny = 1;
  Uint Nz = 1;

  if (init_mode == "sheet_xy" ||
      init_mode == "sheet_xz" ||
      init_mode == "sheet_yz"){
    vector<array<double, 3>> pos_init;

    double n_x = 0.;
    double n_y = 0.;
    double n_z = 0.;
    double ta_x = 0.;
    double ta_y = 0.;
    double ta_z = 0.;
    double tb_x = 0.;
    double tb_y = 0.;
    double tb_z = 0.;
    if (init_mode == "sheet_xy"){
      n_z = 1.0;
      ta_x = 1.0;
      tb_y = 1.0;
    }
    if (init_mode == "sheet_xz"){
      n_y = 1.0;
      ta_x = 1.0;
      tb_z = 1.0;
    }
    if (init_mode == "sheet_yz"){
      n_x = 1.0;
      ta_y = 1.0;
      tb_z = 1.0;
    }

    double x00 = x0 - La/2*ta_x - Lb/2*tb_x;
    double y00 = y0 - La/2*ta_y - Lb/2*tb_y;
    double z00 = z0 - La/2*ta_z - Lb/2*tb_z;

    double x01 = x0 + La/2*ta_x + Lb/2*tb_x;
    double y01 = y0 + La/2*ta_y - Lb/2*tb_y;
    double z01 = z0 - La/2*ta_z - Lb/2*tb_z;

    double x10 = x0 + La/2*ta_x + Lb/2*tb_x;
    double y10 = y0 + La/2*ta_y + Lb/2*tb_y;
    double z10 = z0 + La/2*ta_z + Lb/2*tb_z;

    double x11 = x0 - La/2*ta_x - Lb/2*tb_x;
    double y11 = y0 - La/2*ta_y + Lb/2*tb_y;
    double z11 = z0 - La/2*ta_z + Lb/2*tb_z;

    if (false){
      cout << n_x << " " << n_y << " " << n_z << endl;
      cout << x00 << " " << y00 << " " << z00 << endl;
      cout << x01 << " " << y01 << " " << z01 << endl;
      cout << x10 << " " << y10 << " " << z10 << endl;
      cout << x11 << " " << y11 << " " << z11 << endl;
    }

    intp.probe(x00, y00, z00);
    bool inside_00 = intp.inside_domain();
    intp.probe(x01, y01, z01);
    bool inside_01 = intp.inside_domain();
    intp.probe(x10, y10, z10);
    bool inside_10 = intp.inside_domain();
    intp.probe(x11, y11, z11);
    bool inside_11 = intp.inside_domain();
    if (inside_00 && inside_01 && inside_10 && inside_11){
      cout << "Sheet inside domain." << endl;
    }
    else {
      cout << "Sheet not inside domain" << endl;
      exit(0);
    }

    pos_init.push_back({x00, y00, z00});
    pos_init.push_back({x01, y01, z01});
    pos_init.push_back({x10, y10, z10});
    pos_init.push_back({x11, y11, z11});

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
    vector<array<double, 3>> pos_init;

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

    Uint ipair=0;
    while (ipair < Npairs){
      if (init_mode == "pairs_xyz_x" ||
          init_mode == "pairs_xy_x" ||
          init_mode == "pairs_xz_x" ||
          init_mode == "pairs_yz_x"){
        x0 = uni_dist_x(gen);
      }
      if (init_mode == "pairs_xyz_y" ||
          init_mode == "pairs_xy_y" ||
          init_mode == "pairs_xz_y" ||
          init_mode == "pairs_yz_y"){
        y0 = uni_dist_y(gen);
      }
      if (init_mode == "pairs_xyz_z" ||
          init_mode == "pairs_xy_z" ||
          init_mode == "pairs_xz_z" ||
          init_mode == "pairs_yz_z"){
        z0 = uni_dist_z(gen);
      }

      double dx, dy, dz;
      if (init_mode == "pair_yz" ||
          init_mode == "pairs_yz_x" ||
          init_mode == "pairs_yz_y" ||
          init_mode == "pairs_yz_z"){
        dx = 0.;
      }
      else {
        dx = rnd_normal(gen);
      }
      if (init_mode == "pair_xz" ||
          init_mode == "pairs_xz_x" ||
          init_mode == "pairs_xz_y" ||
          init_mode == "pairs_xz_z"){
        dy = 0.;
      }
      else {
        dy = rnd_normal(gen);
      }
      if (init_mode == "pair_xy" ||
          init_mode == "pairs_xy_x" ||
          init_mode == "pairs_xy_y" ||
          init_mode == "pairs_xy_z"){
        dz = 0.;
      }
      else {
        dz = rnd_normal(gen);
      }
      double ds_tentative = sqrt(pow(dx, 2)+pow(dy, 2)+pow(dz, 2));
      dx *= 0.5*ds/ds_tentative;
      dy *= 0.5*ds/ds_tentative;
      dz *= 0.5*ds/ds_tentative;

      double x_a = x0 + dx;
      double y_a = y0 + dy;
      double z_a = z0 + dz;
      intp.probe(x_a, y_a, z_a);
      bool inside_a = intp.inside_domain();
      double x_b = x0 - dx;
      double y_b = y0 - dy;
      double z_b = z0 - dz;
      intp.probe(x_b, y_b, z_b);
      bool inside_b = intp.inside_domain();
      if (inside_a && inside_b){
        cout << "INSIDE" << endl;
        pos_init.push_back({x_a, y_a, z_a});
        pos_init.push_back({x_b, y_b, z_b});
        double ds0 = sqrt(pow(x_a-x_b, 2) + pow(y_a-y_b, 2) + pow(z_a-z_b, 2));
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
  vector<array<double, 3>> pos;
  for (Uint ix=0; ix<Nx; ++ix){
    for (Uint iy=0; iy<Ny; ++iy){
      for (Uint iz=0; iz<Nz; ++iz){
        x = x0;
        y = y0;
        z = z0;
        if (init_rand_x) x = (ix+0.5)*dx;
        if (init_rand_y) y = (iy+0.5)*dy;
        if (init_rand_z) z = (iz+0.5)*dz;
        intp.probe(x, y, z);
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
        pos.push_back({x, y, z});
      }
    }
  }
  uniform_real_distribution<> uni_dist_dx(-0.5*dx, 0.5*dx);
  uniform_real_distribution<> uni_dist_dy(-0.5*dy, 0.5*dy);
  uniform_real_distribution<> uni_dist_dz(-0.5*dz, 0.5*dz);
  discrete_distribution<Uint> discrete_dist(wei.begin(), wei.end());

  vector<array<double, 3>> pos_init;
  for (Uint irw=0; irw<Nrw; ++irw){
    do {
      Uint ind = discrete_dist(gen);
      x = pos[ind][0];
      y = pos[ind][1];
      z = pos[ind][2];
      if (init_rand_x) x += uni_dist_dx(gen);
      if (init_rand_y) y += uni_dist_dy(gen);
      if (init_rand_z) z += uni_dist_dz(gen);
      intp.probe(x, y, z);
    } while (!intp.inside_domain());
    pos_init.push_back({x, y, z});
  }

  sort(pos_init.begin(), pos_init.end());

  for (Uint irw=1; irw < Nrw; ++irw){
    double ds0 = dist(pos_init[irw-1], pos_init[irw]);
    if (ds0 < 2)  // 2 lattice units
      edges.push_back({{irw-1, irw}, ds0});
    // Needs customization for 2D/3D applications
  }

  return pos_init;
}

#endif
