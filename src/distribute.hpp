#ifndef __DISTRIBUTE_HPP
#define __DISTRIBUTE_HPP

#include <iomanip>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include "Interpol.hpp"
#include "utils.hpp"


struct less_than_op {
  inline bool operator() (const Vector3d &a, const Vector3d &b){
    return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]) || (a[0] == b[0] && a[1] == b[1] && a[2] < b[2]);
  }
};

void load_positions(std::string input_file,
                    std::vector<Vector3d> &pos_init){
  std::ifstream infile(input_file);
  double x, y, z;
  while (infile >> x >> y >> z){
    pos_init.push_back({x, y, z});
  }
  infile.close();
}

void load_positions(std::string input_file,
                    std::vector<Vector3d> &pos_init,
                    const Uint Nrw){
  std::ifstream infile(input_file);
  double x, y, z;
  while (infile >> x >> y >> z){
    pos_init.push_back({x, y, z});
  }
  infile.close();
  if (Nrw != pos_init.size()){
    std::cout << "Wrong dimensions..." << std::endl;
    exit(0);
  }
}

void dump_positions(std::string output_file,
                    std::vector<Vector3d>& x_rw,
                    const Uint Nrw){
  std::ofstream outfile(output_file);
  for (Uint irw=0; irw<Nrw; ++irw){
    outfile << std::setprecision(12)
            << x_rw[irw][0] << " "
            << x_rw[irw][1] << " "
            << x_rw[irw][2] << std::endl;
  }
  outfile.close();
}

void dump_positions(const std::string output_file,
                    const std::vector<Vector3d> &pos){
  std::ofstream outfile(output_file);
  for (std::vector<Vector3d>::const_iterator posit=pos.begin();
       posit != pos.end(); ++posit){
    outfile << std::setprecision(12)
            << (*posit)[0] << " "
            << (*posit)[1] << " "
            << (*posit)[2] << std::endl;
  }
  outfile.close();
}

void load_faces(std::string input_file,
                FacesType& faces){
  std::ifstream infile(input_file);
  Uint first, second, third;
  double fourth;
  while (infile >> first >> second >> third >> fourth){
    faces.push_back({{first, second, third}, fourth});
  }
  infile.close();
}

void load_edges(std::string input_file,
                EdgesType &edges){
  std::ifstream infile(input_file);
  Uint first, second;
  double third;
  while (infile >> first >> second >> third){
    edges.push_back({{first, second}, third});
  }
  infile.close();
}

void load_list(std::string input_file,
               std::list<Uint> &li){
  std::ifstream infile(input_file);
  Uint a;
  while (infile >> a){
    li.push_back(a);
  }
  infile.close();
}

void dump_list(std::string output_file,
               const std::list<Uint> &li){
  std::ofstream outfile(output_file);
  for (std::list<Uint>::const_iterator lit=li.begin();
       lit != li.end(); ++lit){
    outfile << *lit << std::endl;
  }
  outfile.close();
}

void dump_faces(std::string output_file,
                const FacesType &faces){
  std::ofstream outfile(output_file);
  for (FacesType::const_iterator faceit = faces.begin();
       faceit != faces.end(); ++faceit){
    outfile << faceit->first[0] << " " << faceit->first[1] << " " << faceit->first[2]
            << " " << faceit->second << std::endl;
  }
  outfile.close();
}

void dump_edges(std::string output_file,
                EdgesType &edges){
  std::ofstream outfile(output_file);
  for (EdgesType::iterator edgeit = edges.begin();
       edgeit != edges.end(); ++edgeit){
    outfile << edgeit->first[0] << " " << edgeit->first[1] << " " << edgeit->second << std::endl;
  }
  outfile.close();
}

void load_colors(std::string input_file,
                 std::vector<double>& c_rw, const Uint Nrw){
  std::ifstream infile(input_file);
  for (Uint irw=0; irw < Nrw; ++irw){
    infile >> c_rw[irw];
  }
  infile.close();
}

void dump_colors(std::string output_file,
                 std::vector<double>& c_rw, const Uint Nrw){
  std::ofstream outfile(output_file);
  for (Uint irw=0; irw < Nrw; ++irw){
    outfile << std::setprecision(12) << c_rw[irw] << std::endl;
  }
  outfile.close();
}

std::vector<Vector3d> initial_positions(const std::string init_mode,
                                        const std::string init_weight,
                                        Uint &Nrw,
                                        const Vector3d &x0,
                                        const double La,
                                        const double Lb,
                                        const double ds,
                                        const double t0,
                                        Interpol *intp,
                                        std::mt19937 &gen,
                                        EdgesType &edges,
                                        FacesType &faces
                                   ){
  intp->update(t0);

  Vector3d x;
  double Lx = intp->get_Lx();
  double Ly = intp->get_Ly();
  double Lz = intp->get_Lz();
  Vector3d x_min = intp->get_x_min();
  Vector3d x_max = intp->get_x_max();

  Uint Nx = 1;
  Uint Ny = 1;
  Uint Nz = 1;

  std::vector<std::string> key = split_string(init_mode, "_");

  if (key.size() == 0){
    std::cout << "init_mode not specified." << std::endl;
    exit(0);
  }

  if (key[0] == "point"){
    std::vector<Vector3d> pos_init;
    intp->probe(x0);
    for (Uint irw=0; irw < Nrw; ++irw){
      if (intp->inside_domain()){
        pos_init.push_back(x0);
      }
    }
    edges.clear();
    return pos_init;
  }

  if (key[0] == "uniform"){
    std::vector<Vector3d> pos_init;

    Vector3d x_a = x0;
    Vector3d x_b = x0;
    if (key[1] == "x"){
      x_a[0] = x_min[0];
      x_b[0] = x_max[0];
    }
    else if (key[1] == "y"){
      x_a[1] = x_min[1];
      x_b[1] = x_max[1];
    }
    else if (key[1] == "z"){
      x_a[2] = x_min[2];
      x_b[2] = x_max[2];
    }
    else {
      std::cout << "Unrecognized initialization..." << std::endl;
      exit(0);
    }
    Vector3d Dx = (x_b - x_a) / (Nrw-1);
    for (Uint irw=0; irw < Nrw; ++irw){
      x = x_a + Dx * irw;
      intp->probe(x);
      if (intp->inside_domain()){
        pos_init.push_back(x);
      }
    }
    for (Uint irw=0; irw < pos_init.size()-1; ++irw){
      if ((pos_init[irw] - pos_init[irw+1]).norm() < 1.5*Dx.norm()){
        edges.push_back({{irw, irw+1}, dist(pos_init[irw], pos_init[irw+1])});
      }
    }
    return pos_init;
  }

  if (key[0] == "sheet"){
    std::vector<Vector3d> pos_init;

    Vector3d n(0., 0., 0.);
    Vector3d ta(0., 0., 0.);
    Vector3d tb(0., 0., 0.);
    if (key[1] == "xy"){
      n[2] = 1.0;
      ta[0] = 1.0;
      tb[1] = 1.0;
    }
    if (key[1] == "xz"){
      n[1] = 1.0;
      ta[0] = 1.0;
      tb[2] = 1.0;
    }
    if (key[1] == "yz"){
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

    intp->probe(x00);
    bool inside_00 = intp->inside_domain();
    intp->probe(x01);
    bool inside_01 = intp->inside_domain();
    intp->probe(x10);
    bool inside_10 = intp->inside_domain();
    intp->probe(x11);
    bool inside_11 = intp->inside_domain();
    if (inside_00 && inside_01 && inside_10 && inside_11){
      std::cout << "Sheet inside domain." << std::endl;
    }
    else {
      std::cout << "Sheet not inside domain" << std::endl;
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
  else if (key[0] == "pair" || key[0] == "pairs"){
    std::vector<Vector3d> pos_init;

    std::uniform_real_distribution<> uni_dist_x(x_min[0], x_max[0]);
    std::uniform_real_distribution<> uni_dist_y(x_min[1], x_max[1]);
    std::uniform_real_distribution<> uni_dist_z(x_min[2], x_max[2]);
    std::normal_distribution<double> rnd_normal(0.0, 1.0);

    //std::cout << "x_min = " << x_min << std::endl;
    //std::cout << "x_max = " << x_max << std::endl;

    // std::uniform_real_distribution<> uni_dist_theta(0., 2*M_PI);
    Uint Npairs = (key[0] == "pair") ? 1 : Nrw/2;

    std::cout << "Npairs = " << Npairs << std::endl;

    Vector3d x0_ = x0;
    Uint ipair=0;
    while (ipair < Npairs){
      if (key[0] == "pairs" && key.size() == 3){
        if (contains(key[2], "x")){
          x0_[0] = uni_dist_x(gen);
        }
        if (contains(key[2], "y")){
          x0_[1] = uni_dist_y(gen);
        }
        if (contains(key[2], "z")){
          x0_[2] = uni_dist_z(gen);
        }
      }
      Vector3d dx(0., 0., 0.);
      if (!contains(key[1], "x")){
        dx[0] = 0.;
      }
      else {
        dx[0] = rnd_normal(gen);
      }
      if (!contains(key[1], "y")){
        dx[1] = 0.;
      }
      else {
        dx[1] = rnd_normal(gen);
      }
      if (!contains(key[1], "z")){
        dx[2] = 0.;
      }
      else {
        dx[2] = rnd_normal(gen);
      }
      dx *= 0.5*ds/dx.norm();

      Vector3d x_a = x0_ + dx;
      intp->probe(x_a);
      bool inside_a = intp->inside_domain();
      Vector3d x_b = x0_ - dx;
      intp->probe(x_b);
      bool inside_b = intp->inside_domain();
      if (inside_a && inside_b){
        //std::cout << "INSIDE" << std::endl;
        // std::cout << "INSIDE: " << x_a << " " << x_b << std::endl; 
        pos_init.push_back(x_a);
        pos_init.push_back(x_b);
        double ds0 = (x_a-x_b).norm();
        edges.push_back({{2*ipair, 2*ipair+1}, ds0});
        ++ipair;
      }
      else if (key[0] == "pair"){
        std::cout << "Pair not inside domain" << std::endl;
        exit(0);
      }
      //std::cout << x_a << " " << x_b << std::endl;
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

  double tol = 1e-12;
  Uint N_est = 1000000;
  double dx_est;
  if (Lx > tol && Ly > tol && Lz > tol){
    dx_est = pow(Lx*Ly*Lz/N_est, 1./3);
  }
  else if (Lx > tol && Ly > tol){
    dx_est = pow(Lx*Ly/N_est, 1./2);
  }
  else if (Lx > tol && Lz > tol){
    dx_est = pow(Lx*Lz/N_est, 1./2);
  }
  else if (Ly > tol && Lz > tol){
    dx_est = pow(Ly*Lz/N_est, 1./2);
  }
  else if (Lx > tol){
    dx_est = Lx/N_est;
  }
  else if (Ly > tol){
    dx_est = Ly/N_est;
  }
  else if (Lz > tol){
    dx_est = Lz/N_est;
  }
  else {
    std::cout << "Something is wrong with the domain!" << std::endl;
    exit(0);
  }
  if (init_rand_x) Nx = Lx/dx_est+1;
  if (init_rand_y) Ny = Ly/dx_est+1;
  if (init_rand_z) Nz = Lz/dx_est+1;
  double dx = Lx/Nx;
  double dy = Ly/Ny;
  double dz = Lz/Nz;

  double ww;

  std::vector<double> wei;
  std::vector<Vector3d> pos;
  for (Uint ix=0; ix<Nx; ++ix){
    for (Uint iy=0; iy<Ny; ++iy){
      for (Uint iz=0; iz<Nz; ++iz){
        x = x0;
        if (init_rand_x) x[0] = x_min[0]+(ix+0.5)*dx;
        if (init_rand_y) x[1] = x_min[1]+(iy+0.5)*dy;
        if (init_rand_z) x[2] = x_min[2]+(iz+0.5)*dz;
        intp->probe(x);
        if (init_weight == "ux"){
          ww = abs(intp->get_ux());
        }
        else if (init_weight == "uy"){
          ww = abs(intp->get_uy());
        }
        else if (init_weight == "uz"){
          ww = abs(intp->get_uz());
        }
        else if (init_weight == "u"){
          ww = sqrt(pow(intp->get_ux(), 2)
                    + pow(intp->get_uy(), 2)
                    + pow(intp->get_uz(), 2));
        }
        else {
          ww = 1.;
        }

        wei.push_back(ww);
        pos.push_back(x);
      }
    }
  }
  std::uniform_real_distribution<> uni_dist_dx(-0.5*dx, 0.5*dx);
  std::uniform_real_distribution<> uni_dist_dy(-0.5*dy, 0.5*dy);
  std::uniform_real_distribution<> uni_dist_dz(-0.5*dz, 0.5*dz);
  std::discrete_distribution<Uint> discrete_dist(wei.begin(), wei.end());

  std::vector<Vector3d> pos_init;
  for (Uint irw=0; irw<Nrw; ++irw){
    do {
      Uint ind = discrete_dist(gen);
      x = pos[ind];
      if (init_rand_x) x[0] += uni_dist_dx(gen);
      if (init_rand_y) x[1] += uni_dist_dy(gen);
      if (init_rand_z) x[2] += uni_dist_dz(gen);
      intp->probe(x);
      //std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
    } while (!intp->inside_domain());
    pos_init.push_back(x);
  }

  sort(pos_init.begin(), pos_init.end(), less_than_op());

  for (Uint irw=1; irw < Nrw; ++irw){
    double ds0 = dist(pos_init[irw-1], pos_init[irw]);
    if (ds0 < 10*ds)  // 2 lattice units (before) --> 10 x ds_max (now)
      edges.push_back({{irw-1, irw}, ds0});
    // Needs customization for 2D/3D applications
  }

  return pos_init;
}

#endif
