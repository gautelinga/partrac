#include <iomanip>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include "Interpol.hpp"

#ifndef __DISTRIBUTE_HPP
#define __DISTRIBUTE_HPP

using namespace std;

void load_positions(string input_file,
		    vector<array<double, 3>> &pos_init,
		    const int Nrw){
  ifstream infile(input_file);
  double x, y, z;
  while (infile >> x >> y >> z){
    pos_init.push_back({x, y, z});
  }
  infile.close();
  if (Nrw != int(pos_init.size())){
    cout << "Wrong dimensions..." << endl;
    exit(0);
  }
}

void dump_positions(string output_file,
		    double* x_rw, double* y_rw, double* z_rw,
		    int Nrw){
  ofstream outfile(output_file);
  for (int irw=0; irw<Nrw; ++irw){
    outfile << std::setprecision(12) << x_rw[irw] << " " << y_rw[irw] << " " << z_rw[irw] << endl;
  }
  outfile.close();
}

void load_edges(string input_file,
		vector<tuple<int, int, double>> &edges){
  ifstream infile(input_file);
  int first, second;
  double third;
  while (infile >> first >> second >> third){
    edges.push_back({first, second, third});
  }
  infile.close();
}

void dump_edges(string output_file,
		vector<tuple<int, int, double>> &edges){
  ofstream outfile(output_file);
  for (vector<tuple<int, int, double>>::iterator edgeit = edges.begin();
       edgeit != edges.end(); ++edgeit){
    outfile << get<0>(*edgeit) << " " << get<1>(*edgeit) << " " << get<2>(*edgeit) << endl;
  }
  outfile.close();
}

void load_colors(string input_file,
		 double* c_rw, const int Nrw){
  ifstream infile(input_file);
  for (int irw=0; irw < Nrw; ++irw){
    infile >> c_rw[irw];
  }
  infile.close();
}

void dump_colors(string output_file,
		 double* c_rw, const int Nrw){
  ofstream outfile(output_file);
  for (int irw=0; irw < Nrw; ++irw){
    outfile << setprecision(12) << c_rw[irw] << endl;
  }
  outfile.close();
}

vector<array<double, 3>> initial_positions(string init_mode,
					   string init_weight,
					   int &Nrw,
					   const double x0,
					   const double y0,
					   const double z0,
					   const double ds,
					   Interpol &intp,
					   mt19937 &gen
					   ){
  double x, y, z;
  double Lx = intp.get_Lx();
  double Ly = intp.get_Ly();
  double Lz = intp.get_Lz();

  int Nx = 1;
  int Ny = 1;
  int Nz = 1;

  if (init_mode == "pair_xy"){
    uniform_real_distribution<> uni_dist_theta(0., 2*M_PI);
    double theta = uni_dist_theta(gen);
    double dx = 0.5*ds*cos(theta);
    double dy = 0.5*ds*sin(theta);

    vector<array<double, 3>> pos_init;
    double x_a = x0 + dx;
    double y_a = y0 + dy;
    double z_a = z0;
    intp.probe(x_a, y_a, z_a);
    bool inside_a = intp.inside_domain();
    double x_b = x0 - dx;
    double y_b = y0 - dy;
    double z_b = z0;
    bool inside_b = intp.inside_domain();
    if (!inside_a || !inside_b){
      cout << "Pair not inside domain" << endl;
      exit(0);
    }
    pos_init.push_back({x_a, y_a, z_a});
    pos_init.push_back({x_b, y_b, z_b});
    Nrw = 2;
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
  //double x, y, z;
  double ww;

  vector<double> wei;
  vector<array<double, 3>> pos;
  for (int ix=0; ix<Nx; ++ix){
    for (int iy=0; iy<Ny; ++iy){
      for (int iz=0; iz<Nz; ++iz){
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
  discrete_distribution<int> dist(wei.begin(), wei.end());

  vector<array<double, 3>> pos_init;
  for (int irw=0; irw<Nrw; ++irw){
    do {
      int ind = dist(gen);
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
  return pos_init;
}

#endif
