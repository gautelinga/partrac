#include <iostream>
#include "io.hpp"
#include "utils.hpp"
#include "Interpol.hpp"
#include "distribute.hpp"
#include <vector>
#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <random>
#include "H5Cpp.h"

using namespace H5;

int main(int argc, char* argv[]){
  // Input parameters
  if (argc < 2) {
    std::cout << "Specify an input timestamps file." << std::endl;
    return 0;
  }
  // Default parameters
  double Dm = 1.0;
  double t0 = 0.;
  double T = 10000000.;
  int Nrw = 100;
  double traj_intv = 1.0;
  double pos_intv = 10.0;
  double stat_intv = 10.0;
  bool dump_traj = true;
  bool verbose = false;
  double U0 = 1.0;
  double x0 = 0.0;
  double y0 = 0.0;
  double z0 = 0.0;
  double dt = 1.0;
  int int_order = 1;
  string init_mode = "line_x";  // what else?
  string init_weight = "none";

  size_t found;
  string argstr, key, val;
  std::cout << "TEST" << std::endl;
  for (int iarg=2; iarg < argc; ++iarg){
    argstr = argv[iarg];
    found = argstr.find('=');
    if (found != string::npos){
      key = argstr.substr(0, found);
      val = argstr.substr(found+1);
      boost::trim(key);
      boost::trim(val);

      if (key == "Dm") Dm = stod(val);
      if (key == "t0") t0 = stod(val);
      if (key == "T") T = stod(val);
      if (key == "dt") dt = stod(val);
      if (key == "Nrw") Nrw = stoi(val);
      if (key == "traj_intv") traj_intv = stod(val);
      if (key == "pos_intv") pos_intv = stod(val);
      if (key == "stat_intv") stat_intv = stod(val);
      if (key == "dump_traj") dump_traj = (val == "true" || val == "True") ? true : false;
      if (key == "verbose") verbose = (val == "true" || val == "True") ? true : false;
      if (key == "U") U0 = stod(val);
      if (key == "x0") x0 = stod(val);
      if (key == "y0") y0 = stod(val);
      if (key == "z0") z0 = stod(val);
      if (key == "int_order") int_order = stoi(val);
      if (key == "init_mode") init_mode = val;
      if (key == "init_weight") init_weight = val;
    }
  }

  if (verbose){
    print_param("Dm         ", Dm);
    print_param("t0         ", t0);
    print_param("T          ", T);
    print_param("dt         ", dt);
    print_param("Nrw        ", Nrw);
    print_param("traj_intv  ", traj_intv);
    print_param("pos_intv   ", pos_intv);
    print_param("stat_intv  ", stat_intv);
    print_param("dump_traj  ", dump_traj ? "true" : "false");
    print_param("verbose    ", verbose ? "true" : "false");
    print_param("U          ", U0);
    print_param("x0         ", x0);
    print_param("y0         ", y0);
    print_param("z0         ", z0);
    print_param("int_order  ", int_order);
    print_param("init_mode  ", init_mode);
    print_param("init_weight", init_mode);
  }
  traj_intv = max(traj_intv, dt);
  pos_intv = max(pos_intv, dt);
  stat_intv = max(pos_intv, dt);

  string infilename = string(argv[1]);

  Interpol intp(infilename);
  string folder = intp.get_folder();
  
  string rwfolder = create_folder(folder + "/RandomWalkers/");
  string newfolder = create_folder(rwfolder +
    "/Dm" + to_string(Dm) + "_U" + to_string(U0) +
    "_dt" + to_string(dt) + "_Nrw" + to_string(Nrw) + "/");
  string trajfolder = create_folder(newfolder + "Trajectories/");
  string posfolder = create_folder(newfolder + "Positions/");

  random_device rd;
  mt19937 gen(rd());
  
  int n = 0;
  double x, y, z;
  double ux, uy, uz;
  double rho, p;
  double divu, vortz;
  
  double Lx = intp.get_Lx();
  double Ly = intp.get_Ly();
  double Lz = intp.get_Lz();
  double U02 = U0*U0;
  
  t0 = max(intp.get_t_min(), t0);  // 587500.0;
  T = min(intp.get_t_max(), T);
  
  cout << "t0 = " << t0 << endl;
  cout << "T  = " << T  << endl;
  
  intp.update(587500.0);
  uniform_real_distribution<> uni_dist_x(0*Lx, Lx);
  uniform_real_distribution<> uni_dist_y(0*Ly, Ly);
  uniform_real_distribution<> uni_dist_z(0*Lz, Lz);
  
  ofstream ofile("test.dat");
  while (n < 100000){
    x = uni_dist_x(gen);  // mod?
    y = uni_dist_y(gen);
    z = uni_dist_z(gen);

    intp.probe(x, y, z);
    if (intp.inside_domain()){
      ux = intp.get_ux();
      uy = intp.get_uy();
      uz = intp.get_uz();
      rho = intp.get_rho();
      p = intp.get_p();
      divu = intp.get_divu();
      vortz = intp.get_vortz();
      ofile << x  << " " << y  << " " << z  << " "
	    << ux << " " << uy << " " << uz << " "
	    << rho << " " << p << " " << divu << " "
	    << vortz << endl;
    }    
    ++n;
  }
  ofile.close();
  
  normal_distribution<double> rnd_normal(0.0, 1.0);

  double* x_rw = new double[Nrw];
  double* y_rw = new double[Nrw];
  double* z_rw = new double[Nrw];
  double* ux_rw = new double[Nrw];
  double* uy_rw = new double[Nrw];
  double* uz_rw = new double[Nrw];
  // Second-order terms
  double* Jux_rw = new double[Nrw];
  double* Juy_rw = new double[Nrw];
  double* Juz_rw = new double[Nrw];
  if (int_order > 2){
    cout << "No support for such high temporal integration order." << endl;
    exit(0);
  }

  vector<array<double, 3>> pos_init = initial_positions(init_mode,
							init_weight,
							Nrw,
							x0, y0, z0,
							intp, gen);
  
  ofstream paramsfile(newfolder + "/params.dat");
  write_param(paramsfile, "Dm", Dm);
  write_param(paramsfile, "t0", t0);
  write_param(paramsfile, "T", T);
  write_param(paramsfile, "dt", dt);
  write_param(paramsfile, "Nrw", Nrw);
  write_param(paramsfile, "traj_intv", traj_intv);
  write_param(paramsfile, "pos_intv", pos_intv);
  write_param(paramsfile, "stat_intv", stat_intv);
  write_param(paramsfile, "dump_traj", dump_traj ? "true" : "false");
  write_param(paramsfile, "verbose", verbose ? "true" : "false");
  write_param(paramsfile, "U", U0);
  write_param(paramsfile, "x0", x0);
  write_param(paramsfile, "y0", y0);
  write_param(paramsfile, "z0", z0);
  write_param(paramsfile, "Lx", Lx);
  write_param(paramsfile, "Ly", Ly);
  write_param(paramsfile, "Lz", Lz);
  write_param(paramsfile, "int_order", int_order);
  write_param(paramsfile, "init_mode", init_mode);
  write_param(paramsfile, "init_weight", init_mode);
  paramsfile.close();
  
  ofstream* traj_outs = new ofstream[Nrw];
  for (int irw=0; irw < Nrw; ++irw){
    // Initial position
    if (dump_traj){
      string filename = trajfolder + "traj_" + to_string(irw) + ".traj";
      traj_outs[irw].open(filename, ofstream::out);
      if (verbose)
	std::cout << "Opened: " << irw << std::endl;
    }
    x = pos_init[irw][0];
    y = pos_init[irw][1];
    z = pos_init[irw][2];

    // cout << x << " " << y << " " << z << endl;

    x_rw[irw] = x;
    y_rw[irw] = y;
    z_rw[irw] = z;

    ux_rw[irw] = U0*intp.get_ux();
    uy_rw[irw] = U0*intp.get_uy();
    uz_rw[irw] = U0*intp.get_uz();
    
    // Second-order terms
    if (int_order >= 2){
      Jux_rw[irw] = U02*intp.get_Jux();
      Juy_rw[irw] = U02*intp.get_Juy();
      Juz_rw[irw] = U02*intp.get_Juz();
    }
  }

  double eta_x, eta_y, eta_z;
  double dx_rw, dy_rw, dz_rw;
  int it = 0;
  double dt2 = dt*dt;


  double sqrt2Dmdt = sqrt(2*Dm*dt);
  if (verbose){
    print_param("sqrt(2*Dm*dt)", sqrt2Dmdt);
    print_param("U*dt", U0*dt);
  }
  
  double x_mean, dx2_mean;
  double y_mean, dy2_mean;
  double z_mean, dz2_mean;
  double t = t0;

  int n_accepted = 0;
  int n_declined = 0;

  H5File posf_h5(newfolder + "/positions.h5", H5F_ACC_TRUNC);
  hsize_t posf_dims[2];
  posf_dims[0] = Nrw;
  posf_dims[1] = 6;
  DataSpace posf_dspace(2, posf_dims);
  double** posdata = new double*[Nrw];
  for (int irw=0; irw<Nrw; ++irw){
    posdata[irw] = new double[6];
  }
  
  ofstream statfile(newfolder + "/tdata.dat");
  ofstream declinedfile(newfolder + "/declinedpos.dat");
  while (t <= T){
    intp.update(t);
    if (it % int(stat_intv/dt) == 0){
      cout << "Time = " << t << endl;
    }
    //
    if (it % int(pos_intv/dt) == 0){
      string posfile = posfolder + "xy_t" + to_string(t) + ".pos";
      ofstream pos_out(posfile);

      for (int irw=0; irw < Nrw; ++irw){
	pos_out << irw << " "
		<< x_rw[irw] << " "
		<< y_rw[irw] << " "
		<< z_rw[irw] << " "
		<< ux_rw[irw] << " "
		<< uy_rw[irw] << " "
		<< uz_rw[irw] << std::endl;
	posdata[irw][0] = x_rw[irw];
	posdata[irw][1] = y_rw[irw];
	posdata[irw][2] = z_rw[irw];
	posdata[irw][3] = ux_rw[irw];
	posdata[irw][4] = uy_rw[irw];
	posdata[irw][5] = uz_rw[irw];
      }
      pos_out.close();

      DataSet posf_dset = posf_h5.createDataSet(to_string(t),
						PredType::NATIVE_DOUBLE,
						posf_dspace
						);
      posf_dset.write(posdata, PredType::NATIVE_DOUBLE);
      // posf_h5.close();
    }    
    if (it % int(stat_intv/dt) == 0){
      x_mean = 0.;
      y_mean = 0.;
      z_mean = 0.;
      dx2_mean = 0.;
      dy2_mean = 0.;
      dz2_mean = 0.;
      for (int irw=0; irw < Nrw; ++irw){
	// Sample mean
	x_mean += x_rw[irw]/Nrw;
	y_mean += y_rw[irw]/Nrw;
	z_mean += z_rw[irw]/Nrw;
      }
      for (int irw=0; irw < Nrw; ++irw){
	// Sample variance
	dx2_mean += pow(x_rw[irw]-x_mean, 2)/(Nrw-1);
	dy2_mean += pow(y_rw[irw]-y_mean, 2)/(Nrw-1);
	dz2_mean += pow(z_rw[irw]-z_mean, 2)/(Nrw-1);
      }
      statfile << t << " "
	       << x_mean << " " << dx2_mean << " "
	       << y_mean << " " << dy2_mean << " "
	       << z_mean << " " << dz2_mean << " "
	       << n_accepted << " " << n_declined << std::endl;
    }
    for (int irw=0; irw < Nrw; ++irw){
      if (dump_traj && it % int(traj_intv/dt) == 0){
	traj_outs[irw] << t << " "
		       << x_rw[irw] << " "
		       << y_rw[irw] << " "
		       << z_rw[irw] << " "
		       << ux_rw[irw] << " "
		       << uy_rw[irw] << " "
		       << uz_rw[irw] << std::endl;
      }
      dx_rw = ux_rw[irw]*dt;
      dy_rw = uy_rw[irw]*dt;
      dz_rw = uz_rw[irw]*dt;

      // Second-order terms
      if (int_order >= 2){
	dx_rw += 0.5*Jux_rw[irw]*dt2;
	dy_rw += 0.5*Juy_rw[irw]*dt2;
	dz_rw += 0.5*Juz_rw[irw]*dt2;
      }
      if (Dm > 0.0){
	eta_x = rnd_normal(gen);
	eta_y = rnd_normal(gen);
	eta_z = rnd_normal(gen);
	dx_rw += sqrt2Dmdt*eta_x;
	dy_rw += sqrt2Dmdt*eta_y;
	dz_rw += sqrt2Dmdt*eta_z;
      }
      intp.probe(x_rw[irw]+dx_rw,
		 y_rw[irw]+dy_rw,
		 z_rw[irw]+dz_rw);
      if (intp.inside_domain()){
	x_rw[irw] += dx_rw;
	y_rw[irw] += dy_rw;
	z_rw[irw] += dz_rw;
	
	ux_rw[irw] = U0*intp.get_ux();
	uy_rw[irw] = U0*intp.get_uy();
	uz_rw[irw] = U0*intp.get_uz();

	// Second-order terms
	if (int_order >= 2){
	  Jux_rw[irw] = U02*intp.get_Jux();
	  Juy_rw[irw] = U02*intp.get_Juy();
	  Juz_rw[irw] = U02*intp.get_Juz();
	}
	
	n_accepted++;
      }
      else {
	n_declined++;
	declinedfile << t << " "
		     << x_rw[irw]+dx_rw << " "
		     << y_rw[irw]+dy_rw << " "
		     << z_rw[irw]+dz_rw << std::endl;
      }
    }
    t += dt;
    it += 1;
  }
  
  if (dump_traj){
    for (int irw=0; irw < Nrw; ++irw){
      traj_outs[irw].close();
    }
  }
  statfile.close();
  declinedfile.close();

  return 0;

}
