#include <iostream>
#include "io.hpp"
#include "utils.hpp"
#include "Parameters.hpp"
#include "Interpol.hpp"
#include "distribute.hpp"
#include <vector>
#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <sstream>
#include <random>
#include <math.h>
#include "H5Cpp.h"

using namespace H5;

double dist(int i1, int i2, double* x_rw, double* y_rw, double* z_rw){
  return sqrt(pow(x_rw[i1]-x_rw[i2], 2)
	      + pow(y_rw[i1]-y_rw[i2], 2)
	      + pow(z_rw[i1]-z_rw[i2], 2));
}

int main(int argc, char* argv[]){
  // Input parameters
  if (argc < 2) {
    std::cout << "Specify an input timestamps file." << std::endl;
    return 0;
  }
  Parameters prm(argc, argv);
  if (prm.restart_folder != ""){
    prm.parse_file(prm.restart_folder + "/Checkpoints/params.dat");
    prm.parse_cmd(argc, argv);
  }
  
  string infilename = string(argv[1]);

  Interpol intp(infilename);
  
  double Dm = prm.Dm;
  double dt = prm.dt;
  int Nrw = prm.Nrw;
  double U0 = prm.U0;

  bool refine = prm.refine;
  double ds_max = prm.ds_max;
  int Nrw_max = prm.Nrw_max;

  string folder = intp.get_folder();
  string rwfolder = create_folder(folder + "/RandomWalkers/");
  string newfolder;
  if (prm.restart_folder != ""){
    newfolder = prm.folder;
  }
  else {
    ostringstream ss_Dm, ss_dt, ss_Nrw;
    ss_Dm << std::scientific << std::setprecision(7) << Dm;
    ss_dt << std::scientific << std::setprecision(7) << dt;
    ss_Nrw << Nrw;
    newfolder = create_folder(rwfolder +
			      "/Dm" + ss_Dm.str() + // "_U" + to_string(prm.U0) +
			      "_dt" + ss_dt.str() +
			      "_Nrw" + ss_Nrw.str() + "/");
  }
  string posfolder = create_folder(newfolder + "Positions/");
  string checkpointsfolder = create_folder(newfolder + "Checkpoints/");
  string histfolder = create_folder(newfolder + "Histograms/");
  prm.folder = newfolder;

  prm.print();
  
  random_device rd;
  mt19937 gen(rd());
  
  double Lx = intp.get_Lx();
  double Ly = intp.get_Ly();
  double Lz = intp.get_Lz();
  prm.Lx = Lx;
  prm.Ly = Ly;
  prm.Lz = Lz;
  prm.nx = intp.get_nx();
  prm.ny = intp.get_ny();
  prm.nz = intp.get_nz();
  
  double U02 = U0*U0;
  
  double t0 = max(intp.get_t_min(), prm.t0);
  double T = min(intp.get_t_max(), prm.T);
  prm.t0 = t0;
  prm.T = T;
  
  if (prm.interpolation_test > 0){
    int n = 0;
    double x, y, z;
    double ux, uy, uz;
    double rho, p;
    double divu, vortz;
    
    intp.update(t0);
    uniform_real_distribution<> uni_dist_x(0, Lx);
    uniform_real_distribution<> uni_dist_y(0, Ly);
    uniform_real_distribution<> uni_dist_z(0, Lz);

    ofstream nodalfile(newfolder + "/nodal_values.dat");
    for (int ix=0; ix<intp.get_nx(); ++ix){
      for (int iy=0; iy<intp.get_ny(); ++iy){
	for (int iz=0; iz<intp.get_nz(); ++iz){
	  bool inside = intp.get_nodal_inside(ix, iy, iz);
	  ux = intp.get_nodal_ux(ix, iy, iz);
	  uy = intp.get_nodal_uy(ix, iy, iz);
	  uz = intp.get_nodal_uz(ix, iy, iz);
	  nodalfile << ix << " " << iy << " " << iz << " " << inside << " "
		    << ux << " " << uy << " " << uz << endl;
	}
      }
    }
    nodalfile.close();
    
    ofstream ofile(newfolder + "/interpolation.dat");
    while (n < prm.interpolation_test){
      x = uni_dist_x(gen);
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
  }
  
  normal_distribution<double> rnd_normal(0.0, 1.0);

  double* x_rw = new double[Nrw_max];
  double* y_rw = new double[Nrw_max];
  double* z_rw = new double[Nrw_max];

  double* c_rw = new double[Nrw_max];
  
  double* ux_rw = new double[Nrw_max];
  double* uy_rw = new double[Nrw_max];
  double* uz_rw = new double[Nrw_max];
  double* rho_rw = new double[Nrw_max];
  double* p_rw = new double[Nrw_max];
  
  // Second-order terms
  double* Jux_rw = new double[Nrw_max];
  double* Juy_rw = new double[Nrw_max];
  double* Juz_rw = new double[Nrw_max];
  if (prm.int_order > 2){
    cout << "No support for such high temporal integration order." << endl;
    exit(0);
  }

  vector<array<double, 3>> pos_init;
  vector<tuple<int, int, double>> edges;
  
  if (prm.restart_folder != ""){
    string posfile = prm.restart_folder + "/Checkpoints/positions.pos";
    load_positions(posfile, pos_init, Nrw);
    string edgefile = prm.restart_folder + "/Checkpoints/edges.edge";
    load_edges(edgefile, edges);
    string colfile = prm.restart_folder + "/Checkpoints/colors.col";
    load_colors(colfile, c_rw, Nrw);
  }
  else {
    pos_init = initial_positions(prm.init_mode,
				 prm.init_weight,
				 Nrw,
				 prm.x0, prm.y0, prm.z0,
				 prm.ds_max,
				 intp, gen);
  }

  double ds0, ds;
  
  for (int irw=0; irw < Nrw; ++irw){
    // Assign initial position
    x_rw[irw] = pos_init[irw][0];
    y_rw[irw] = pos_init[irw][1];
    z_rw[irw] = pos_init[irw][2];

    if (prm.restart_folder == ""){
      c_rw[irw] = double(irw)/(Nrw-1);
      if (irw > 0){
	ds0 = dist(irw-1, irw, x_rw, y_rw, z_rw);
	if (ds0 < 2)  // 2 lattice units
	  edges.push_back({irw-1, irw, ds0});
	// Needs customization for 2D/3D applications
      }
    }
    intp.update(t0);
    intp.probe(x_rw[irw],
	       y_rw[irw],
	       z_rw[irw]);
    ux_rw[irw] = U0*intp.get_ux();
    uy_rw[irw] = U0*intp.get_uy();
    uz_rw[irw] = U0*intp.get_uz();

    rho_rw[irw] = intp.get_rho();
    p_rw[irw] = intp.get_p();
    
    // Second-order terms
    if (prm.int_order >= 2){
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
  if (prm.verbose){
    print_param("sqrt(2*Dm*dt)", sqrt2Dmdt);
    print_param("U*dt", U0*dt);
  }
  
  double x_mean, dx2_mean;
  double y_mean, dy2_mean;
  double z_mean, dz2_mean;
  double t = t0;
  if (prm.restart_folder != ""){
    t = prm.t;
  }
  
  prm.dump(newfolder, t);
  
  unsigned long int n_accepted = prm.n_accepted;
  unsigned long int n_declined = prm.n_declined;
  
  const int NCOLS = 9;

  H5File* posf_h5 = new H5File(newfolder + "/data_from_t" + to_string(t) + ".h5", H5F_ACC_TRUNC);
  
  int int_stat_intv = int(prm.stat_intv/dt);
  int int_dump_intv = int(prm.dump_intv/dt);
  int int_checkpoint_intv = int(prm.checkpoint_intv/dt);
  int int_chunk_intv = int_dump_intv*prm.dump_chunk_size;
  int int_refine_intv = int(prm.refine_intv/dt);
  int int_hist_intv = int_stat_intv*prm.hist_chunk_size;
  
  double s;
  int n_too_long;
  double logdens_mean;
  double logdens_var;
  double logdens;
  
  string write_mode = prm.write_mode;
  
  ofstream statfile(newfolder + "/tdata_from_t" + to_string(t) + ".dat");
  ofstream declinedfile(newfolder + "/declinedpos_from_t" + to_string(t) + ".dat");
  while (t <= T){
    intp.update(t);
    if (it % int_stat_intv == 0){
      cout << "Time = " << t << endl;

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
      s = 0.;
      n_too_long = 0;
      double logdens_wsum = 0.;
      double s0 = 0.;
      vector<array<double,2>> logdens_vec;
      for (vector<tuple<int, int, double>>::iterator edgeit = edges.begin();
	   edgeit != edges.end(); ++edgeit){
	int first = get<0>(*edgeit);
	int second = get<1>(*edgeit);
	double ds0 = get<2>(*edgeit);
	ds = dist(first, second, x_rw, y_rw, z_rw);
	if (ds > ds_max){
	  ++n_too_long;
	}
	logdens = log(ds/ds0);
	logdens_wsum += logdens*ds;
	logdens_vec.push_back({logdens, ds});
	s += ds;
	s0 += ds0;
      }
      logdens_mean = logdens_wsum/s;
      double logdens_var_wsum = 0.;
      for (vector<array<double, 2>>::iterator lit = logdens_vec.begin();
	   lit != logdens_vec.end(); ++lit){
	logdens_var_wsum += pow((*lit)[0]-logdens_mean, 2)*(*lit)[1];
      }
      logdens_var = logdens_var_wsum/s;
      double s_ = 0.;
      if (int_hist_intv > 0 && it % int_hist_intv == 0){
	sort(logdens_vec.begin(), logdens_vec.end());
	ofstream logdensfile(histfolder + "/logdens_t" + to_string(t) + ".hist");
	for (vector<array<double, 2>>::iterator lit = logdens_vec.begin();
	     lit != logdens_vec.end(); ++lit){
	  s_ += (*lit)[1];
	  logdensfile << (*lit)[0] << " " << (*lit)[1] << " " << s_ << " " << s_/s << endl;
	}
	logdensfile.close();
      }
      
      statfile << t << " "  // 1
	       << x_mean << " "  // 2
	       << dx2_mean << " "
	       << y_mean << " "
	       << dy2_mean << " "
	       << z_mean << " "
	       << dz2_mean << " "
	       << n_accepted << " "
	       << n_declined << " "
	       << s << " "
	       << n_too_long << " "
	       << Nrw << " "
	       << logdens_mean << " "
	       << logdens_var << " "
	       << s0
	       << std::endl;
    }
    //
    if (it % int_checkpoint_intv == 0 || t+dt > T){
       prm.t = t;
       prm.n_accepted = n_accepted;
       prm.n_declined = n_declined;
       prm.Nrw = Nrw;
       prm.dump(checkpointsfolder);
       dump_positions(checkpointsfolder + "/positions.pos",
		      x_rw, y_rw, z_rw, Nrw);
       dump_edges(checkpointsfolder + "/edges.edge",
		  edges);
       dump_colors(checkpointsfolder + "/colors.col",
		  c_rw, Nrw);
    }
    // Refinement
    if (refine && it % int_refine_intv == 0){
      for (vector<tuple<int, int, double>>::iterator edgeit = edges.begin();
	   edgeit != edges.end(); ++edgeit){
	int first = get<0>(*edgeit);
	int second = get<1>(*edgeit);
	double ds0 = get<2>(*edgeit);
	ds = dist(first, second, x_rw, y_rw, z_rw);
	if (ds > ds_max && Nrw < Nrw_max){
	  get<1>(*edgeit) = Nrw;
	  get<2>(*edgeit) = ds0/2;
	  edgeit = edges.insert(edgeit+1, {Nrw, second, ds0/2});
	  edgeit -= 2;

	  x_rw[Nrw] = 0.5*(x_rw[first]+x_rw[second]);
	  y_rw[Nrw] = 0.5*(y_rw[first]+y_rw[second]);
	  z_rw[Nrw] = 0.5*(z_rw[first]+z_rw[second]);

	  c_rw[Nrw] = 0.5*(c_rw[first]+c_rw[second]);
	  
	  ux_rw[Nrw] = 0.5*(ux_rw[first]+ux_rw[second]);
	  uy_rw[Nrw] = 0.5*(uy_rw[first]+uy_rw[second]);
	  uz_rw[Nrw] = 0.5*(uz_rw[first]+uz_rw[second]);
	  if (it % int_dump_intv == 0){
	    rho_rw[Nrw] = 0.5*(rho_rw[first]+rho_rw[second]);
	    p_rw[Nrw] = 0.5*(p_rw[first]+p_rw[second]);
	  }
	  ++Nrw;
	}
      }
    }
    //
    if (it % int_dump_intv == 0){
      if (write_mode == "text"){
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
	}
	pos_out.close();
      }
      else if (write_mode == "hdf5"){
	// Clear file if it exists, otherwise create
	if (it % int_chunk_intv == 0 && it > 0){
	  posf_h5->close();
	  posf_h5 = new H5File(newfolder + "/data_from_t" + to_string(t) + ".h5", H5F_ACC_TRUNC);
	}
	hsize_t posf_dims[2];
	posf_dims[0] = Nrw;
	posf_dims[1] = NCOLS;
	DataSpace posf_dspace(2, posf_dims);
	double* posdata = new double[Nrw*NCOLS];
	for (int irw=0; irw < Nrw; ++irw){
	  posdata[irw*NCOLS+0] = x_rw[irw];
	  posdata[irw*NCOLS+1] = y_rw[irw];
	  posdata[irw*NCOLS+2] = z_rw[irw];
	  posdata[irw*NCOLS+3] = ux_rw[irw];
	  posdata[irw*NCOLS+4] = uy_rw[irw];
	  posdata[irw*NCOLS+5] = uz_rw[irw];
	  posdata[irw*NCOLS+6] = rho_rw[irw];
	  posdata[irw*NCOLS+7] = p_rw[irw];
	  posdata[irw*NCOLS+8] = c_rw[irw];
	}
	DataSet posf_dset = posf_h5->createDataSet(to_string(t),
						   PredType::NATIVE_DOUBLE,
						   posf_dspace);
	posf_dset.write(posdata, PredType::NATIVE_DOUBLE);
      }
    }
    for (int irw=0; irw < Nrw; ++irw){
      dx_rw = ux_rw[irw]*dt;
      dy_rw = uy_rw[irw]*dt;
      dz_rw = uz_rw[irw]*dt;

      // Second-order terms
      if (prm.int_order >= 2){
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

	if ((it+1) % int_dump_intv == 0){
	  rho_rw[irw] = intp.get_rho();
	  p_rw[irw] = intp.get_p();
	}
	
	// Second-order terms
	if (prm.int_order >= 2){
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
  
  statfile.close();
  declinedfile.close();

  return 0;

}
