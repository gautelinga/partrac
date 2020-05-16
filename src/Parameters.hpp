#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include "io.hpp"

#ifndef __PARAMETERS_HPP
#define __PARAMETERS_HPP

using namespace std;

class Parameters {
public:
  Parameters(int argc, char* argv[]);
  void parse_cmd(int argc, char* argv[]);
  void parse_file(string infile);
  void set_param(string key, string val);
  void print();
  void dump(string folder);
  void dump(string folder, const double t);
  // Default parameters
  string folder = "";
  string restart_folder = "";
  double Dm = 1.0;
  double t0 = 0.;
  double t = 0.;
  double T = 10000000.;
  int Nrw = 100;
  double dump_intv = 100.0;
  double stat_intv = 100.0;
  double checkpoint_intv = 1000.0;
  int dump_chunk_size = 50;
  bool verbose = false;
  double U0 = 1.0;
  //
  double x0 = 0.0;
  double y0 = 0.0;
  double z0 = 0.0;
  //
  double La = 0.0;
  double Lb = 0.0;
  //
  double Lx = 0.0;
  double Ly = 0.0;
  double Lz = 0.0;
  int nx = 0;
  int ny = 0;
  int nz = 0;
  double dt = 1.0;
  int int_order = 1;
  string init_mode = "line_x";  // what else?
  string init_weight = "none";
  string write_mode = "hdf5";  // or text
  int interpolation_test = 0;
  long int n_accepted = 0;
  long int n_declined = 0;
  bool refine = false;
  double refine_intv = 100.0;
  int hist_chunk_size = 10;
  double ds_max = 1.0;
  int Nrw_max = -1;
private:
  void write_params_to_file(string);
};

Parameters::Parameters(int argc, char* argv[]){
  parse_cmd(argc, argv);
}

void Parameters::parse_cmd(int argc, char* argv[]){
  size_t found;
  string argstr, key, val;
  for (int iarg=2; iarg < argc; ++iarg){
    argstr = argv[iarg];
    found = argstr.find('=');
    if (found != string::npos){
      key = argstr.substr(0, found);
      val = argstr.substr(found+1);
      boost::trim(key);
      boost::trim(val);
      set_param(key, val);
    }
  }
  dump_intv = max(dump_intv, dt);
  stat_intv = max(stat_intv, dt);
  Nrw_max = max(Nrw_max, Nrw);
}

void Parameters::parse_file(string infile){
  ifstream input(infile);
  if (!input){
    cout << "File " << infile <<" doesn't exist." << endl;
    exit(0);
  }
  size_t found;
  string key, val;
  for (string line; getline(input, line); ){
    found = line.find('=');
    if (found != string::npos){
      key = line.substr(0, found);
      val = line.substr(found+1);
      boost::trim(key);
      boost::trim(val);
      set_param(key, val);
    }
  }
}

void Parameters::set_param(string key, string val){
  if (key == "Dm") Dm = stod(val);
  if (key == "t0") t0 = stod(val);
  if (key == "T") T = stod(val);
  if (key == "dt") dt = stod(val);
  if (key == "Nrw") Nrw = stoi(val);
  if (key == "dump_intv") dump_intv = stod(val);
  if (key == "stat_intv") stat_intv = stod(val);
  if (key == "checkpoint_intv") checkpoint_intv = stod(val);
  if (key == "verbose") verbose = (val == "true" || val == "True") ? true : false;
  if (key == "U") U0 = stod(val);
  if (key == "x0") x0 = stod(val);
  if (key == "y0") y0 = stod(val);
  if (key == "z0") z0 = stod(val);
  if (key == "La") La = stod(val);
  if (key == "Lb") Lb = stod(val);
  if (key == "int_order") int_order = stoi(val);
  if (key == "init_mode") init_mode = val;
  if (key == "init_weight") init_weight = val;
  if (key == "dump_mode") write_mode = val;
  if (key == "interpolation_test") interpolation_test = stoi(val);
  if (key == "dump_chunk_size") dump_chunk_size = stoi(val);
  if (key == "n_accepted") n_accepted = stoi(val);
  if (key == "n_declined") n_declined = stoi(val);
  if (key == "restart_folder") restart_folder = val;
  if (key == "folder") folder = val;
  if (key == "t") t = stod(val);

  if (key == "refine") refine = (val == "true" || val == "True") ? true : false;
  if (key == "refine_intv") refine_intv = stod(val);
  if (key == "hist_chunk_size") hist_chunk_size = stoi(val);
  if (key == "ds_max") ds_max = stod(val);
  if (key == "Nrw_max") Nrw_max = stoi(val);
}

void Parameters::print(){
  if (verbose){
    print_param("Dm             ", Dm);
    print_param("t0             ", t0);
    print_param("T              ", T);
    print_param("dt             ", dt);
    print_param("Nrw            ", Nrw);
    print_param("dump_intv      ", dump_intv);
    print_param("stat_intv      ", stat_intv);
    print_param("checkpoint_intv", checkpoint_intv);
    print_param("verbose        ", verbose ? "true" : "false");
    print_param("U              ", U0);
    print_param("x0             ", x0);
    print_param("y0             ", y0);
    print_param("z0             ", z0);
    print_param("La             ", La);
    print_param("Lb             ", Lb);
    print_param("int_order      ", int_order);
    print_param("init_mode      ", init_mode);
    print_param("init_weight    ", init_weight);
    print_param("dump_mode      ", write_mode);
    print_param("restart_folder ", restart_folder);
    print_param("folder         ", folder);

    print_param("refine         ", refine ? "true" : "false");
    print_param("refine_intv    ", refine_intv);
    print_param("hist_chunk_size", hist_chunk_size);
    print_param("ds_max         ", ds_max);
    print_param("Nrw_max        ", Nrw_max);
  }
}

void Parameters::dump(string dumpfolder){
  string filename = dumpfolder + "/params.dat";
  write_params_to_file(filename);
}

void Parameters::dump(string dumpfolder, const double t){
  string filename = dumpfolder + "/params_from_t" + to_string(t) + ".dat";
  write_params_to_file(filename);
}

void Parameters::write_params_to_file(string filename){
  ofstream paramsfile(filename);
  write_param(paramsfile, "Dm", Dm);
  write_param(paramsfile, "t0", t0);
  write_param(paramsfile, "t", t);
  write_param(paramsfile, "T", T);
  write_param(paramsfile, "dt", dt);
  write_param(paramsfile, "Nrw", Nrw);
  write_param(paramsfile, "dump_intv", dump_intv);
  write_param(paramsfile, "stat_intv", stat_intv);
  write_param(paramsfile, "checkpoint_intv", checkpoint_intv);
  write_param(paramsfile, "verbose", verbose ? "true" : "false");
  write_param(paramsfile, "U", U0);
  write_param(paramsfile, "x0", x0);
  write_param(paramsfile, "y0", y0);
  write_param(paramsfile, "z0", z0);
  write_param(paramsfile, "La", La);
  write_param(paramsfile, "Lb", Lb);
  write_param(paramsfile, "Lx", Lx);
  write_param(paramsfile, "Ly", Ly);
  write_param(paramsfile, "Lz", Lz);
  write_param(paramsfile, "nx", nx);
  write_param(paramsfile, "ny", ny);
  write_param(paramsfile, "nz", nz);
  write_param(paramsfile, "int_order", int_order);
  write_param(paramsfile, "init_mode", init_mode);
  write_param(paramsfile, "init_weight", init_weight);
  write_param(paramsfile, "write_mode", write_mode);
  write_param(paramsfile, "interpolation_test", interpolation_test);
  write_param(paramsfile, "dump_chunk_size", dump_chunk_size);
  write_param(paramsfile, "n_accepted", n_accepted);
  write_param(paramsfile, "n_declined", n_declined);
  write_param(paramsfile, "restart_folder", restart_folder);
  write_param(paramsfile, "folder", folder);

  write_param(paramsfile, "refine", refine ? "true" : "false");
  write_param(paramsfile, "refine_intv", refine_intv);
  write_param(paramsfile, "hist_chunk_size", hist_chunk_size);
  write_param(paramsfile, "ds_max", ds_max);
  write_param(paramsfile, "Nrw_max", Nrw_max);
  paramsfile.close();
}

#endif
