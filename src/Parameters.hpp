#ifndef __PARAMETERS_HPP
#define __PARAMETERS_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include "io.hpp"

double stodouble(const std::string val){
  std::stringstream ss(val);
  double d;
  ss >> d;
  if (ss.fail()){
    std::cout << "Unable to format " + val + " as a double." << std::endl;
    exit(0);
  }
  return d;
}

int stoint(const std::string val){
  double d = stodouble(val);
  int i = int(d);
  if (d - i != 0){
    std::cout << "The value " + val + " is not an int!" << std::endl;
    exit(0);
  }
  return i;
}

bool stobool(const std::string val){
  return (val == "true" || val == "True") ? true : false;
}

std::string bool2string(const bool val){
  return (val ? "true" : "false");
}

class Parameters {
public:
  Parameters(int argc, char* argv[]);
  void parse_cmd(int argc, char* argv[]);
  void parse_file(std::string infile);
  void set_param(std::string key, std::string val);
  void print();
  void dump(std::string folder);
  void dump(std::string folder, const double t);
  // Default parameters
  std::string mode = "structured";
  std::string folder = "";
  std::string restart_folder = "";
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
  //
  bool filter = false;
  double filter_intv = 0.0;
  int filter_target = 0;
  //
  bool frozen_fields = false;
  bool local_dt = false;
  double dl_max = 1.0;
  double u_eps = 1e-7;
  double t_frozen = 0.;
  bool cut_if_stuck = true;
  //
  double dt = 1.0;
  int int_order = 1;
  std::string init_mode = "line_x";  // what else?
  std::string init_weight = "none";
  std::string write_mode = "hdf5";  // or text
  int interpolation_test = 0;
  long int n_accepted = 0;
  long int n_declined = 0;
  bool refine = false;
  bool coarsen = false;
  double refine_intv = 100.0;
  double coarsen_intv = 1000.0;
  int hist_chunk_size = 10;
  double ds_max = 1.0;
  double ds_min = 0.1;
  double ds_init = 1.0;
  int Nrw_max = -1;
  double curv_refine_factor = 0.0;
  bool output_all_props = true;
  bool minimal_output = false;
  bool inject = false;
  double inject_intv = 0.0;
  double T_inject = 1e10;
  bool inject_edges = true;
  //
  bool clear_initial_edges = false;
  //
  int seed = 0;
  //
  std::string tag = "";
  //
  bool random = true;
private:
  void write_params_to_file(std::string);
};

Parameters::Parameters(int argc, char* argv[]){
  parse_cmd(argc, argv);
}

void Parameters::parse_cmd(int argc, char* argv[]){
  size_t found;
  std::string argstr, key, val;
  for (int iarg=2; iarg < argc; ++iarg){
    argstr = argv[iarg];
    found = argstr.find('=');
    if (found != std::string::npos){
      key = argstr.substr(0, found);
      val = argstr.substr(found+1);
      boost::trim(key);
      boost::trim(val);
      set_param(key, val);
    }
  }
  dump_intv = std::max(dump_intv, dt);
  stat_intv = std::max(stat_intv, dt);
  Nrw_max = std::max(Nrw_max, Nrw);
}

void Parameters::parse_file(std::string infile){
  std::ifstream input(infile);
  if (!input){
    std::cout << "File " << infile <<" doesn't exist." << std::endl;
    exit(0);
  }
  size_t found;
  std::string key, val;
  for (std::string line; getline(input, line); ){
    found = line.find('=');
    if (found != std::string::npos){
      key = line.substr(0, found);
      val = line.substr(found+1);
      boost::trim(key);
      boost::trim(val);
      set_param(key, val);
    }
  }
}

void Parameters::set_param(std::string key, std::string val){
  if (key == "Dm") Dm = stodouble(val);
  if (key == "t0") t0 = stodouble(val);
  if (key == "T") T = stodouble(val);
  if (key == "dt") dt = stodouble(val);
  if (key == "Nrw") Nrw = stoint(val);
  if (key == "dump_intv") dump_intv = stodouble(val);
  if (key == "stat_intv") stat_intv = stodouble(val);
  if (key == "checkpoint_intv") checkpoint_intv = stodouble(val);
  if (key == "verbose") verbose = stobool(val);
  if (key == "U") U0 = stodouble(val);
  if (key == "x0") x0 = stodouble(val);
  if (key == "y0") y0 = stodouble(val);
  if (key == "z0") z0 = stodouble(val);
  if (key == "La") La = stodouble(val);
  if (key == "Lb") Lb = stodouble(val);
  if (key == "int_order") int_order = stoint(val);
  if (key == "init_mode") init_mode = val;
  if (key == "init_weight") init_weight = val;
  if (key == "dump_mode") write_mode = val;
  if (key == "interpolation_test") interpolation_test = stoint(val);
  if (key == "dump_chunk_size") dump_chunk_size = stoint(val);
  if (key == "n_accepted") n_accepted = stoint(val);  // may be too large?
  if (key == "n_declined") n_declined = stoint(val);
  if (key == "restart_folder") restart_folder = val;
  if (key == "mode") mode = val;
  if (key == "folder") folder = val;
  if (key == "t") t = stodouble(val);

  if (key == "refine") refine = stobool(val);
  if (key == "refine_intv") refine_intv = stodouble(val);
  if (key == "coarsen") coarsen = stobool(val);
  if (key == "coarsen_intv") coarsen_intv = stodouble(val);
  if (key == "hist_chunk_size") hist_chunk_size = stoint(val);
  if (key == "ds_max") ds_max = stodouble(val);
  if (key == "ds_min") ds_min = stodouble(val);
  if (key == "ds_init") ds_init = stodouble(val);
  if (key == "Nrw_max") Nrw_max = stoint(val);
  if (key == "curv_refine_factor") curv_refine_factor = stodouble(val);
  if (key == "output_all_props") output_all_props = stobool(val);
  if (key == "minimal_output") minimal_output = stobool(val);

  if (key == "filter") filter = stobool(val);
  if (key == "filter_intv") filter_intv = stodouble(val);
  if (key == "filter_target") filter_target = stoint(val);

  if (key == "frozen_fields") frozen_fields = stobool(val);
  if (key == "local_dt") local_dt = stobool(val);
  if (key == "dl_max") dl_max = stodouble(val);
  if (key == "u_eps") u_eps = stodouble(val);
  if (key == "t_frozen") t_frozen = stodouble(val);
  if (key == "cut_if_stuck") cut_if_stuck = stobool(val);

  if (key == "inject") inject = stobool(val);
  if (key == "inject_intv") inject_intv = stodouble(val);
  if (key == "T_inject") T_inject = stodouble(val);
  if (key == "inject_edges") inject_edges = stobool(val);

  if (key == "clear_initial_edges") clear_initial_edges = stobool(val);
  if (key == "seed") seed = stoint(val);

  if (key == "tag") tag = val;

  if (key == "random") random = stobool(val);
}

void Parameters::print(){
  if (verbose){
    print_param("Dm                 ", Dm);
    print_param("t0                 ", t0);
    print_param("T                  ", T);
    print_param("dt                 ", dt);
    print_param("Nrw                ", Nrw);
    print_param("dump_intv          ", dump_intv);
    print_param("stat_intv          ", stat_intv);
    print_param("checkpoint_intv    ", checkpoint_intv);
    print_param("verbose            ", bool2string(verbose));
    print_param("U                  ", U0);
    print_param("x0                 ", x0);
    print_param("y0                 ", y0);
    print_param("z0                 ", z0);
    print_param("La                 ", La);
    print_param("Lb                 ", Lb);
    print_param("int_order          ", int_order);
    print_param("init_mode          ", init_mode);
    print_param("init_weight        ", init_weight);
    print_param("dump_mode          ", write_mode);
    print_param("restart_folder     ", restart_folder);
    print_param("mode               ", mode);
    print_param("folder             ", folder);

    print_param("refine             ", bool2string(refine));
    print_param("refine_intv        ", refine_intv);
    print_param("coarsen            ", bool2string(coarsen));
    print_param("coarsen_intv       ", coarsen_intv);
    print_param("hist_chunk_size    ", hist_chunk_size);
    print_param("ds_max             ", ds_max);
    print_param("ds_min             ", ds_min);
    print_param("ds_init            ", ds_init);
    print_param("Nrw_max            ", Nrw_max);
    print_param("curv_refine_factor ", curv_refine_factor);
    print_param("output_all_props   ", bool2string(output_all_props));
    print_param("minimal_output     ", bool2string(minimal_output));

    print_param("filter             ", bool2string(filter));
    print_param("filter_intv        ", filter_intv);
    print_param("filter_target      ", filter_target);

    print_param("frozen_fields      ", bool2string(frozen_fields));
    print_param("local_dt           ", bool2string(local_dt));
    print_param("dl_max             ", dl_max);
    print_param("u_eps              ", u_eps);
    print_param("t_frozen           ", t_frozen);
    print_param("cut_if_stuck", bool2string(cut_if_stuck));

    print_param("inject             ", bool2string(inject));
    print_param("inject_intv        ", inject_intv);
    print_param("T_inject           ", T_inject);
    print_param("inject_edges       ", bool2string(inject_edges));

    print_param("clear_initial_edges", bool2string(clear_initial_edges));
    print_param("seed               ", seed);

    print_param("tag                ", tag);
    print_param("random             ", bool2string(random));
  }
}

void Parameters::dump(std::string dumpfolder){
  std::string filename = dumpfolder + "/params.dat";
  write_params_to_file(filename);
}

void Parameters::dump(std::string dumpfolder, const double t){
  std::string filename = dumpfolder + "/params_from_t" + std::to_string(t) + ".dat";
  write_params_to_file(filename);
}

void Parameters::write_params_to_file(std::string filename){
  std::ofstream paramsfile(filename);
  write_param(paramsfile, "Dm", Dm);
  write_param(paramsfile, "t0", t0);
  write_param(paramsfile, "t", t);
  write_param(paramsfile, "T", T);
  write_param(paramsfile, "dt", dt);
  write_param(paramsfile, "Nrw", Nrw);
  write_param(paramsfile, "dump_intv", dump_intv);
  write_param(paramsfile, "stat_intv", stat_intv);
  write_param(paramsfile, "checkpoint_intv", checkpoint_intv);
  write_param(paramsfile, "verbose", bool2string(verbose));
  write_param(paramsfile, "U", U0);
  write_param(paramsfile, "x0", x0);
  write_param(paramsfile, "y0", y0);
  write_param(paramsfile, "z0", z0);
  write_param(paramsfile, "La", La);
  write_param(paramsfile, "Lb", Lb);
  write_param(paramsfile, "Lx", Lx);
  write_param(paramsfile, "Ly", Ly);
  write_param(paramsfile, "Lz", Lz);
  write_param(paramsfile, "int_order", int_order);
  write_param(paramsfile, "init_mode", init_mode);
  write_param(paramsfile, "init_weight", init_weight);
  write_param(paramsfile, "write_mode", write_mode);
  write_param(paramsfile, "interpolation_test", interpolation_test);
  write_param(paramsfile, "dump_chunk_size", dump_chunk_size);
  write_param(paramsfile, "n_accepted", n_accepted);
  write_param(paramsfile, "n_declined", n_declined);
  write_param(paramsfile, "restart_folder", restart_folder);
  write_param(paramsfile, "mode", mode);
  write_param(paramsfile, "folder", folder);

  write_param(paramsfile, "refine", bool2string(refine));
  write_param(paramsfile, "refine_intv", refine_intv);
  write_param(paramsfile, "coarsen", bool2string(coarsen));
  write_param(paramsfile, "coarsen_intv", coarsen_intv);
  write_param(paramsfile, "hist_chunk_size", hist_chunk_size);
  write_param(paramsfile, "ds_max", ds_max);
  write_param(paramsfile, "ds_min", ds_min);
  write_param(paramsfile, "ds_init", ds_init);
  write_param(paramsfile, "Nrw_max", Nrw_max);
  write_param(paramsfile, "curv_refine_factor", curv_refine_factor);
  write_param(paramsfile, "output_all_props", bool2string(output_all_props));
  write_param(paramsfile, "minimal_output", bool2string(minimal_output));

  write_param(paramsfile, "filter", bool2string(filter));
  write_param(paramsfile, "filter_intv", filter_intv);
  write_param(paramsfile, "filter_target", filter_target);

  write_param(paramsfile, "frozen_fields", bool2string(frozen_fields));
  write_param(paramsfile, "local_dt", bool2string(local_dt));
  write_param(paramsfile, "dl_max", dl_max);
  write_param(paramsfile, "u_eps", u_eps);
  write_param(paramsfile, "t_frozen", t_frozen);
  write_param(paramsfile, "cut_if_stuck", bool2string(cut_if_stuck));

  write_param(paramsfile, "inject", bool2string(inject));
  write_param(paramsfile, "inject_intv", inject_intv);
  write_param(paramsfile, "T_inject", T_inject);
  write_param(paramsfile, "inject_edges", bool2string(inject_edges));

  write_param(paramsfile, "clear_initial_edges", bool2string(clear_initial_edges));
  write_param(paramsfile, "seed", seed);

  write_param(paramsfile, "tag", tag);
  write_param(paramsfile, "random", bool2string(random));

  paramsfile.close();
}

#endif
