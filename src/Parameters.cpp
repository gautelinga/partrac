#include "Parameters.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

Parameters::Parameters(int argc, char* argv[])
{
  parse_cmd(argc, argv);
}

void Parameters::set_traj_defaults()
{
  // Default parameters
  set("mode", "structured");
  set("folder", "");
  set("restart_folder", "");
  set("Dm", 1.0);
  set("t0", 0.0);
  set("t", 0.0);
  set("T", 10000000.0);
  set("Nrw", 100);
  set("dump_intv", 100.0);
  set("stat_intv", 100.0);
  set("checkpoint_intv", 1000.0);
  set("dump_chunk_size", 50);
  set("verbose", false);
  set("U0", 1.0);
  //
  set("x0", 0.0);
  set("y0", 0.0);
  set("z0", 0.0);
  //
  set("La", 0.0);
  set("Lb", 0.0);
  //
  set("Lx", 0.0);
  set("Ly", 0.0);
  set("Lz", 0.0);
  set("nx", 0);
  set("ny", 0);
  set("nz", 0);
  //
  set("filter", false);
  set("filter_intv", 0.0);
  set("filter_target", 0);
  //
  set("dt", 1.0);
  set("int_order", 1);
  set("init_mode", "line_x");  // what else?
  set("init_weight", "none");
  set("write_mode", "hdf5");  // or text
  set("interpolation_test", 0);
  set("n_accepted", 0);
  set("n_declined", 0);
  set("refine", false);
  set("coarsen", false);
  set("refine_intv", 100.0);
  set("coarsen_intv", 1000.0);
  set("hist_chunk_size", 10);
  set("ds_max", 1.0);
  set("ds_min", 0.1);
  set("ds_init", 1.0);
  set("Nrw_max", -1);
  set("curv_refine_factor", 0.0);
  set("output_all_props", true);
  set("minimal_output", false);
  set("inject", false);
  set("inject_intv", 0.0);
  set("T_inject", 1e10);
  set("inject_edges", true);
  //
  set("clear_initial_edges", false);
  //
  set("seed", 0);
  //
  set("tag", "");
  //
  set("random", true);
}

void Parameters::parse_cmd(int argc, char* argv[])
{
  std::size_t found;
  std::string argstr, key, val;
  for (int iarg=2; iarg < argc; ++iarg){
    argstr = argv[iarg];
    found = argstr.find('=');
    if (found != std::string::npos){
      key = argstr.substr(0, found);
      val = argstr.substr(found+1);
      //boost::trim(key);
      //boost::trim(val);
      trim(key);
      trim(val);
      std::cout << "set " << key << " to " << val << '\n';
      set(key, val);
    }
  }
  set("dump_intv", std::max(getd("dump_intv"), getd("dt")));
  set("stat_intv", std::max(getd("stat_intv"), getd("dt")));
  set("Nrw_max", std::max(getui("Nrw_max"), getui("Nrw")));
}

void Parameters::parse_file(std::string infile)
{
  std::ifstream input(infile);
  if (!input.good()){
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
      // boost::trim(key);
      // boost::trim(val);
      trim(key);
      trim(val);
      set(key, val);
    }
  }
}

void Parameters::print(std::ostream& out) const
{
  // FIXME only print if verbose?
  for (const auto& [k, v] : dict_){
    out << std::left << std::setw(18) << k<< " = ";
    std::visit([&out] (auto&& x)
    {
      using T = std::decay_t<decltype(x)>;
      if constexpr (std::is_same_v<bool, T>)
                     out << (x ? "true" : "false");
      else
        out << x;
    }, v);
    out << '\n';
  }
}

void Parameters::dump(std::string dumpfolder) const
{
  std::string filename = dumpfolder + "/params.dat";
  write_params_to_file(filename);
}

void Parameters::dump(std::string dumpfolder, const double t) const
{
  std::string filename = dumpfolder + "/params_from_t" + std::to_string(t) + ".dat";
  write_params_to_file(filename);
}

bool Parameters::contains(std::string k) const
{
  return dict_.count(k) == 1;
}

void Parameters::write_params_to_file(std::string filename) const
{
  std::ofstream paramsfile(filename);
  assert(paramsfile.good());
  print(paramsfile);
  paramsfile.close();
}

void Parameters::ltrim(std::string& s)
{
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
    return !std::isspace(ch);
  }));

}

void Parameters::rtrim(std::string& s)
{
  s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
    return !std::isspace(ch);
  }).base(), s.end());
}

void Parameters::trim(std::string& s)
{
  rtrim(s);
  ltrim(s);
}
