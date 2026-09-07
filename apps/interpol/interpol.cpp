#include <iostream>
#include <vector>
#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <sstream>
#include <random>
#include <cmath>
#include <set>
#include <iterator>
//#include "H5Cpp.h"
//#include "hdf5.h"
#include <ctime>

#include "io.hpp"
#include "rng.hpp"
#include "utils.hpp"
#include "Params.hpp"

#include "ParticleSet.hpp"
#include "Topology.hpp"
#include "Integrator.hpp"
#include "ExplicitIntegrator.hpp"
#include "RKIntegrator.hpp"
#include "Initializer.hpp"
#include "helpers.hpp"

#include "interpol_schema.hpp"

int main(int argc, char* argv[])
{

    std::cout << "Initialized Interpolator." << std::endl;

  // Input parameters
  if (argc < 2) {
    std::cout << "Specify an input file." << std::endl;
    return 1;
  }
  partrac::Params prm = partrac::parse_or_exit(interpol_schema(), argc, argv);

  if (prm.get<int>("num_threads") > 0){
      omp_set_dynamic(0);
      omp_set_num_threads(prm.get<int>("num_threads"));
  }

  std::string infilename = prm.input_file();

  std::shared_ptr<Interpol> intp;
  set_interpolate_mode(intp, prm.get<std::string>("mode"), infilename);
  intp->set_U0(prm.get<double>("U"));
  intp->set_int_order(prm.get<int>("int_order"));

  std::string folder = intp->get_folder();
  std::string rwfolder = folder + "/Interpolation/"; 
    create_folder(rwfolder);
  std::string newfolder;
  if (prm.get<std::string>("restart_folder") != ""){
    newfolder = prm.get<std::string>("folder");
  }
  else {
    newfolder = get_newfoldername(rwfolder, prm);
      create_folder(newfolder);
  }

  // Parallel generators
  std::vector<std::mt19937> gens = make_generators(prm);

  double t0 = std::max(intp->get_t_min(), prm.get<double>("t0"));
  std::cout << "Testing interpolation..." << std::endl;
  test_interpolation(prm.get<Uint>("Nrw"), intp, newfolder, t0, gens);

  return 0;
}
