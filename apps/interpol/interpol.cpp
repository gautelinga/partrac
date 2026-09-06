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
#include "utils.hpp"
#include "Params.hpp"

#include "ParticleSet.hpp"
#include "Topology.hpp"
#include "Integrator.hpp"
#include "ExplicitIntegrator.hpp"
#include "RKIntegrator.hpp"
#include "Initializer.hpp"
#include "helpers.hpp"
#include "MPIwrap.hpp"

// Parameters accepted by this app
partrac::Schema interpol_schema(){
  partrac::Schema s("interpol");
  s.require<std::string>("mode", "interpolator type");
  s.require<Uint>("Nrw", "number of probe points");
  s.require<int>("int_order", "interpolation order");
  s.opt<double>("Dm", 0.0, "diffusivity, enters the folder name only");
  s.opt<double>("dt", 1.0, "timestep, enters the folder name only");
  s.opt<double>("t0", 0.0, "time to probe the field at");
  s.opt<double>("U", 1.0, "velocity scale");
  s.opt<int>("seed", 0, "random seed");
  s.opt<int>("num_threads", 0, "OpenMP threads, 0 = leave alone");
  s.opt<bool>("random", true, "draw the seed randomly");
  s.opt<std::string>("tag", "", "appended to the folder name");
  s.opt<std::string>("restart_folder", "", "folder to restart from");
  s.runtime<std::string>("folder", "", "output folder");
  s.choices("mode", {"analytic", "structured", "lbm", "felbm", "fenics",
                     "tet", "triangle", "trianglefreq", "xdmftriangle", "xdmftet"});
  return s;
}

int main(int argc, char* argv[])
{
  MPIwrap mpi(argc, argv);

  if (mpi.rank() == 0)
    std::cout << "Initialized Interpolator with " << mpi.size() << " processes." << std::endl;

  // Input parameters
  if (argc < 2 && mpi.rank() == 0) {
    std::cout << "Specify an input file." << std::endl;
    return 0;
  }
  partrac::Params prm = partrac::parse_or_exit(interpol_schema(), argc, argv);

  if (prm.get<int>("num_threads") > 0){
      omp_set_dynamic(0);
      omp_set_num_threads(prm.get<int>("num_threads"));
  }

  std::string infilename = std::string(argv[1]);

  std::shared_ptr<Interpol> intp;
  set_interpolate_mode(intp, prm.get<std::string>("mode"), infilename);
  intp->set_U0(prm.get<double>("U"));
  intp->set_int_order(prm.get<int>("int_order"));

  std::string folder = intp->get_folder();
  std::string rwfolder = folder + "/Interpolation/"; 
  if (mpi.rank() == 0)
    create_folder(rwfolder);
  std::string newfolder;
  if (prm.get<std::string>("restart_folder") != ""){
    newfolder = prm.get<std::string>("folder");
  }
  else {
    newfolder = get_newfoldername(rwfolder, prm);
    mpi.barrier();
    if (mpi.rank() == 0)
      create_folder(newfolder);
    mpi.barrier();
  }
  if (mpi.rank() == 0){

  }

  // Parallel generators
  std::vector<std::mt19937> gens;
  for (int i=0, N=omp_get_max_threads(); i<N; ++i) {
    std::mt19937 gen;
    if (prm.get<bool>("random")) {
        std::random_device rd;
        gen.seed(rd());
    }
    else {
        std::seed_seq rd{prm.get<int>("seed") + i};
        gen.seed(rd);
    }
    gens.emplace_back(gen);
  }

  double t0 = std::max(intp->get_t_min(), prm.get<double>("t0"));
  if (mpi.rank() == 0){
    std::cout << "Testing interpolation..." << std::endl;
    test_interpolation(prm.get<Uint>("Nrw"), intp, newfolder, t0, gens);
  }

  return 0;
}
