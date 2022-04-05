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
#include "H5Cpp.h"
//#include "hdf5.h"
#include <ctime>

#include "io.hpp"
#include "utils.hpp"
#include "Parameters.hpp"

#include "ParticleSet.hpp"
#include "Topology.hpp"
#include "Integrator.hpp"
#include "ExplicitIntegrator.hpp"
#include "RKIntegrator.hpp"
#include "Initializer.hpp"
#include "helpers.hpp"
#include "MPIwrap.hpp"

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
  Parameters prm(argc, argv);

  std::string infilename = std::string(argv[1]);

  std::shared_ptr<Interpol> intp;
  set_interpolate_mode(intp, prm.mode, infilename);
  intp->set_U0(prm.U0);
  intp->set_int_order(prm.int_order);

  std::string folder = intp->get_folder();
  std::string rwfolder = folder + "/Interpolation/"; 
  if (mpi.rank() == 0)
    create_folder(rwfolder);
  std::string newfolder;
  if (prm.restart_folder != ""){
    newfolder = prm.folder;
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
  // prm.print();

  std::mt19937 gen;
  if (prm.random) {
    std::random_device rd;
    gen.seed(rd());
  }
  else {
    std::seed_seq rd{prm.seed + mpi.rank()};
    gen.seed(rd);
  }

  double t0 = std::max(intp->get_t_min(), prm.t0);
  if (mpi.rank() == 0){
    std::cout << "Testing interpolation..." << std::endl;
    test_interpolation(prm.Nrw, intp, newfolder, t0, gen);
  }

  return 0;
}
