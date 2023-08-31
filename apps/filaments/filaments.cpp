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
#include "MPIwrap.hpp"
#include "helpers.hpp"

int main(int argc, char* argv[])
{
  MPIwrap mpi(argc, argv);

  if (mpi.rank() == 0)
    std::cout << "Initialized FILAMENTS with " << mpi.size() << " processes." << std::endl;
  mpi.barrier();

  // Input parameters
  if (argc < 2 && mpi.rank() == 0) {
    std::cout << "Specify an input file." << std::endl;
    return 0;
  }
  Parameters prm(argc, argv);
  if (prm.restart_folder != ""){
    prm.parse_file(prm.restart_folder + "/Checkpoints/params.dat");
    prm.parse_cmd(argc, argv);
  }

  std::string infilename = std::string(argv[1]);

  std::shared_ptr<Interpol> intp;
  set_interpolate_mode(intp, prm.mode, infilename);
  intp->set_U0(prm.U0);
  intp->set_int_order(prm.int_order);

  double Dm = prm.Dm;
  double dt = prm.dt;

  bool resize = prm.resize;

  bool frozen_fields = prm.frozen_fields;
  bool local_dt = prm.local_dt;

  std::string folder = intp->get_folder();
  std::string rwfolder = folder + "/Filaments/"; 
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
  newfolder = newfolder + "" + std::to_string(mpi.rank()) + "/";
  std::string posfolder = newfolder + "Positions/";
  std::string checkpointsfolder = newfolder + "Checkpoints/";
  create_folder(newfolder);
  create_folder(posfolder);
  create_folder(checkpointsfolder);
  prm.folder = newfolder;

  if (mpi.rank() == 0)
    prm.print();

  std::mt19937 gen;
  if (prm.random) {
    std::random_device rd;
    gen.seed(rd());
  }
  else {
    std::seed_seq rd{prm.seed + mpi.rank()};
    gen.seed(rd);
  }

  // TODO: These should not be stored in particle tracker parameters.
  prm.Lx = intp->get_Lx();
  prm.Ly = intp->get_Ly();
  prm.Lz = intp->get_Lz();

  double t0 = std::max(intp->get_t_min(), prm.t0);
  double T = std::min(intp->get_t_max(), prm.T);
  if (frozen_fields)
    T = prm.T;
  prm.t0 = t0;
  prm.T = T;

  if (prm.inject && prm.filter){
    if (mpi.rank() == 0)
      std::cout << "Cannot inject and filter at the same time (yet)." << std::endl;
    exit(0);
  }

  // Higher-order time integration?
  if (prm.int_order > 2){
    if (mpi.rank() == 0)
      std::cout << "No support for such high temporal integration order." << std::endl;
    exit(0);
  }
  if (prm.interpolation_test > 0 && mpi.rank() == 0){
    std::cout << "Testing interpolation..." << std::endl;
    test_interpolation(prm.interpolation_test, intp, newfolder, t0, gen);
  }

  if (frozen_fields)
    intp->update(prm.t_frozen);
  else
    intp->update(t0);

  std::shared_ptr<Integrator> integrator;
  if (prm.scheme == "explicit")
    integrator = std::make_shared<ExplicitIntegrator>(Dm, prm.int_order, gen);
  else if (prm.scheme == "RK4")
    integrator = std::make_shared<RK4Integrator>();
  else {
    std::cout << "Unrecognized (ODE integration) scheme: " << prm.scheme << std::endl;
    exit(0);
  }

  ParticleSet ps(intp, prm.Nrw_max, mpi);
  Topology mesh(ps, prm, mpi);

  if (prm.inject){
    std::vector<std::string> key = split_string(prm.init_mode, "_");
    if (key[0] == "uniform"){
      std::cout << "Injection activated!" << std::endl;
    }
    else {
      std::cout << "init_mode " << prm.init_mode << " incompatible with injection." << std::endl;
      exit(0);
    }
  }

  if (prm.restart_folder != ""){
    mesh.load_checkpoint(prm.restart_folder + "/Checkpoints", prm);
  }
  else {
    std::shared_ptr<Initializer> init_state;
    std::vector<std::string> key = split_string(prm.init_mode, "_");
    if (key.size() == 0){
      std::cout << "init_mode not specified." << std::endl;
      exit(0);
    }
    else if (key[0] == "pair" || key[0] == "pairs"){
      init_state = std::make_shared<RandomPairsInitializer>(key, intp, prm, mpi, gen);
    }
    else if (key[0] == "points"){
      init_state = std::make_shared<RandomPointsInitializer>(key, intp, prm, mpi, gen);
    }
    else {
      std::cout << "Unknown init_mode: " << prm.init_mode << std::endl;
      exit(0);
    }
    mesh.load_initial_state(init_state);
  }

  mesh.compute_maps();

  int it = 0;
  double t = t0;
  if (prm.restart_folder != ""){
    t = prm.t;
  }

  prm.dump(newfolder, t);

  std::string h5fname = newfolder + "/data_from_t" + std::to_string(t) + ".h5";
  H5File h5f(h5fname.c_str(), H5F_ACC_TRUNC);

  Uint int_stat_intv = int(prm.stat_intv/dt);
  Uint int_dump_intv = int(prm.dump_intv/dt);
  Uint int_checkpoint_intv = int(prm.checkpoint_intv/dt);
  Uint int_chunk_intv = int_dump_intv*prm.dump_chunk_size;
  Uint int_resize_intv = int(prm.resize_intv/dt);

  std::map<std::string, bool> output_fields;
  output_fields["u"] = !prm.minimal_output;
  output_fields["c"] = !prm.minimal_output;
  output_fields["p"] = !prm.minimal_output && prm.output_all_props;
  output_fields["rho"] = !prm.minimal_output && prm.output_all_props;        
  output_fields["H"] = !prm.minimal_output && mesh.dim() > 0;
  output_fields["n"] = !prm.minimal_output && mesh.dim() > 1;

  std::ofstream statfile(newfolder + "/tdata_from_t" + std::to_string(t) + ".dat");
  write_stats_header(mpi, statfile, mesh.dim());
  
  std::ofstream declinedfile(newfolder + "/declinedpos_from_t" + std::to_string(t) + ".dat");

  // Simulation start
  std::clock_t clock_0 = std::clock();
  while (t <= T){
    if (!frozen_fields)
      intp->update(t);

    // Statistics
    if (it % int_stat_intv == 0){
      std::cout << "Time = " << t << std::endl;
      mesh.write_statistics(statfile, t, prm.ds_max, *integrator);
    }

    // Checkpoint
    if (it % int_checkpoint_intv == 0){
      mesh.write_checkpoint(checkpointsfolder, t, prm);
    }

    // Resizing
    if (resize && it % int_resize_intv == 0){
      bool resized = mesh.resize(prm.ds_max);
      if (prm.verbose && resized)
        std::cout << "Resized edges." << std::endl;
    }

    // Dump detailed data
    if (it % int_dump_intv == 0){
      ps.update_fields(t, output_fields);

      std::string groupname = std::to_string(t);

      // Clear file if it exists, otherwise create
      if (int_chunk_intv > 0 && it % int_chunk_intv == 0 && it > 0){
        h5fname = newfolder + "/data_from_t" + std::to_string(t) + ".h5";
        h5f.openFile(h5fname.c_str(), H5F_ACC_TRUNC);
      }
      else {
        h5f.openFile(h5fname.c_str(), H5F_ACC_RDWR);
      }
      h5f.createGroup(groupname + "/");
      mesh.dump_hdf5(h5f, groupname, output_fields);
      h5f.close();
    }

    auto outside_nodes = integrator->step(ps, t, dt);

    if (outside_nodes.size() > 0)
      std::cout << "Some nodes are outside." << std::endl;

    t += dt;
    it += 1;
  }
  std::clock_t clock_1 = std::clock();
  double duration = (clock_1-clock_0) / (double) CLOCKS_PER_SEC;
  std::cout << "Total simulation time: " << duration << " seconds" << std::endl;

  mesh.write_checkpoint(checkpointsfolder, t, prm);

  // Close files
  statfile.close();
  declinedfile.close();

  return 0;
}
