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

#include "filaments_schema.hpp"

int main(int argc, char* argv[])
{

    std::cout << "Initialized FILAMENTS." << std::endl;

  // Input parameters
  if (argc < 2) {
    std::cout << "Specify an input file." << std::endl;
    return 1;
  }
  partrac::Params prm = partrac::parse_or_exit(filaments_schema(), argc, argv);

  std::string infilename = prm.input_file();

  std::shared_ptr<Interpol> intp;
  set_interpolate_mode(intp, prm.get<std::string>("mode"), infilename);
  intp->set_U0(prm.get<double>("U"));
  intp->set_int_order(prm.get<int>("int_order"));

  double Dm = prm.get<double>("Dm");
  double dt = prm.get<double>("dt");


  bool frozen_fields = prm.get<bool>("frozen_fields");

  std::string folder = intp->get_folder();
  std::string rwfolder = folder + "/Filaments/"; 
    create_folder(rwfolder);
  std::string newfolder;
  if (prm.get<std::string>("restart_folder") != ""){
    newfolder = prm.get<std::string>("folder");
  }
  else {
    newfolder = get_newfoldername(rwfolder, prm);
      create_folder(newfolder);
  }
  newfolder = newfolder + "" + "0" + "/";
  std::string posfolder = newfolder + "Positions/";
  std::string checkpointsfolder = newfolder + "Checkpoints/";
  create_folder(newfolder);
  create_folder(posfolder);
  create_folder(checkpointsfolder);
  prm.set<std::string>("folder", newfolder);

  if (prm.get<bool>("verbose"))
    prm.print();

  // Parallel generators
  std::vector<std::mt19937> gens = make_generators(prm);

  // TODO: These should not be stored in particle tracker parameters.
  prm.set<double>("Lx", intp->get_Lx());
  prm.set<double>("Ly", intp->get_Ly());
  prm.set<double>("Lz", intp->get_Lz());

  double t0 = std::max(intp->get_t_min(), prm.get<double>("t0"));
  double T = std::min(intp->get_t_max(), prm.get<double>("T"));
  if (frozen_fields)
    T = prm.get<double>("T");
  prm.set<double>("t0", t0);
  prm.set<double>("T", T);

  if (prm.get<bool>("inject") && prm.get<bool>("filter")){
      std::cout << "Cannot inject and filter at the same time (yet)." << std::endl;
    exit(1);
  }

  // Higher-order time integration?
  if (prm.get<int>("int_order") > 2){
      std::cout << "No support for such high temporal integration order." << std::endl;
    exit(1);
  }
  //if (prm.interpolation_test > 0){
  //  std::cout << "Testing interpolation..." << std::endl;
  //  test_interpolation(prm.interpolation_test, intp, newfolder, t0, gens[0]);
  //}

  if (frozen_fields)
    intp->update(prm.get<double>("t_frozen"));
  else
    intp->update(t0);

  std::shared_ptr<Integrator> integrator;
  if (prm.get<std::string>("scheme") == "explicit")
    integrator = std::make_shared<ExplicitIntegrator>(Dm, prm.get<int>("int_order"), gens);
  else if (prm.get<std::string>("scheme") == "RK4")
    integrator = std::make_shared<RK4Integrator>();
  else {
    std::cout << "Unrecognized (ODE integration) scheme: " << prm.get<std::string>("scheme") << std::endl;
    exit(1);
  }

  ParticleSet ps(intp, prm.get<Uint>("Nrw_max"));
  Topology mesh(ps, prm);

  if (prm.get<bool>("inject")){
    std::vector<std::string> key = split_string(prm.get<std::string>("init_mode"), "_");
    if (key[0] == "uniform"){
      std::cout << "Injection activated!" << std::endl;
    }
    else {
      std::cout << "init_mode " << prm.get<std::string>("init_mode") << " incompatible with injection." << std::endl;
      exit(1);
    }
  }

  if (prm.get<std::string>("restart_folder") != ""){
    mesh.load_checkpoint(prm.get<std::string>("restart_folder") + "/Checkpoints", prm);
  }
  else {
    std::shared_ptr<Initializer> init_state;
    std::vector<std::string> key = split_string(prm.get<std::string>("init_mode"), "_");
    if (key.size() == 0){
      std::cout << "init_mode not specified." << std::endl;
      exit(1);
    }
    else if (key[0] == "pair" || key[0] == "pairs"){
      init_state = std::make_shared<RandomPairsInitializer>(key, intp, prm, gens[0]);
    }
    else if (key[0] == "points"){
      init_state = std::make_shared<RandomPointsInitializer>(key, intp, prm, gens[0]);
    }
    else {
      std::cout << "Unknown init_mode: " << prm.get<std::string>("init_mode") << std::endl;
      exit(1);
    }
    mesh.load_initial_state(init_state, prm);
  }

  mesh.compute_maps();

  int it = 0;
  double t = t0;
  if (prm.get<std::string>("restart_folder") != ""){
    t = prm.get<double>("t");
  }

  prm.dump(newfolder, t);

  std::string h5fname = newfolder + "/data_from_t" + std::to_string(t) + ".h5";
  // no dump file at all when dumping is off
  H5::H5File h5f;
  if (prm.get<double>("dump_intv") > 0.){
    { H5::H5File create(h5fname.c_str(), H5F_ACC_TRUNC); }
    h5f.openFile(h5fname.c_str(), H5F_ACC_RDWR);
  }

  const double chunk_intv = prm.get<double>("dump_intv")*prm.get<int>("dump_chunk_size");

  std::map<std::string, bool> output_fields;
  output_fields["u"] = !prm.get<bool>("minimal_output");
  output_fields["c"] = !prm.get<bool>("minimal_output");
  output_fields["p"] = !prm.get<bool>("minimal_output") && prm.get<bool>("output_all_props");
  output_fields["rho"] = !prm.get<bool>("minimal_output") && prm.get<bool>("output_all_props");        
  output_fields["H"] = !prm.get<bool>("minimal_output") && mesh.dim() > 0;
  output_fields["n"] = !prm.get<bool>("minimal_output") && mesh.dim() > 1;

  std::ofstream statfile;
  if (prm.get<double>("stat_intv") > 0.){
    statfile.open(newfolder + "/tdata_from_t" + std::to_string(t) + ".dat");
    write_stats_header(statfile, mesh.dim());
  }
  

  // Simulation start
  std::clock_t clock_0 = std::clock();
  while (t <= T + dt/2){
    if (!frozen_fields)
      intp->update(t);

    // Statistics
    if (at_interval(it, prm.get<double>("stat_intv"), dt)){
      std::cout << "Time = " << t << std::endl;
      mesh.write_statistics(statfile, t, prm.get<double>("ds_max"), *integrator);
    }

    // Checkpoint
    if (at_interval(it, prm.get<double>("checkpoint_intv"), dt)){
      mesh.write_checkpoint(checkpointsfolder, t, prm);
    }

    // Resizing
    if (at_interval(it, prm.get<double>("resize_intv"), dt)){
      bool resized = mesh.resize(prm.get<double>("ds_max"));
      if (prm.get<bool>("verbose") && resized)
        std::cout << "Resized edges." << std::endl;
    }

    // Dump detailed data
    if (at_interval(it, prm.get<double>("dump_intv"), dt)){
      ps.update_fields(t, output_fields);

      std::string groupname = std::to_string(t);

      // Clear file if it exists, otherwise create
      if (at_interval(it, chunk_intv, dt) && it > 0){
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

  return 0;
}
