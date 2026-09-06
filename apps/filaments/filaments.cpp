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
#include "Params.hpp"

#include "ParticleSet.hpp"
#include "Topology.hpp"
#include "Integrator.hpp"
#include "ExplicitIntegrator.hpp"
#include "RKIntegrator.hpp"
#include "Initializer.hpp"
#include "MPIwrap.hpp"
#include "helpers.hpp"

// Parameters accepted by this app
partrac::Schema filaments_schema(){
  partrac::Schema s("filaments");
  s.require<std::string>("mode", "interpolator type");
  s.require<std::string>("init_mode", "initial distribution");
  s.require<double>("Dm", "molecular diffusivity");
  s.require<double>("dt", "timestep");
  s.require<double>("T", "final time");
  s.require<Uint>("Nrw", "number of particles");
  s.require<Uint>("Nrw_max", "max number of particles");
  s.require<int>("int_order", "integration order");
  s.require<double>("ds_max", "max edge length");
  s.require<double>("ds_min", "min edge length");
  s.require<double>("ds_init", "initial edge length");
  // only RandomPointsInitializer weights the sampling
  s.require_if<std::string>("init_weight",
                            [](const partrac::Params& p){
                              return p.get<std::string>("init_mode").rfind("points", 0) == 0;
                            },
                            "init_mode is a points distribution",
                            "sampling weight");
  s.opt<double>("t0", 0.0, "start time");
  s.opt<double>("U", 1.0, "velocity scale");
  s.opt<double>("x0", 0.0, "initial position");
  s.opt<double>("y0", 0.0, "initial position");
  s.opt<double>("z0", 0.0, "initial position");
  s.opt<double>("dump_intv", 100.0, "dump interval");
  s.opt<double>("stat_intv", 100.0, "statistics interval");
  s.opt<double>("checkpoint_intv", 1000.0, "checkpoint interval");
  s.opt<double>("resize_intv", 0.0, "resize interval");
  s.opt<double>("t_frozen", 0.0, "time to freeze the fields at");
  s.opt<double>("curv_refine_factor", 0.0, "curvature refinement factor");
  s.opt<int>("seed", 0, "random seed");
  s.opt<int>("dump_chunk_size", 0, "particles per dump chunk");
  s.opt<int>("filter_target", 0, "filter target");
  s.opt<bool>("verbose", false, "print the parameters");
  s.opt<bool>("random", true, "draw the seed randomly");
  s.opt<bool>("resize", false, "resize the filament");
  s.opt<bool>("filter", false, "filter the filament");
  s.opt<bool>("inject", false, "inject new particles");
  s.opt<bool>("inject_edges", true, "inject edges too");
  s.opt<bool>("frozen_fields", false, "freeze the velocity field");
  s.opt<bool>("local_dt", false, "use a local timestep");
  s.opt<bool>("cut_if_stuck", true, "cut edges that get stuck");
  s.opt<bool>("clear_initial_edges", false, "drop the initial edges");
  s.opt<bool>("output_all_props", true, "dump all properties");
  s.opt<bool>("minimal_output", false, "dump less");
  s.opt<std::string>("scheme", "explicit", "ODE integration scheme");
  s.opt<std::string>("tag", "", "appended to the folder name");
  s.opt<std::string>("restart_folder", "", "folder to restart from");
  s.runtime<std::string>("folder", "", "output folder");
  s.runtime<double>("t", 0.0, "current time");
  s.runtime<double>("Lx", 0.0, "domain size, from the interpolator");
  s.runtime<double>("Ly", 0.0, "domain size, from the interpolator");
  s.runtime<double>("Lz", 0.0, "domain size, from the interpolator");
  s.choices("mode", {"analytic", "structured", "lbm", "felbm", "fenics",
                     "tet", "triangle", "trianglefreq", "xdmftriangle", "xdmftet"});
  s.choices("scheme", {"explicit", "RK4"});
  s.check([](const partrac::Params& p){ return p.get<int>("int_order") <= 2; },
          "int_order must be 1 or 2");
  // RK4Integrator takes no arguments, so Dm is dropped
  s.warn([](const partrac::Params& p){
           return p.get<std::string>("scheme") == "RK4" && p.get<double>("Dm") != 0.0;
         },
         "scheme=RK4 ignores Dm");
  // dump_intv and stat_intv become integer step counts, so they must not round
  // down to zero
  s.finalize([](partrac::Params& p){
    const double dt = p.get<double>("dt");
    p.set<double>("dump_intv", std::max(p.get<double>("dump_intv"), dt));
    p.set<double>("stat_intv", std::max(p.get<double>("stat_intv"), dt));
    p.set<Uint>("Nrw_max", std::max(p.get<Uint>("Nrw_max"), p.get<Uint>("Nrw")));
  });
  return s;
}

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
  partrac::Params prm = partrac::parse_or_exit(filaments_schema(), argc, argv);

  std::string infilename = std::string(argv[1]);

  std::shared_ptr<Interpol> intp;
  set_interpolate_mode(intp, prm.get<std::string>("mode"), infilename);
  intp->set_U0(prm.get<double>("U"));
  intp->set_int_order(prm.get<int>("int_order"));

  double Dm = prm.get<double>("Dm");
  double dt = prm.get<double>("dt");

  bool resize = prm.get<bool>("resize");

  bool frozen_fields = prm.get<bool>("frozen_fields");
  bool local_dt = prm.get<bool>("local_dt");

  std::string folder = intp->get_folder();
  std::string rwfolder = folder + "/Filaments/"; 
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
  newfolder = newfolder + "" + std::to_string(mpi.rank()) + "/";
  std::string posfolder = newfolder + "Positions/";
  std::string checkpointsfolder = newfolder + "Checkpoints/";
  create_folder(newfolder);
  create_folder(posfolder);
  create_folder(checkpointsfolder);
  prm.set<std::string>("folder", newfolder);

  if (mpi.rank() == 0 && prm.get<bool>("verbose"))
    prm.print();

  // Parallel generators
  std::vector<std::mt19937> gens;
  for (int i=0, N=omp_get_max_threads(); i<N; ++i) {
    std::mt19937 gen;
    if (prm.get<bool>("random")) {
        std::random_device rd;
        gen.seed(rd());
    }
    else {
        std::seed_seq rd{prm.get<int>("seed") + omp_get_thread_num() };
        gen.seed(rd);
    }
    gens.emplace_back(gen);
  }

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
    if (mpi.rank() == 0)
      std::cout << "Cannot inject and filter at the same time (yet)." << std::endl;
    exit(0);
  }

  // Higher-order time integration?
  if (prm.get<int>("int_order") > 2){
    if (mpi.rank() == 0)
      std::cout << "No support for such high temporal integration order." << std::endl;
    exit(0);
  }
  //if (prm.interpolation_test > 0 && mpi.rank() == 0){
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
    exit(0);
  }

  ParticleSet ps(intp, prm.get<Uint>("Nrw_max"), mpi);
  Topology mesh(ps, prm, mpi);

  if (prm.get<bool>("inject")){
    std::vector<std::string> key = split_string(prm.get<std::string>("init_mode"), "_");
    if (key[0] == "uniform"){
      std::cout << "Injection activated!" << std::endl;
    }
    else {
      std::cout << "init_mode " << prm.get<std::string>("init_mode") << " incompatible with injection." << std::endl;
      exit(0);
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
      exit(0);
    }
    else if (key[0] == "pair" || key[0] == "pairs"){
      init_state = std::make_shared<RandomPairsInitializer>(key, intp, prm, mpi, gens[0]);
    }
    else if (key[0] == "points"){
      init_state = std::make_shared<RandomPointsInitializer>(key, intp, prm, mpi, gens[0]);
    }
    else {
      std::cout << "Unknown init_mode: " << prm.get<std::string>("init_mode") << std::endl;
      exit(0);
    }
    mesh.load_initial_state(init_state);
  }

  mesh.compute_maps();

  int it = 0;
  double t = t0;
  if (prm.get<std::string>("restart_folder") != ""){
    t = prm.get<double>("t");
  }

  prm.dump(newfolder, t);

  std::string h5fname = newfolder + "/data_from_t" + std::to_string(t) + ".h5";
  H5File h5f(h5fname.c_str(), H5F_ACC_TRUNC);

  Uint int_stat_intv = int(prm.get<double>("stat_intv")/dt);
  Uint int_dump_intv = int(prm.get<double>("dump_intv")/dt);
  Uint int_checkpoint_intv = int(prm.get<double>("checkpoint_intv")/dt);
  Uint int_chunk_intv = int_dump_intv*prm.get<int>("dump_chunk_size");
  Uint int_resize_intv = int(prm.get<double>("resize_intv")/dt);

  std::map<std::string, bool> output_fields;
  output_fields["u"] = !prm.get<bool>("minimal_output");
  output_fields["c"] = !prm.get<bool>("minimal_output");
  output_fields["p"] = !prm.get<bool>("minimal_output") && prm.get<bool>("output_all_props");
  output_fields["rho"] = !prm.get<bool>("minimal_output") && prm.get<bool>("output_all_props");        
  output_fields["H"] = !prm.get<bool>("minimal_output") && mesh.dim() > 0;
  output_fields["n"] = !prm.get<bool>("minimal_output") && mesh.dim() > 1;

  std::ofstream statfile(newfolder + "/tdata_from_t" + std::to_string(t) + ".dat");
  write_stats_header(mpi, statfile, mesh.dim());
  
  std::ofstream declinedfile(newfolder + "/declinedpos_from_t" + std::to_string(t) + ".dat");

  // Simulation start
  std::clock_t clock_0 = std::clock();
  while (t <= T + dt/2){
    if (!frozen_fields)
      intp->update(t);

    // Statistics
    if (it % int_stat_intv == 0){
      std::cout << "Time = " << t << std::endl;
      mesh.write_statistics(statfile, t, prm.get<double>("ds_max"), *integrator);
    }

    // Checkpoint
    if (it % int_checkpoint_intv == 0){
      mesh.write_checkpoint(checkpointsfolder, t, prm);
    }

    // Resizing
    if (resize && it % int_resize_intv == 0){
      bool resized = mesh.resize(prm.get<double>("ds_max"));
      if (prm.get<bool>("verbose") && resized)
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
