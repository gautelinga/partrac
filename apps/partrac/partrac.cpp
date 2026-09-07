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
#include <omp.h>

#include "io.hpp"
#include "rng.hpp"
#include "utils.hpp"
#include "Params.hpp"

#include "ParticleSet.hpp"
#include "Topology.hpp"

#include "Integrator.hpp"
#include "ExplicitIntegrator.hpp"
#include "RKIntegrator.hpp"
#include "helpers.hpp"

#include "partrac_schema.hpp"

int main(int argc, char* argv[])
{

    std::cout << "Initialized Partrac." << std::endl;

  // Input parameters
  if (argc < 2) {
    std::cout << "Specify an input file." << std::endl;
    return 1;
  }
  partrac::Params prm = partrac::parse_or_exit(partrac_schema(), argc, argv);

  if (prm.get<int>("num_threads") > 0){
      omp_set_dynamic(0);
      omp_set_num_threads(prm.get<int>("num_threads"));
  }

  std::string infilename = prm.input_file();

  std::cout << "Setting interpolator..." << std::endl;

  std::shared_ptr<Interpol> intp;
  set_interpolate_mode(intp, prm.get<std::string>("mode"), infilename);
  
  intp->set_U0(prm.get<double>("U"));
  intp->set_int_order(prm.get<int>("int_order"));

  double Dm = prm.get<double>("Dm");
  double dt = prm.get<double>("dt");

  bool refine = prm.get<bool>("refine");
  bool coarsen = prm.get<bool>("coarsen");
  bool filter = prm.get<bool>("filter");

  bool frozen_fields = prm.get<bool>("frozen_fields");
  bool local_dt = prm.get<bool>("local_dt");
  //double dl_max = prm.dl_max;

  std::cout << "Creating folders..." << std::endl;

  // --check must not leave anything behind
  const bool dry_run = prm.check_only();

  std::string folder = intp->get_folder();
  std::string rwfolder = folder + "/RandomWalkers/";
  if (!dry_run)
    create_folder(rwfolder);
  std::string newfolder;
  if (prm.get<std::string>("restart_folder") != ""){
    newfolder = prm.get<std::string>("folder");
  }
  else {
    newfolder = get_newfoldername(rwfolder, prm);
    if (!dry_run)
      create_folder(newfolder);
  }
  newfolder = newfolder + "" + "0" + "/";
  std::string posfolder = newfolder + "Positions/";
  std::string checkpointsfolder = newfolder + "Checkpoints/";
  //std::string histfolder = newfolder + "Histograms/";
  if (!dry_run){
    create_folder(newfolder);
    create_folder(posfolder);
    create_folder(checkpointsfolder);
    //create_folder(histfolder);
  }
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
    if (key[0] == "uniform" || key[0] == "point"){
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
    set_initial_state(init_state, intp, prm, gens[0]);
    mesh.load_initial_state(init_state, prm);
  }

  mesh.compute_maps();

  // Everything that reads parameters has now been constructed, so stop here
  if (dry_run){
      std::cout << "Check OK: " << ps.N() << " particles, dim = " << mesh.dim() << std::endl;
    return 0;
  }

  // Initial refinement
  if (refine && !prm.get<bool>("inject") && mesh.dim() > 0){
    std::cout << "Initial refinement" << std::endl;
    Uint n_add = mesh.refine();
    if (prm.get<bool>("verbose"))
      std::cout << "Added " << n_add << " edges." << std::endl;
  }
  if (coarsen && !prm.get<bool>("inject") && mesh.dim() > 0){
    std::cout << "Initial coarsening" << std::endl;
    Uint n_rem = mesh.coarsen();
    if (prm.get<bool>("verbose"))
      std::cout << "Removed " << n_rem << " edges." << std::endl;
  }

  /*compute_interior_prop(interior_ang, mixed_areas, face_normals,
                        faces, edges, edge2faces, ps);
  compute_mean_curv(faces, edges,
                    edge2faces, node2edges,
                    ps, interior_ang, mixed_areas, face_normals);
  */
  mesh.compute_interior();

  //Vector3d dx_rw;
  int it = 0;
  //double dt2 = dt*dt;

  double sqrt2Dmdt = sqrt(2*Dm*dt);
  if (prm.get<bool>("verbose")){
    print_param("sqrt(2*Dm*dt)", sqrt2Dmdt);
    print_param("U*dt         ", prm.get<double>("U")*dt);
  }

  double t = t0;
  if (prm.get<std::string>("restart_folder") != ""){
    t = prm.get<double>("t");
  }

  prm.dump(newfolder, t);

  // Should not be taken from parameters
  //Uint n_accepted = prm.n_accepted;
  //Uint n_declined = prm.n_declined;

  std::string h5fname = newfolder + "/data_from_t" + std::to_string(t) + ".h5";
  H5::H5File h5f(h5fname.c_str(), H5F_ACC_TRUNC);
  //h5f->openFile(h5fname.c_str(), H5F_ACC_TRUNC);
  //H5wrap h5file();
  //h5file.open(h5fname, "w");

  Uint int_stat_intv = int(prm.get<double>("stat_intv")/dt);
  Uint int_dump_intv = int(prm.get<double>("dump_intv")/dt);
  Uint int_checkpoint_intv = int(prm.get<double>("checkpoint_intv")/dt);
  Uint int_chunk_intv = int_dump_intv*prm.get<int>("dump_chunk_size");
  Uint int_refine_intv = int(prm.get<double>("refine_intv")/dt);
  Uint int_coarsen_intv = int(prm.get<double>("coarsen_intv")/dt);
  //Uint int_hist_intv = int_stat_intv*prm.hist_chunk_size;
  Uint int_inject_intv = int(prm.get<double>("inject_intv")/dt);
  Uint int_filter_intv = int(prm.get<double>("filter_intv")/dt);
  //Uint int_resize_intv = int(prm.resize_intv/dt);
  Uint int_tau_intv = int(prm.get<double>("tau_intv")/dt);

  std::map<std::string, bool> output_fields;
  output_fields["u"] = !prm.get<bool>("minimal_output");
  output_fields["c"] = !prm.get<bool>("minimal_output") || local_dt;
  output_fields["p"] = !prm.get<bool>("minimal_output") && prm.get<bool>("output_all_props");
  output_fields["rho"] = !prm.get<bool>("minimal_output"); // && prm.output_all_props;   
  output_fields["H"] = !prm.get<bool>("minimal_output") && mesh.dim() > 0;
  output_fields["n"] = !prm.get<bool>("minimal_output") && mesh.dim() > 1;
  output_fields["t_loc"] = local_dt;
  output_fields["tau"] = prm.get<bool>("integrate_tau");

  if (local_dt && !frozen_fields){
    std::cout << "Error: local_dt=true requires the use of frozen_fields=true!" << std::endl;
    exit(1);
  }
  else if (local_dt && Dm > 0.0){
    std::cout << "Error: local_dt=true requires the use of Dm=0.0!" << std::endl;
  }
  else if (local_dt){
    std::cout << "Note: Using local time steps. Time t should now be considered only as a parametrizing variable." << std::endl;
  }

  bool any_exit_plane = (prm.get<std::string>("exit_plane") == "x" || prm.get<std::string>("exit_plane") == "y" || prm.get<std::string>("exit_plane") == "z") && prm.get<double>("Ln") > 0;
  int exit_dim = prm.get<std::string>("exit_plane") == "x" ? 0 : (prm.get<std::string>("exit_plane") == "y" ? 1 : 2);

  //std::string write_mode = prm.write_mode;

  std::ofstream statfile;
  {
    statfile.open(newfolder + "/tdata_from_t" + std::to_string(t) + ".dat");
    write_stats_header(statfile, mesh.dim());
  }
  std::ofstream declinedfile(newfolder + "/declinedpos_from_t" + std::to_string(t) + ".dat");

  // Simulation start
  std::clock_t clock_0 = std::clock();
  while (t < T + dt/2){
    if (!frozen_fields)
      intp->update(t);

    // Injection
    if (prm.get<bool>("inject") && it > 0 && it % int_inject_intv == 0 && t <= prm.get<double>("T_inject")){
      mesh.inject();
    }
    // Curvature computation
    if ((refine && it % int_refine_intv == 0) || (coarsen && it % int_coarsen_intv == 0) || it % int_dump_intv == 0){
      mesh.compute_interior();
    }

    // Refinement
    if (refine && it % int_refine_intv == 0 && it > 0){
      Uint n_add = mesh.refine();
      /*Uint n_add = refinement(faces, edges, edge2faces, node2edges, edges_inlet,
                              ps, ds_max,
                              prm.curv_refine_factor, prm.cut_if_stuck);*/
      if (prm.get<bool>("verbose"))
        std::cout << "Added " << n_add << " edges." << std::endl;
    }
    // Coarsening
    if (coarsen && it % int_coarsen_intv == 0){
      Uint n_rem = mesh.coarsen();
      /*Uint n_rem = coarsening(faces, edges,
                              edge2faces, node2edges,
                              edges_inlet, nodes_inlet,
                              ps, ds_min,
                              prm.curv_refine_factor);*/
      if (prm.get<bool>("verbose"))
        std::cout << "Removed " << n_rem << " edges." << std::endl;
    }
    // Filtering
    if (filter && it % int_filter_intv == 0){
      bool filtered = mesh.filter();
      /*bool filtered = filtering(faces, edges,
                                edge2faces, node2edges,
                                ps, prm.filter_target);*/
      if (prm.get<bool>("verbose") && filtered)
        std::cout << "Filtered edges." << std::endl;
    }
    // Resizing
    /*if (resize && it % int_resize_intv == 0){
      bool resized = mesh.resize();
      if (prm.verbose && resized)
        std::cout << "Resized edges." << std::endl;
    }*/
    // Removal
    if (any_exit_plane && it % int_filter_intv == 0){
      Uint n_rem = mesh.remove_beyond(exit_dim, prm.get<double>("Ln"));
      if (prm.get<bool>("verbose"))
        std::cout << "Removed " << n_rem << " nodes that were beyond." << std::endl;
    }

    // Tau integration
    if (prm.get<bool>("integrate_tau") && it % int_tau_intv == 0){
      mesh.integrate_tau(dt * int_tau_intv, prm.get<double>("tau_max"));
    }

    // Update fields if needed
    if (it % int_dump_intv == 0 || it % int_stat_intv == 0){
      ps.update_fields(t, output_fields);
    }

    // Statistics
    if (it % int_stat_intv == 0){
      std::cout << "Time = " << t << std::endl;
      mesh.write_statistics(statfile, t, prm.get<double>("ds_max"), *integrator);
    }

    // Checkpoint
    if (it % int_checkpoint_intv == 0){
      mesh.write_checkpoint(checkpointsfolder, t, prm);
    }

    // Dump detailed data
    if (it % int_dump_intv == 0){
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

    //if (outside_nodes.size() > 0 && prm.verbose)
    //  std::cout << "Some nodes are outside." << std::endl;

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

  //delete intp;
  return 0;
}
