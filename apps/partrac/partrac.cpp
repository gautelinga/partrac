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
  //double dl_max = prm.dl_max;

  std::cout << "Creating folders..." << std::endl;

  // --check must not leave anything behind
  const bool dry_run = prm.check_only();

  std::string folder = intp->get_folder();
  RunFolders out = make_run_folders(folder, "RandomWalkers", prm, dry_run ? DryRun : DefaultLayout);
  const std::string& newfolder = out.run;
  const std::string& checkpointsfolder = out.checkpoints;

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

  if (prm.get<bool>("inject"))
    std::cout << "Injection activated!" << std::endl;

  const bool restarting = prm.get<std::string>("restart_folder") != "";

  if (restarting){
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

  // Initial refinement, not on a restart
  if (refine && !restarting && !prm.get<bool>("inject") && mesh.dim() > 0){
    std::cout << "Initial refinement" << std::endl;
    Uint n_add = mesh.refine();
    if (prm.get<bool>("verbose"))
      std::cout << "Added " << n_add << " edges." << std::endl;
  }
  if (coarsen && !restarting && !prm.get<bool>("inject") && mesh.dim() > 0){
    std::cout << "Initial coarsening" << std::endl;
    Uint n_rem = mesh.coarsen(true);
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
  if (restarting){
    t = prm.get<double>("t");
    // Resume the step count, so the intervals keep their phase
    it = static_cast<int>(prm.get<Uint>("it"));
  }

  prm.dump(newfolder, t);

  // Should not be taken from parameters
  //Uint n_accepted = prm.n_accepted;
  //Uint n_declined = prm.n_declined;

  std::string h5fname = newfolder + "/data_from_t" + std::to_string(t) + ".h5";
  // no dump file at all when dumping is off
  H5::H5File h5f;
  if (prm.get<double>("dump_intv") > 0.){
    { H5::H5File create(h5fname.c_str(), H5F_ACC_TRUNC); }
    h5f.openFile(h5fname.c_str(), H5F_ACC_RDWR);
  }
  //h5f->openFile(h5fname.c_str(), H5F_ACC_TRUNC);
  //H5wrap h5file();
  //h5file.open(h5fname, "w");

  const double chunk_intv = prm.get<double>("dump_intv")*prm.get<int>("dump_chunk_size");
  //Uint int_hist_intv = int_stat_intv*prm.hist_chunk_size;
  //Uint int_resize_intv = int(prm.resize_intv/dt);

  std::map<std::string, bool> output_fields;
  output_fields["u"] = !prm.get<bool>("minimal_output");
  output_fields["c"] = !prm.get<bool>("minimal_output");
  output_fields["p"] = !prm.get<bool>("minimal_output") && prm.get<bool>("output_all_props");
  output_fields["rho"] = !prm.get<bool>("minimal_output"); // && prm.output_all_props;   
  // H and n are only computed with the curvature
  output_fields["H"] = !prm.get<bool>("minimal_output") && mesh.dim() > 0 && mesh.computes_curvature();
  output_fields["n"] = !prm.get<bool>("minimal_output") && mesh.dim() > 1 && mesh.computes_curvature();
  output_fields["tau"] = prm.get<bool>("integrate_tau");

  // Coarsening runs on every step its interval names, whether or not coarsening
  // was asked for -- with it off the threshold drops to what is numerically
  // zero. Refinement is what raises a zero-length median, so with coarsening
  // off the cleanup follows the refinement interval: coarsen_intv is one such
  // a run had no reason to set, and its default would leave the mesh degenerate
  const double coarsen_intv = coarsen ? prm.get<double>("coarsen_intv")
                                      : prm.get<double>("refine_intv");

  bool any_exit_plane = (prm.get<std::string>("exit_plane") == "x" || prm.get<std::string>("exit_plane") == "y" || prm.get<std::string>("exit_plane") == "z") && prm.get<double>("Ln") > 0;
  int exit_dim = prm.get<std::string>("exit_plane") == "x" ? 0 : (prm.get<std::string>("exit_plane") == "y" ? 1 : 2);

  //std::string write_mode = prm.write_mode;

  std::ofstream statfile;
  if (prm.get<double>("stat_intv") > 0.){
    statfile.open(newfolder + "/tdata_from_t" + std::to_string(t) + ".dat");
    write_stats_header(statfile, mesh.stats_header_columns(prm.get<double>("ds_max")));
  }

  // Simulation start
  std::clock_t clock_0 = std::clock();
  while (t < T + dt/2){
    if (!frozen_fields)
      intp->update(t);

    // Injection
    if (prm.get<bool>("inject") && it > 0 && at_interval(it, prm.get<double>("inject_intv"), dt) && t <= prm.get<double>("T_inject")){
      mesh.inject();
    }
    // Curvature computation
    if ((refine && at_interval(it, prm.get<double>("refine_intv"), dt)) || (coarsen && at_interval(it, prm.get<double>("coarsen_intv"), dt)) || at_interval(it, prm.get<double>("dump_intv"), dt)){
      mesh.compute_interior();
    }

    // Refinement
    if (refine && at_interval(it, prm.get<double>("refine_intv"), dt) && it > 0){
      Uint n_add = mesh.refine();
      /*Uint n_add = refinement(faces, edges, edge2faces, node2edges, edges_inlet,
                              ps, ds_max,
                              prm.curv_refine_factor, prm.cut_if_stuck);*/
      if (prm.get<bool>("verbose"))
        std::cout << "Added " << n_add << " edges." << std::endl;
    }
    // Coarsening
    if (at_interval(it, coarsen_intv, dt)){
      Uint n_rem = mesh.coarsen(coarsen);
      /*Uint n_rem = coarsening(faces, edges,
                              edge2faces, node2edges,
                              edges_inlet, nodes_inlet,
                              ps, ds_min,
                              prm.curv_refine_factor);*/
      if (prm.get<bool>("verbose"))
        std::cout << "Removed " << n_rem << " edges." << std::endl;
    }
    // Filtering
    if (filter && at_interval(it, prm.get<double>("filter_intv"), dt)){
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
    if (any_exit_plane && at_interval(it, prm.get<double>("filter_intv"), dt)){
      Uint n_rem = mesh.remove_beyond(exit_dim, prm.get<double>("Ln"));
      if (prm.get<bool>("verbose"))
        std::cout << "Removed " << n_rem << " nodes that were beyond." << std::endl;
    }

    // Update fields if needed
    if (at_interval(it, prm.get<double>("dump_intv"), dt) || at_interval(it, prm.get<double>("stat_intv"), dt)){
      ps.update_fields(t, output_fields);
    }

    // Statistics
    if (at_interval(it, prm.get<double>("stat_intv"), dt)){
      std::cout << "Time = " << t << std::endl;
      mesh.write_statistics(statfile, t, prm.get<double>("ds_max"), *integrator);
    }

    // Checkpoint
    if (at_interval(it, prm.get<double>("checkpoint_intv"), dt)){
      prm.set<Uint>("it", it);
      mesh.write_checkpoint(checkpointsfolder, t, prm);
    }

    // Dump detailed data
    if (at_interval(it, prm.get<double>("dump_intv"), dt)){
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

    // Tau integration, after the step whose interval it covers
    if (prm.get<bool>("integrate_tau") && at_interval(it + 1, prm.get<double>("tau_intv"), dt)){
      mesh.integrate_tau(dt * steps_per(prm.get<double>("tau_intv"), dt), prm.get<double>("tau_max"));
    }

    // Nodes that could not move: they stay, but where they pile up is useful
    if (outside_nodes.size() > 0 && prm.get<bool>("verbose")){
      Vector3d x_stuck = {0., 0., 0.};
      for (const Uint i : outside_nodes)
        x_stuck += ps.x(i);
      x_stuck /= outside_nodes.size();
      std::cout << outside_nodes.size() << " nodes could not move at t = " << t
                << ", centred on (" << x_stuck[0] << ", " << x_stuck[1] << ", "
                << x_stuck[2] << ")" << std::endl;
    }

    t += dt;
    it += 1;
  }
  std::clock_t clock_1 = std::clock();
  double duration = (clock_1-clock_0) / (double) CLOCKS_PER_SEC;
  std::cout << "Total simulation time: " << duration << " seconds" << std::endl;

  prm.set<Uint>("it", it);
  mesh.write_checkpoint(checkpointsfolder, t, prm);

  // Close files
  statfile.close();

  //delete intp;
  return 0;
}
