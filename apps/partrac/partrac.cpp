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
#include "helpers.hpp"
#include "MPIwrap.hpp"

int main(int argc, char* argv[])
{
  MPIwrap mpi(argc, argv);

  if (mpi.rank() == 0)
    std::cout << "Initialized Partrac with " << mpi.size() << " processes." << std::endl;
  mpi.barrier();
  std::cout << "This is process " << mpi.rank() << " out of " << mpi.size() << "." << std::endl;
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

  std::cout << "Setting interpolator..." << std::endl;

  std::shared_ptr<Interpol> intp;
  set_interpolate_mode(intp, prm.mode, infilename);
  
  intp->set_U0(prm.U0);
  intp->set_int_order(prm.int_order);

  double Dm = prm.Dm;
  double dt = prm.dt;

  bool refine = prm.refine;
  bool coarsen = prm.coarsen;
  bool filter = prm.filter;

  bool frozen_fields = prm.frozen_fields;
  bool local_dt = prm.local_dt;
  //double dl_max = prm.dl_max;

  std::cout << "Creating folders..." << std::endl;

  std::string folder = intp->get_folder();
  std::string rwfolder = folder + "/RandomWalkers/"; 
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
  //std::string histfolder = newfolder + "Histograms/";
  //if (mpi.rank() == 0){ // Might change in the future!
  {
    create_folder(newfolder);
    create_folder(posfolder);
    create_folder(checkpointsfolder);
    //create_folder(histfolder);
  }
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

  /*
  if (prm.resize && (prm.refine || prm.coarsen || prm.filter)){
    if (prm.refine){
      if (mpi.rank() == 0)
        std::cout << "Cannot resize and refine at the same time." << std::endl;
    }
    if (prm.coarsen){
      if (mpi.rank() == 0)
        std::cout << "Cannot resize and coarsen at the same time." << std::endl;
    }
    if (prm.filter){
      if (mpi.rank() == 0)
        std::cout << "Cannot resize and filter at the same time." << std::endl;
    }
    exit(0);
  }*/

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
    if (key[0] == "uniform" || key[0] == "point"){
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
    set_initial_state(init_state, intp, mpi, prm, gen);
    mesh.load_initial_state(init_state);
  }

  mesh.compute_maps();

  // Initial refinement
  if (refine && !prm.inject && mesh.dim() > 0){
    std::cout << "Initial refinement" << std::endl;
    Uint n_add = mesh.refine();
    if (prm.verbose && mpi.rank() == 0)
      std::cout << "Added " << n_add << " edges." << std::endl;
  }
  if (coarsen && !prm.inject && mesh.dim() > 0){
    std::cout << "Initial coarsening" << std::endl;
    Uint n_rem = mesh.coarsen();
    if (prm.verbose && mpi.rank() == 0)
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
  if (prm.verbose && mpi.rank() == 0){
    print_param("sqrt(2*Dm*dt)", sqrt2Dmdt);
    print_param("U*dt         ", prm.U0*dt);
  }

  double t = t0;
  if (prm.restart_folder != ""){
    t = prm.t;
  }

  //if (mpi.rank() == 0)
  prm.dump(newfolder, t);

  // Should not be taken from parameters
  //Uint n_accepted = prm.n_accepted;
  //Uint n_declined = prm.n_declined;

  std::string h5fname = newfolder + "/data_from_t" + std::to_string(t) + ".h5";
  H5File h5f(h5fname.c_str(), H5F_ACC_TRUNC);
  //h5f->openFile(h5fname.c_str(), H5F_ACC_TRUNC);
  //H5wrap h5file(mpi);
  //h5file.open(h5fname, "w");

  Uint int_stat_intv = int(prm.stat_intv/dt);
  Uint int_dump_intv = int(prm.dump_intv/dt);
  Uint int_checkpoint_intv = int(prm.checkpoint_intv/dt);
  Uint int_chunk_intv = int_dump_intv*prm.dump_chunk_size;
  Uint int_refine_intv = int(prm.refine_intv/dt);
  Uint int_coarsen_intv = int(prm.coarsen_intv/dt);
  //Uint int_hist_intv = int_stat_intv*prm.hist_chunk_size;
  Uint int_inject_intv = int(prm.inject_intv/dt);
  Uint int_filter_intv = int(prm.filter_intv/dt);
  //Uint int_resize_intv = int(prm.resize_intv/dt);
  Uint int_tau_intv = int(prm.tau_intv/dt);

  std::map<std::string, bool> output_fields;
  output_fields["u"] = !prm.minimal_output;
  output_fields["c"] = !prm.minimal_output || local_dt;
  output_fields["p"] = !prm.minimal_output && prm.output_all_props;
  output_fields["rho"] = !prm.minimal_output && prm.output_all_props;        
  output_fields["H"] = !prm.minimal_output && mesh.dim() > 0;
  output_fields["n"] = !prm.minimal_output && mesh.dim() > 1;
  output_fields["t_loc"] = local_dt;
  output_fields["tau"] = prm.integrate_tau;

  if (local_dt && !frozen_fields){
    std::cout << "Error: local_dt=true requires the use of frozen_fields=true!" << std::endl;
    exit(0);
  }
  else if (local_dt && Dm > 0.0){
    std::cout << "Error: local_dt=true requires the use of Dm=0.0!" << std::endl;
  }
  else if (local_dt){
    std::cout << "Note: Using local time steps. Time t should now be considered only as a parametrizing variable." << std::endl;
  }

  bool any_exit_plane = (prm.exit_plane == "x" || prm.exit_plane == "y" || prm.exit_plane == "z") && prm.Ln > 0;
  int exit_dim = prm.exit_plane == "x" ? 0 : (prm.exit_plane == "y" ? 1 : 2);

  //std::string write_mode = prm.write_mode;

  std::ofstream statfile;
  //if (mpi.rank() == 0){
  {
    statfile.open(newfolder + "/tdata_from_t" + std::to_string(t) + ".dat");
    write_stats_header(mpi, statfile, mesh.dim());
  }
  std::ofstream declinedfile(newfolder + "/declinedpos_from_t" + std::to_string(t) + ".dat");

  // Simulation start
  std::clock_t clock_0 = std::clock();
  while (t < T + dt/2){
    if (!frozen_fields)
      intp->update(t);
    // Update fields if needed
    if (it % int_dump_intv == 0 || it % int_stat_intv == 0){
      ps.update_fields(t, output_fields);
    }

    // Statistics
    if (it % int_stat_intv == 0){
      std::cout << "Time = " << t << std::endl;
      mesh.write_statistics(statfile, t, prm.ds_max, *integrator);
    }
    // Checkpoint
    if (it % int_checkpoint_intv == 0){
      mesh.write_checkpoint(checkpointsfolder, t, prm);
    }
    // Injection
    if (prm.inject && it > 0 && it % int_inject_intv == 0 && t <= prm.T_inject){
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
      if (prm.verbose)
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
      if (prm.verbose)
        std::cout << "Removed " << n_rem << " edges." << std::endl;
    }
    // Filtering
    if (filter && it % int_filter_intv == 0){
      bool filtered = mesh.filter();
      /*bool filtered = filtering(faces, edges,
                                edge2faces, node2edges,
                                ps, prm.filter_target);*/
      if (prm.verbose && filtered)
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
      Uint n_rem = mesh.remove_beyond(exit_dim, prm.Ln);
      if (prm.verbose)
        std::cout << "Removed " << n_rem << " nodes that were beyond." << std::endl;
    }
    // Tau integration
    if (prm.integrate_tau && it % int_tau_intv == 0){
      mesh.integrate_tau(dt * int_tau_intv, prm.tau_max);
    }

    // Dump detailed data
    if (it % int_dump_intv == 0){
      std::string groupname = std::to_string(t);
      //if (mpi.rank() == 0) {
      //{
        // Clear file if it exists, otherwise create
      if (int_chunk_intv > 0 && it % int_chunk_intv == 0 && it > 0){
        h5fname = newfolder + "/data_from_t" + std::to_string(t) + ".h5";
        // h5f = std::make_shared<H5::H5File>(h5fname.c_str(), H5F_ACC_TRUNC);
        h5f.openFile(h5fname.c_str(), H5F_ACC_TRUNC);
        //h5f->openFile(h5fname.c_str(), H5F_ACC_TRUNC);
        //h5file.open(h5fname, "w");
      }
      else {
        //h5f->openFile(h5fname.c_str(), H5F_ACC_RDWR);
        h5f.openFile(h5fname.c_str(), H5F_ACC_RDWR);
        //h5file.open(h5fname, "a");
      }
      //h5f->createGroup(groupname + "/");
      h5f.createGroup(groupname + "/");
      //}

      mesh.dump_hdf5(h5f, groupname, output_fields);

      h5f.close();
      //h5f->close();
      //h5file.close();
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
  //if (mpi.rank() == 0){
  statfile.close();
  declinedfile.close();

  //delete intp;
  return 0;
}
