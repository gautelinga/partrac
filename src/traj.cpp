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
#include <ctime>

#include "io.hpp"
#include "utils.hpp"
#include "Parameters.hpp"
#include "Interpol.hpp"
#include "StructuredInterpol.hpp"
#include "AnalyticInterpol.hpp"
#ifdef USE_DOLFIN
#include "DolfInterpol.hpp"
#include "TetInterpol.hpp"
#include "TriangleInterpol.hpp"
#endif
#include "ParticleSet.hpp"
#include "Topology.hpp"
#include "Integrator.hpp"
#include "ExplicitIntegrator.hpp"
#include "distribute.hpp"
#include "stats.hpp"
#include "helpers.hpp"


int main(int argc, char* argv[])
{
  // Input parameters
  if (argc < 2) {
    std::cout << "Specify an input file." << std::endl;
    return 0;
  }
  Parameters prm(argc, argv);
  if (prm.restart_folder != ""){
    prm.parse_file(prm.restart_folder + "/Checkpoints/params.dat");
    prm.parse_cmd(argc, argv);
  }

  std::string infilename = std::string(argv[1]);

  Interpol *intp;
  std::string mode = prm.mode;
  if (mode == "analytic"){
    std::cout << "AnalyticInterpol initiated." << std::endl;
    intp = new AnalyticInterpol(infilename);
  }
  else if (mode == "unstructured" || mode == "fenics" || mode == "xdmf" ||
           mode == "tet" || mode == "triangle"){
#ifdef USE_DOLFIN
    if (mode == "tet")
      intp = new TetInterpol(infilename);
    else if (mode == "triangle")
      intp = new TriangleInterpol(infilename);
    else if (mode == "fenics")
      intp = new DolfInterpol(infilename);
    else if (mode == "xdmf"){
      std::cout << "XDMF format is not implemented yet." << std::endl;
      exit(0);
    }
    else {
      std::cout << "Mode should be 'fenics', 'tet' or 'triangle'." << std::endl;
      exit(1);
    }
#else
    std::cout << "You have to compile with PARTRAC_ENABLE_FENICS=ON." << std::endl;
    exit(0);
#endif
  }
  else if (mode == "structured" || mode == "lbm" || mode == "felbm"){
    intp = new StructuredInterpol(infilename);
  }
  else {
    std::cout << "Mode not supported." << std::endl;
    exit(0);
  }
  intp->set_U0(prm.U0);

  double Dm = prm.Dm;
  double dt = prm.dt;

  bool refine = prm.refine;
  bool coarsen = prm.coarsen;
  bool filter = prm.filter;

  bool frozen_fields = prm.frozen_fields;
  bool local_dt = prm.local_dt;
  double dl_max = prm.dl_max;

  std::string folder = intp->get_folder();
  std::string rwfolder = create_folder(folder + "/RandomWalkers/");
  std::string newfolder;
  if (prm.restart_folder != ""){
    newfolder = prm.folder;
  }
  else {
    std::string newfoldername = get_newfoldername(rwfolder, prm);
    newfolder = create_folder(newfoldername);
  }
  std::string posfolder = create_folder(newfolder + "Positions/");
  std::string checkpointsfolder = create_folder(newfolder + "Checkpoints/");
  std::string histfolder = create_folder(newfolder + "Histograms/");
  prm.folder = newfolder;

  prm.print();

  std::mt19937 gen;
  if (prm.random) {
    std::random_device rd;
    gen.seed(rd());
  }
  else {
    std::seed_seq rd{prm.seed};
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
    std::cout << "Cannot inject and filter at the same time (yet)." << std::endl;
    exit(0);
  }

  // Higher-order time integration?
  if (prm.int_order > 2){
    std::cout << "No support for such high temporal integration order." << std::endl;
    exit(0);
  }
  if (prm.interpolation_test > 0){
    std::cout << "Testing interpolation..." << std::endl;
    intp->set_int_order(2);
    test_interpolation(prm.interpolation_test, intp, newfolder, t0, gen);
  }
  intp->set_int_order(prm.int_order);

  if (frozen_fields)
    intp->update(prm.t_frozen);
  else
    intp->update(t0);

  ParticleSet ps(intp, prm.Nrw_max);
  //ps.Nrw = prm.Nrw;  // to remove?

  Topology mesh(ps, prm);

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
    std::vector<Vector3d> pos_init;
    pos_init = initial_positions(prm.init_mode,
                                 prm.init_weight,
                                 prm.Nrw,
                                 {prm.x0, prm.y0, prm.z0},
                                 prm.La, prm.Lb,
                                 prm.ds_init, t0,
                                 intp, gen, mesh.edges, mesh.faces);
    if (prm.clear_initial_edges){
      std::cout << "Clearing initial edges!" << std::endl;
      mesh.clear();
    }
    if (prm.inject){
      mesh.pos_inj = pos_init;
      mesh.edges_inj = mesh.edges;
      for (Uint i=0; i<pos_init.size(); ++i){
        mesh.nodes_inlet.push_back(i);
      }
      for (Uint i=0; i<mesh.edges.size(); ++i){
        mesh.edges_inlet.push_back(i);
      }
    }
    ps.add(pos_init, 0);
  }

  // std::cout << "aaa" << std::endl;
  /*add_particles(pos_init, intp,
                x_rw, u_rw, c_rw, tau_rw, rho_rw, p_rw, a_rw,
                U0, prm.restart_folder, prm.int_order, 0);*/
  // Nrw = pos_init.size();
  mesh.compute_maps();

  // Initial refinement
  if (refine && !prm.inject && mesh.dim() > 0){
    std::cout << "Initial refinement" << std::endl;
    /*Uint n_add = refinement(faces, edges,
                            edge2faces, node2edges,
                            edges_inlet,
                            ps, ds_max,
                            prm.curv_refine_factor,
                            prm.cut_if_stuck);*/
    Uint n_add = mesh.refine();

    std::cout << "Initial coarsening" << std::endl;
    /*Uint n_rem = coarsening(faces, edges,
                            edge2faces, node2edges,
                            edges_inlet, nodes_inlet,
                            ps, ds_min,
                            prm.curv_refine_factor);*/
    Uint n_rem = mesh.coarsen();

    if (prm.verbose)
      std::cout << "Added " << n_add << " edges and removed " << n_rem << " edges." << std::endl;
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
  if (prm.verbose){
    print_param("sqrt(2*Dm*dt)", sqrt2Dmdt);
    print_param("U*dt         ", prm.U0*dt);
  }

  double t = t0;
  if (prm.restart_folder != ""){
    t = prm.t;
  }

  Integrator* integrator;
  integrator = new ExplicitIntegrator(intp, Dm, prm.int_order, gen);
  ps.attach_integrator(integrator);

  prm.dump(newfolder, t);

  // Should not be taken from parameters
  Uint n_accepted = prm.n_accepted;
  Uint n_declined = prm.n_declined;

  std::string h5fname = newfolder + "/data_from_t" + std::to_string(t) + ".h5";
  H5FilePtr h5f = std::make_shared<H5::H5File>(h5fname.c_str(), H5F_ACC_TRUNC);
  // h5f->openFile(h5fname.c_str(), H5F_ACC_TRUNC);
  
  Uint int_stat_intv = int(prm.stat_intv/dt);
  Uint int_dump_intv = int(prm.dump_intv/dt);
  Uint int_checkpoint_intv = int(prm.checkpoint_intv/dt);
  Uint int_chunk_intv = int_dump_intv*prm.dump_chunk_size;
  Uint int_refine_intv = int(prm.refine_intv/dt);
  Uint int_coarsen_intv = int(prm.coarsen_intv/dt);
  Uint int_hist_intv = int_stat_intv*prm.hist_chunk_size;
  Uint int_inject_intv = int(prm.inject_intv/dt);
  Uint int_filter_intv = int(prm.filter_intv/dt);

  std::map<std::string, bool> output_fields;
  output_fields["u"] = !prm.minimal_output;
  output_fields["c"] = !prm.minimal_output || local_dt;
  output_fields["p"] = !prm.minimal_output && prm.output_all_props;
  output_fields["rho"] = !prm.minimal_output && prm.output_all_props;        
  output_fields["H"] = !prm.minimal_output && mesh.dim() > 0;
  output_fields["n"] = !prm.minimal_output && mesh.dim() > 1;
  output_fields["tau"] = local_dt;

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

  //std::string write_mode = prm.write_mode;

  std::ofstream statfile(newfolder + "/tdata_from_t" + std::to_string(t) + ".dat");
  write_stats_header(statfile, mesh.faces, mesh.edges);
  std::ofstream declinedfile(newfolder + "/declinedpos_from_t" + std::to_string(t) + ".dat");

  // Simulation start
  std::clock_t clock_0 = std::clock();
  while (t <= T){
    if (!frozen_fields)
      intp->update(t);
    // Statistics
    if (it % int_stat_intv == 0){
      std::cout << "Time = " << t << std::endl;
      bool do_dump_hist = (int_hist_intv > 0 && it % int_hist_intv == 0);
      write_stats(statfile, t, ps, mesh.faces, mesh.edges,
                  prm.ds_max,
                  do_dump_hist,
                  histfolder,
                  n_accepted, n_declined);
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
    // Dump detailed data
    if (it % int_dump_intv == 0){
      ps.update_fields(t, output_fields);
      {
        // Clear file if it exists, otherwise create
        if (int_chunk_intv > 0 && it % int_chunk_intv == 0 && it > 0){
          //h5f->close();
          h5fname = newfolder + "/data_from_t" + std::to_string(t) + ".h5";
          h5f = std::make_shared<H5::H5File>(h5fname.c_str(), H5F_ACC_TRUNC);
        }
        else {
          h5f->openFile(h5fname.c_str(), H5F_ACC_RDWR);
        }
        std::string groupname = std::to_string(t);
        h5f->createGroup(groupname + "/");

        /*if (faces.size() == 0){
          for (EdgesType::const_iterator edgeit = edges.begin();
               edgeit != edges.end(); ++edgeit){
            int inode = edgeit->first[0];
            int jnode = edgeit->first[1];
            double ds0 = edgeit->second;
            double ds = ps.dist(inode, jnode);

            if (ps.e_rw[inode] <= 0.) ps.e_rw[inode] = ds/ds0;
            else ps.e_rw[inode] = 0.5*(ps.e_rw[inode] + ds/ds0);
            if (ps.e_rw[jnode] <= 0.) ps.e_rw[jnode] = ds/ds0;
            else ps.e_rw[jnode] = 0.5*(ps.e_rw[jnode] + ds/ds0);
          }
        }*/
        ps.dump_hdf5(h5f, groupname, output_fields);
        mesh.dump_hdf5(h5f, groupname);

        h5f->close();
      }
    }
    /*
    if (!local_dt){
      for (Uint irw=0; irw < ps.Nrw; ++irw){
        dx_rw = ps.u_rw[irw]*dt;

        // Set elongation
        //ps.e_rw[irw] = 0.;

        // Second-order terms
        if (prm.int_order >= 2){
          dx_rw += 0.5*ps.a_rw[irw]*dt2;
        }
        if (Dm > 0.0){
          Vector3d eta = {rnd_normal(gen),
                          rnd_normal(gen),
                          rnd_normal(gen)};
          dx_rw += sqrt2Dmdt*eta;
        }
        intp->probe(ps.x_rw[irw]+dx_rw);
        if (intp->inside_domain()){
          ps.x_rw[irw] += dx_rw;
          ps.u_rw[irw] = intp->get_u();

          if ((it+1) % int_dump_intv == 0){
            ps.rho_rw[irw] = intp->get_rho();
            ps.p_rw[irw] = intp->get_p();
          }

          // Second-order terms
          if (prm.int_order >= 2){
            ps.a_rw[irw] = intp->get_Ju() + intp->get_a();
          }
          n_accepted++;
        }
        else {
          n_declined++;
          declinedfile << t << " "
                       << ps.x_rw[irw][0]+dx_rw[0] << " "
                       << ps.x_rw[irw][1]+dx_rw[1] << " "
                       << ps.x_rw[irw][2]+dx_rw[2] << std::endl;
        }
      }
    }
    else {
      double dtau;
      double u_abs;
      std::set<Uint> nodes_to_remove;
      for (Uint irw=0; irw < ps.Nrw; ++irw){
        u_abs = ps.u_rw[irw].norm();
        if (u_abs < prm.u_eps){
          //
          nodes_to_remove.insert(irw);
        }
        dtau = dl_max/(u_abs+prm.u_eps);
        ps.tau_rw[irw] += dtau;
        dx_rw = ps.u_rw[irw]*dtau;

        ps.e_rw[irw] = 0.; //elongation - why??

        if (prm.int_order >= 2){
          dx_rw += 0.5 * ps.a_rw[irw] * dtau * dtau;
        }
        intp->probe(ps.x_rw[irw]+dx_rw);
        if (intp->inside_domain()){
          ps.x_rw[irw] += dx_rw;
          ps.u_rw[irw] = intp->get_u();
          // other scalars?
          if (prm.int_order >= 2){
            ps.a_rw[irw] = intp->get_Ju() + intp->get_a();
          }
          n_accepted++;
        }
        else {
          n_declined++;
          //log??
        }
      }
      if (nodes_to_remove.size() > 0){
        if (!prm.cut_if_stuck){
          std::cout << "Node is stuck! Turn on cut_if_stuck to continue." << std::endl;
          exit(0);
        }
        bool do_output_all = prm.output_all_props && !prm.minimal_output && (it+1) % int_dump_intv == 0;
        std::vector<bool> node_isactive(ps.Nrw, true);
        for (std::set<Uint>::const_iterator sit = nodes_to_remove.begin();
             sit != nodes_to_remove.end(); ++sit){
          node_isactive[*sit] = false;
        }
        remove_nodes_safely(faces, edges,
                            edge2faces, node2edges,
                            edges_inlet, nodes_inlet,
                            node_isactive,
                            ps);
      }
    }
    */
    bool all_inside = ps.integrate(t, dt);
    
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
  //if (write_mode == "hdf5"){
  h5f->close();
  //}

  delete intp;
  return 0;
}
