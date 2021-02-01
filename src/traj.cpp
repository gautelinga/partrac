#include <iostream>
#include "io.hpp"
#include "utils.hpp"
#include "Parameters.hpp"
#include "Interpol.hpp"
#include "StructuredInterpol.hpp"
#include "AnalyticInterpol.hpp"
#ifdef USE_DOLFIN
#include "DolfInterpol.hpp"
#include "TetInterpol.hpp"
#endif
#include "distribute.hpp"
#include "mesh.hpp"
#include "stats.hpp"
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

using namespace H5;


int main(int argc, char* argv[]){
  // Input parameters
  if (argc < 2) {
    std::cout << "Specify an input timestamps file." << std::endl;
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
  else if (mode == "unstructured" || mode == "fenics" || mode == "tet"){
#ifdef USE_DOLFIN
    std::cout << "FEniCS/XDMF format is not implemented yet." << std::endl;
    if (mode == "tet")
      intp = new TetInterpol(infilename);
    else if (mode == "fenics")
      intp = new DolfInterpol(infilename);
    else {
      std::cout << "mode should be fenics or tet" << std::endl;
      exit(1);
    }
#else
    std::cout << "You have to compile with PARTRAC_ENABLE_FENICS=ON." << std::endl;
    exit(0);
#endif
  }
  else if (mode == "structured" || mode == "lbm"){
    intp = new StructuredInterpol(infilename);
  }
  else {
    std::cout << "Mode not supported." << std::endl;
    exit(0);
  }

  double Dm = prm.Dm;
  double dt = prm.dt;
  Uint Nrw = prm.Nrw;
  double U0 = prm.U0;

  bool refine = prm.refine;
  bool coarsen = prm.coarsen;
  double ds_max = prm.ds_max;
  double ds_min = prm.ds_min;
  Uint Nrw_max = prm.Nrw_max;

  std::string folder = intp->get_folder();
  std::string rwfolder = create_folder(folder + "/RandomWalkers/");
  std::string newfolder;
  if (prm.restart_folder != ""){
    newfolder = prm.folder;
  }
  else {
    std::ostringstream ss_Dm, ss_dt, ss_Nrw, ss_seed;
    ss_Dm << std::scientific << std::setprecision(7) << Dm;
    ss_dt << std::scientific << std::setprecision(7) << dt;
    ss_Nrw << Nrw;
    ss_seed << prm.seed;
    newfolder = create_folder(rwfolder +
                              "/Dm" + ss_Dm.str() + // "_U" + std::to_string(prm.U0) +
                              "_dt" + ss_dt.str() +
                              "_Nrw" + ss_Nrw.str() +
                              "_seed" + ss_seed.str() +
                              "_tag" + prm.tag + 
                              "/");
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
  
  double Lx = intp->get_Lx();
  double Ly = intp->get_Ly();
  double Lz = intp->get_Lz();
  prm.Lx = Lx;
  prm.Ly = Ly;
  prm.Lz = Lz;
  //prm.nx = intp->get_nx();
  //prm.ny = intp->get_ny();
  //prm.nz = intp->get_nz();

  double U02 = U0*U0;

  double t0 = std::max(intp->get_t_min(), prm.t0);
  double T = std::min(intp->get_t_max(), prm.T);
  prm.t0 = t0;
  prm.T = T;

  if (prm.interpolation_test > 0){
    std::cout << "Testing interpolation..." << std::endl;
    intp->set_int_order(2);
    test_interpolation(prm.interpolation_test, intp, newfolder, t0, gen);
  }
  intp->set_int_order(prm.int_order);

  std::normal_distribution<double> rnd_normal(0.0, 1.0);

  std::vector<Vector3d> x_rw(Nrw_max);

  std::vector<double> c_rw(Nrw_max);
  std::vector<double> e_rw(Nrw_max);
  std::vector<double> H_rw(Nrw_max);

  std::vector<Vector3d> n_rw(Nrw_max);

  std::vector<Vector3d> u_rw(Nrw_max);

  std::vector<double> rho_rw(Nrw_max);
  std::vector<double> p_rw(Nrw_max);

  if (prm.inject && prm.filter){
    std::cout << "Cannot inject and filter at the same time (yet)." << std::endl;
    exit(0);
  }

  // Second-order terms
  std::vector<Vector3d> a_rw(Nrw_max);
  if (prm.int_order > 2){
    std::cout << "No support for such high temporal integration order." << std::endl;
    exit(0);
  }

  std::vector<Vector3d> pos_init;
  EdgesType edges;
  FacesType faces;

  std::vector<Vector3d> pos_inj;
  EdgesType edges_inj;
  EdgesListType edges_inlet;
  NodesListType nodes_inlet;

  if (prm.restart_folder != ""){
    std::string posfile = prm.restart_folder + "/Checkpoints/positions.pos";
    load_positions(posfile, pos_init, Nrw);
    std::string facefile = prm.restart_folder + "/Checkpoints/faces.face";
    load_faces(facefile, faces);
    std::string edgefile = prm.restart_folder + "/Checkpoints/edges.edge";
    load_edges(edgefile, edges);
    std::string colfile = prm.restart_folder + "/Checkpoints/colors.col";
    load_colors(colfile, c_rw, Nrw);
    if (prm.inject){
      std::string posinjfile = prm.restart_folder + "/Checkpoints/positions_inj.pos";
      load_positions(posinjfile, pos_inj);
      std::string edgesinjfile = prm.restart_folder + "/Checkpoints/edges_inj.edge";
      load_edges(edgesinjfile, edges_inj);
      std::string edgesinletfile = prm.restart_folder + "/Checkpoints/edges_inlet.list";
      std::string nodesinletfile = prm.restart_folder + "/Checkpoints/nodes_inlet.list";
      load_list(edgesinletfile, edges_inlet);
      load_list(nodesinletfile, nodes_inlet);
    }
  }
  else {
    pos_init = initial_positions(prm.init_mode,
                                 prm.init_weight,
                                 Nrw,
                                 {prm.x0, prm.y0, prm.z0},
                                 prm.La, prm.Lb,
                                 prm.ds_init, t0,
                                 intp, gen, edges, faces);
    if (prm.clear_initial_edges){
      std::cout << "Clearing initial edges!" << std::endl;
      edges.clear();
      faces.clear();
    }
    if (prm.inject){
     pos_inj = pos_init;
      edges_inj = edges;
      for (Uint i=0; i<pos_init.size(); ++i){
        nodes_inlet.push_back(i);
      }
      for (Uint i=0; i<edges.size(); ++i){
        edges_inlet.push_back(i);
      }
    }
  }

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
  // std::cout << "aaa" << std::endl;

  // Compute edge2faces map
  Edge2FacesType edge2faces;
  compute_edge2faces(edge2faces, faces, edges);

  // Compute node2edges map
  Node2EdgesType node2edges;
  compute_node2edges(node2edges, edges, Nrw);

  intp->update(t0);
  add_particles(pos_init, intp,
                x_rw, u_rw, c_rw, rho_rw, p_rw, a_rw,
                U0, prm.restart_folder, prm.int_order, 0);
  // Nrw = pos_init.size();

  // Initial refinement
  if (refine && !prm.inject && edges.size() > 0){
    std::cout << "Initial refinement" << std::endl;

    bool do_output_all = prm.output_all_props && !prm.minimal_output;
    Uint n_add = refinement(faces, edges,
                            edge2faces, node2edges,
                            edges_inlet,
                            x_rw,
                            u_rw,
                            rho_rw, p_rw, c_rw,
                            H_rw, n_rw,
                            a_rw,
                            Nrw, Nrw_max, ds_max,
                            do_output_all, intp,
                            U0, prm.int_order,
                            prm.curv_refine_factor);

    Uint n_rem = coarsening(faces, edges,
                            edge2faces, node2edges,
                            edges_inlet, nodes_inlet,
                            x_rw,
                            u_rw,
                            rho_rw, p_rw, c_rw,
                            H_rw, n_rw,
                            a_rw,
                            Nrw, ds_min,
                            do_output_all, intp,
                            U0, prm.int_order,
                            prm.curv_refine_factor);
    if (prm.verbose)
      std::cout << "Added " << n_add << " edges and removed " << n_rem << " edges." << std::endl;

    // std::vector<bool> face_isactive(faces.size(), true);
    // for (Uint i=0; i<faces.size(); ++i){
    //   face_isactive[i] = x_rw[edges[faces[i].first[0]].first[0]] < 30.0;
    // }
    // remove_faces(faces, face_isactive);
    // remove_unused_edges(faces, edges);
    // remove_unused_nodes(edges, x_rw, y_rw, z_rw, Nrw);
    // compute_edge2faces(edge2faces, faces, edges);
  }
  // Initial curvature?
  InteriorAnglesType interior_ang;
  std::vector<double> mixed_areas;
  std::vector<Vector3d> face_normals;
  compute_interior_prop(interior_ang, mixed_areas, face_normals,
                        faces, edges, edge2faces,
                        x_rw, Nrw);
  compute_mean_curv(H_rw, n_rw,
                    faces, edges,
                    edge2faces, node2edges,
                    x_rw, Nrw,
                    interior_ang, mixed_areas, face_normals);

  Vector3d dx_rw;
  int it = 0;
  double dt2 = dt*dt;

  double sqrt2Dmdt = sqrt(2*Dm*dt);
  if (prm.verbose){
    print_param("sqrt(2*Dm*dt)", sqrt2Dmdt);
    print_param("U*dt         ", U0*dt);
  }

  double t = t0;
  if (prm.restart_folder != ""){
    t = prm.t;
  }

  prm.dump(newfolder, t);

  unsigned long int n_accepted = prm.n_accepted;
  unsigned long int n_declined = prm.n_declined;

  H5File* h5f = new H5File(newfolder + "/data_from_t" + std::to_string(t) + ".h5", H5F_ACC_TRUNC);

  Uint int_stat_intv = int(prm.stat_intv/dt);
  Uint int_dump_intv = int(prm.dump_intv/dt);
  Uint int_checkpoint_intv = int(prm.checkpoint_intv/dt);
  Uint int_chunk_intv = int_dump_intv*prm.dump_chunk_size;
  Uint int_refine_intv = int(prm.refine_intv/dt);
  Uint int_coarsen_intv = int(prm.coarsen_intv/dt);
  Uint int_hist_intv = int_stat_intv*prm.hist_chunk_size;
  Uint int_inject_intv = int(prm.inject_intv/dt);

  bool filter = prm.filter;
  Uint int_filter_intv = int(prm.filter_intv/dt);

  std::string write_mode = prm.write_mode;

  std::ofstream statfile(newfolder + "/tdata_from_t" + std::to_string(t) + ".dat");
  write_stats_header(statfile, faces, edges);
  std::ofstream declinedfile(newfolder + "/declinedpos_from_t" + std::to_string(t) + ".dat");
  while (t <= T){
    intp->update(t);
    // Statistics
    if (it % int_stat_intv == 0){
      std::cout << "Time = " << t << std::endl;
      bool do_dump_hist = (int_hist_intv > 0 && it % int_hist_intv == 0);
      write_stats(statfile, t,
                  x_rw,
                  u_rw,
                  Nrw, faces, edges,
                  ds_max,
                  do_dump_hist,
                  histfolder,
                  n_accepted, n_declined);
    }
    // Checkpoint
    if (it % int_checkpoint_intv == 0){
       prm.t = t;
       prm.n_accepted = n_accepted;
       prm.n_declined = n_declined;
       prm.Nrw = Nrw;
       prm.dump(checkpointsfolder);
       dump_positions(checkpointsfolder + "/positions.pos", x_rw, Nrw);
       dump_faces(checkpointsfolder + "/faces.face", faces);
       dump_edges(checkpointsfolder + "/edges.edge", edges);
       dump_colors(checkpointsfolder + "/colors.col", c_rw, Nrw);
       if (prm.inject){
         dump_positions(checkpointsfolder + "/positions_inj.pos", pos_inj);
         dump_edges(checkpointsfolder + "/edges_inj.edge", edges_inj);
         dump_list(checkpointsfolder + "/edges_inlet.list", edges_inlet);
         dump_list(checkpointsfolder + "/nodes_inlet.list", nodes_inlet);
       }
    }
   // Injection
    if (prm.inject && it > 0 && it % int_inject_intv == 0 && t <= prm.T_inject){
      if (Nrw + pos_inj.size() <= Nrw_max){
        Uint irw0 = Nrw;
        Uint iedge0 = edges.size();
        Uint iface0 = faces.size();
        add_particles(pos_inj, intp,
                      x_rw, u_rw, c_rw, rho_rw, p_rw, a_rw,
                      U0, prm.restart_folder, prm.int_order, irw0);
        for (Uint irw=0; irw < pos_inj.size(); ++irw){
          node2edges.push_back({});
        }
        Nrw += pos_inj.size();
        if (prm.verbose){
          std::cout << "Added " << pos_inj.size() << " nodes." << std::endl;
          /*
          std::cout << "at ..." << std::endl;
          for (Uint irw=irw0; irw<Nrw; ++irw){
            std::cout << x_rw[irw][0] << " " << x_rw[irw][1] << " " << x_rw[irw][2] << std::endl;
          }
          */
        }
        if (prm.inject_edges){
          Uint iedge = iedge0;
          for (EdgesType::const_iterator edgeit = edges_inj.begin();
               edgeit != edges_inj.end(); ++edgeit){
            Uint inode = irw0 + edgeit->first[0];
            Uint jnode = irw0 + edgeit->first[1];
            edges.push_back({{inode, jnode}, edgeit->second});
            edge2faces.push_back({});
            node2edges[inode].push_back(iedge);
            node2edges[jnode].push_back(iedge);
            ++iedge;
          }
          // Straight edges
          for (Uint irw=0; irw < pos_inj.size(); ++irw){
            Uint inode = irw0+irw; // New node
            Uint jnode = *std::next(nodes_inlet.begin(), irw); // Old node. do better
            double ds0 = dist(inode, jnode, x_rw);
            edges.push_back({{inode, jnode}, ds0});
            edge2faces.push_back({});
            node2edges[inode].push_back(iedge);
            node2edges[jnode].push_back(iedge);
            ++iedge;
          }
          // Diagonal edges
          Uint iface = iface0;
          for (Uint jedge=0; jedge < edges_inj.size(); ++jedge){
            Uint kedge = *std::next(edges_inlet.begin(), jedge); // do better!
            Uint ledge = iedge0+jedge;
            Uint inode = edges[ledge].first[0]; // New edge
            Uint lnode = edges[ledge].first[1]; //
            Uint jnode = edges[kedge].first[1]; // Old edge
            Uint knode = edges[kedge].first[0]; //
            double ds0 = dist(inode, jnode, x_rw);
            edges.push_back({{inode, jnode}, ds0});
            edge2faces.push_back({});
            node2edges[inode].push_back(iedge);
            node2edges[jnode].push_back(iedge);
            long double dA01 = area(iedge, kedge, x_rw, edges);
            Uint medge = get_common_entry(node2edges[knode], node2edges[inode]);
            faces.push_back({{iedge, kedge, medge}, dA01});
            edge2faces[iedge].push_back(iface);
            edge2faces[kedge].push_back(iface);
            edge2faces[medge].push_back(iface);
            ++iface;
            long double dA02 = area(iedge, ledge, x_rw, edges);
            Uint nedge = get_common_entry(node2edges[lnode], node2edges[jnode]);
            faces.push_back({{iedge, ledge, nedge}, dA02});
            edge2faces[iedge].push_back(iface);
            edge2faces[ledge].push_back(iface);
            edge2faces[nedge].push_back(iface);
            ++iface;
            ++iedge;
          }
          assert(iedge == edges.size());
          if (prm.verbose){
            std::cout << "Added " << iedge-iedge0 << " edges." << std::endl;
            std::cout << "Added " << faces.size()-iface0 << " faces." << std::endl;
          }
          edges_inlet.clear();
          for (Uint jedge=0; jedge < edges_inj.size(); ++jedge){
            edges_inlet.push_back({iedge0+jedge});
          }
        }
        nodes_inlet.clear();
        for (Uint irw=0; irw < pos_inj.size(); ++irw){
          nodes_inlet.push_back(irw0+irw);
        }
      }
    }
    // Curvature computation
    if ((refine && it % int_refine_intv == 0) || // Consider other interval
        (coarsen && it % int_coarsen_intv == 0)){
      compute_interior_prop(interior_ang, mixed_areas, face_normals,
                            faces, edges, edge2faces,
                            x_rw, Nrw);
      compute_mean_curv(H_rw, n_rw, faces, edges,
                        edge2faces, node2edges,
                        x_rw, Nrw,
                        interior_ang, mixed_areas, face_normals);
    }

   // Refinement
    if (refine && it % int_refine_intv == 0 && it > 0){
      bool do_output_all = prm.output_all_props && !prm.minimal_output && it % int_dump_intv == 0;
      Uint n_add = refinement(faces, edges, edge2faces, node2edges, edges_inlet,
                              x_rw,
                              u_rw,
                              rho_rw, p_rw, c_rw,
                              H_rw, n_rw,
                              a_rw,
                              Nrw, Nrw_max, ds_max,
                              do_output_all, intp,
                              U0, prm.int_order,
                              prm.curv_refine_factor);
      if (prm.verbose)
        std::cout << "Added " << n_add << " edges." << std::endl;
    }
    // Coarsening
    if (coarsen && it % int_coarsen_intv == 0){
      bool do_output_all = prm.output_all_props && !prm.minimal_output && it % int_dump_intv == 0;
      Uint n_rem = coarsening(faces, edges,
                              edge2faces, node2edges,
                              edges_inlet, nodes_inlet,
                              x_rw,
                              u_rw,
                              rho_rw, p_rw, c_rw,
                              H_rw, n_rw,
                              a_rw,
                              Nrw, ds_min,
                              do_output_all, intp,
                              U0, prm.int_order,
                              prm.curv_refine_factor);
      if (prm.verbose)
        std::cout << "Removed " << n_rem << " edges." << std::endl;
    }
    // Filtering
    if (filter && it % int_filter_intv == 0){
      bool do_output_all = prm.output_all_props && !prm.minimal_output && it % int_dump_intv == 0;
      bool filtered = filtering(faces, edges,
                                edge2faces, node2edges,
                                x_rw, u_rw,
                                rho_rw, p_rw, c_rw,
                                H_rw, n_rw,
                                a_rw, Nrw,
                                do_output_all, prm.int_order, prm.filter_target);
      if (prm.verbose && filtered)
        std::cout << "Filtered edges." << std::endl;
    }
    // Dump detailed data
    if (it % int_dump_intv == 0){
      if (write_mode == "text"){
        std::string posfile = posfolder + "xy_t" + std::to_string(t) + ".pos";
        std::ofstream pos_out(posfile);
        posdata2txt(pos_out, x_rw, u_rw, Nrw);
        pos_out.close();
      }
      else if (write_mode == "hdf5"){
        // Clear file if it exists, otherwise create
        if (it % int_chunk_intv == 0 && it > 0){
          h5f->close();
          h5f = new H5File(newfolder + "/data_from_t" + std::to_string(t) + ".h5", H5F_ACC_TRUNC);
        }
        std::string groupname = std::to_string(t);
        h5f->createGroup(groupname + "/");

        if (faces.size() == 0){
          for (EdgesType::const_iterator edgeit = edges.begin();
               edgeit != edges.end(); ++edgeit){
            int inode = edgeit->first[0];
            int jnode = edgeit->first[1];
            double ds0 = edgeit->second;
            double ds = dist(inode, jnode, x_rw);

            if (e_rw[inode] <= 0.) e_rw[inode] = ds/ds0;
            else e_rw[inode] = 0.5*(e_rw[inode] + ds/ds0);
            if (e_rw[jnode] <= 0.) e_rw[jnode] = ds/ds0;
            else e_rw[jnode] = 0.5*(e_rw[jnode] + ds/ds0);
          }
        }

        vector2hdf5(h5f, groupname + "/points", x_rw, Nrw);
        if (!prm.minimal_output)
          vector2hdf5(h5f, groupname + "/u", u_rw, Nrw);
        if (prm.output_all_props && !prm.minimal_output){
          scalar2hdf5(h5f, groupname + "/rho", rho_rw, Nrw);
          scalar2hdf5(h5f, groupname + "/p", p_rw, Nrw);
        }
        if (!prm.minimal_output)
          scalar2hdf5(h5f, groupname + "/c", c_rw, Nrw);
        if (!prm.minimal_output && edges.size() > 0)
          scalar2hdf5(h5f, groupname + "/H", H_rw, Nrw);
        if (!prm.minimal_output && faces.size() > 0)
          vector2hdf5(h5f, groupname + "/n", n_rw, Nrw);
        if (faces.size() == 0)
          scalar2hdf5(h5f, groupname + "/e", e_rw, Nrw);
        if (edges.size() > 0)
          mesh2hdf(h5f, groupname, x_rw, faces, edges);
      }
    }
    for (Uint irw=0; irw < Nrw; ++irw){
      dx_rw = u_rw[irw]*dt;

      // Set elongation
      e_rw[irw] = 0.;

      // Second-order terms
      if (prm.int_order >= 2){
        dx_rw += 0.5*a_rw[irw]*dt2;
      }
      if (Dm > 0.0){
        Vector3d eta = {rnd_normal(gen),
                        rnd_normal(gen),
                        rnd_normal(gen)};
        dx_rw += sqrt2Dmdt*eta;
      }
      intp->probe(x_rw[irw]+dx_rw);
     if (intp->inside_domain()){
        x_rw[irw] += dx_rw;
        u_rw[irw] = U0*intp->get_u();

        if ((it+1) % int_dump_intv == 0){
          rho_rw[irw] = intp->get_rho();
          p_rw[irw] = intp->get_p();
        }

        // Second-order terms
        if (prm.int_order >= 2){
          a_rw[irw] = U02*intp->get_Ju() + U0*intp->get_a();
        }
        n_accepted++;
      }
      else {
        n_declined++;
        declinedfile << t << " "
                     << x_rw[irw][0]+dx_rw[0] << " "
                     << x_rw[irw][1]+dx_rw[1] << " "
                     << x_rw[irw][2]+dx_rw[2] << std::endl;
      }
    }
    t += dt;
    it += 1;
  }

 // Final checkpoint
  prm.t = t;
  prm.n_accepted = n_accepted;
  prm.n_declined = n_declined;
  prm.Nrw = Nrw;
  prm.dump(checkpointsfolder);
  dump_positions(checkpointsfolder + "/positions.pos", x_rw, Nrw);
  dump_faces(checkpointsfolder + "/faces.face", faces);
  dump_edges(checkpointsfolder + "/edges.edge", edges);
  dump_colors(checkpointsfolder + "/colors.col", c_rw, Nrw);
  if (prm.inject){
    dump_positions(checkpointsfolder + "/positions_inj.pos", pos_inj);
    dump_edges(checkpointsfolder + "/edges_inj.edge", edges_inj);
    dump_list(checkpointsfolder + "/edges_inlet.list", edges_inlet);
    dump_list(checkpointsfolder + "/nodes_inlet.list", nodes_inlet);
  }

  // Close files
  statfile.close();
  declinedfile.close();
  if (write_mode == "hdf5"){
    h5f->close();
  }

  delete intp;
  return 0;
}
