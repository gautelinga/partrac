#ifndef __HELPERS_HPP
#define __HELPERS_HPP

#include <random>
#include "H5Cpp.h"

#include "typedefs.hpp"

#include "Interpol.hpp"
#include "StructuredInterpol.hpp"
#include "StructuredConstInterpol.hpp"
#include "AnalyticInterpol.hpp"
#ifdef USE_DOLFIN
#include "DolfInterpol.hpp"
#include "TetInterpol.hpp"
#include "TriangleInterpol.hpp"
#include "TriangleFreqInterpol.hpp"
#include "XDMFTriangleInterpol.hpp"
#endif
#include "Initializer.hpp"


void set_interpolate_mode(std::shared_ptr<Interpol>& intp, const std::string& mode, const std::string& infilename){
  if (mode == "analytic"){
    std::cout << "AnalyticInterpol initiated." << std::endl;
    intp = std::make_shared<AnalyticInterpol>(infilename);
  }
  else if (mode == "unstructured" || mode == "fenics" || mode == "xdmf" ||
           mode == "tet" || mode == "triangle" || mode == "trianglefreq" || mode == "xdmftriangle"){
#ifdef USE_DOLFIN
    if (mode == "tet"){
      intp = std::make_shared<TetInterpol>(infilename);
    }
    else if (mode == "triangle"){
      intp = std::make_shared<TriangleInterpol>(infilename);
    }
    else if (mode == "trianglefreq"){
      intp = std::make_shared<TriangleFreqInterpol>(infilename);
    }
    else if (mode == "xdmftriangle"){
      intp = std::make_shared<XDMFTriangleInterpol>(infilename);
    }
    else if (mode == "fenics"){
      intp = std::make_shared<DolfInterpol>(infilename);
    }
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
    intp = std::make_shared<StructuredInterpol>(infilename);
  }
  else {
    std::cout << "Mode not supported." << std::endl;
    exit(0);
  }
}

void set_initial_state(std::shared_ptr<Initializer>& init_state, std::shared_ptr<Interpol> intp, MPIwrap& mpi, Parameters& prm, std::mt19937& gen){
  std::vector<std::string> key = split_string(prm.init_mode, "_");
  if (key.size() == 0){
    std::cout << "init_mode not specified." << std::endl;
    exit(0);
  }
  else if (key[0] == "point"){
    //init_state = new PointInitializer(key, intp, prm, mpi);
    init_state = std::make_shared<PointInitializer>(key, intp, prm, mpi);
  }
  else if (key[0] == "uniform"){
    init_state = std::make_shared<UniformInitializer>(key, intp, prm, mpi);
  }
  else if (key[0] == "strip"){
    init_state = std::make_shared<StripInitializer>(key, intp, prm, mpi);
  }
  else if (key[0] == "sheet"){
    init_state = std::make_shared<SheetInitializer>(key, intp, prm, mpi);
  }
  else if (key[0] == "ellipsoid"){
    init_state = std::make_shared<EllipsoidInitializer>(key, intp, prm, mpi);
  }
  else if (key[0] == "pair" || key[0] == "pairs"){
    init_state = std::make_shared<RandomPairsInitializer>(key, intp, prm, mpi, gen);
  }
  else if (key[0] == "points"){
    init_state = std::make_shared<RandomPointsInitializer>(key, intp, prm, mpi, gen);
  }
  else if (key[0] == "randomgaussianstrip"){
    init_state = std::make_shared<RandomGaussianStripInitializer>(key, intp, prm, mpi, gen);
  }
  else if (key[0] == "randomgaussiancircle"){
    init_state = std::make_shared<RandomGaussianCircleInitializer>(key, intp, prm, mpi, gen);
  }
  else {
    std::cout << "Unknown init_mode: " << prm.init_mode << std::endl;
    exit(0);
  }
}

static void test_interpolation(Uint num_points, std::shared_ptr<Interpol> intp,
                        const std::string &newfolder, const double t0,
                        std::mt19937 &gen){
  Uint n = 0;
  Uint n_inside = 0;

  Vector3d x_min = intp->get_x_min();
  Vector3d x_max = intp->get_x_max();

  double Lx = intp->get_Lx();
  double Ly = intp->get_Ly();
  double Lz = intp->get_Lz();

  std::cout << "Lx=" << Lx << ", Ly=" << Ly << ", Lz=" << Lz << std::endl;

  intp->update(t0);

  // std::ofstream nodalfile(newfolder + "/nodal_values.dat");
  // for (Uint ix=0; ix<intp->get_nx(); ++ix){
  //   for (Uint iy=0; iy<intp->get_ny(); ++iy){
  //     for (Uint iz=0; iz<intp->get_nz(); ++iz){
  //       bool inside = intp->get_nodal_inside(ix, iy, iz);
  //       Vector3d u(intp->get_nodal_ux(ix, iy, iz),
  //                  intp->get_nodal_uy(ix, iy, iz),
  //                  intp->get_nodal_uz(ix, iy, iz));
  //       nodalfile << ix << " " << iy << " " << iz << " " << inside << " "
  //                 << u[0] << " " << u[1] << " " << u[2] << std::endl;
  //     }
  //   }
  // }
  // nodalfile.close();

  std::uniform_real_distribution<> uni_dist_x(x_min[0], x_max[0]);
  std::uniform_real_distribution<> uni_dist_y(x_min[1], x_max[1]);
  std::uniform_real_distribution<> uni_dist_z(x_min[2], x_max[2]);

  std::ofstream ofile(newfolder + "/interpolation.txt");

  std::string sep = ",";
  ofile << "x" << sep << "y" << sep << "z" << sep
        << "ux" << sep << "uy" << sep << "uz" << sep
        << "rho" << sep << "p" << sep << "divu" << sep
        << "vortz" << sep
        << "uxx" << sep << "uxy" << sep << "uxz" << sep
        << "uyx" << sep << "uyy" << sep << "uyz" << sep
        << "uzx" << sep << "uzy" << sep << "uzz"
        << std::endl;

  while (n < num_points){
    Vector3d x(uni_dist_x(gen), uni_dist_y(gen), uni_dist_z(gen));
    intp->probe(x, t0);
    if (intp->inside_domain()){
      Vector3d u = intp->get_u();
      double rho = intp->get_rho();
      double p = intp->get_p();
      double divu = intp->get_divu();
      double vortz = intp->get_vortz();
      ofile << x[0] << sep << x[1] << sep << x[2] << sep
            << u[0] << sep << u[1] << sep << u[2] << sep
            << rho << sep << p << sep << divu << sep
            << vortz << sep
            << intp->get_uxx() << sep << intp->get_uxy() << sep << intp->get_uxz() << sep
            << intp->get_uyx() << sep << intp->get_uyy() << sep << intp->get_uyz() << sep
            << intp->get_uzx() << sep << intp->get_uzy() << sep << intp->get_uzz()
            << std::endl;
      ++n_inside;
    }
    ++n;
  }
  std::cout << "Inside: " << n_inside << "/" << n << std::endl;
  std::cout << "Approximate volume: " << (n_inside*Lx*Ly*Lz)/n << std::endl;
  std::cout << "Approximate area:   " << (n_inside*Lx*Ly)/n << std::endl;
  ofile.close();
}

/*std::vector<double> areas(const std::vector<Uint> &kfaces,
                          std::vector<Vector3d>& x_rw,
                          const FacesType &faces, const EdgesType &edges){
  // To be decommisioned??
  
  std::vector<double> a;
  for (std::vector<Uint>::const_iterator faceit=kfaces.begin();
       faceit != kfaces.end(); ++faceit){
    a.push_back(area(*faceit, x_rw, faces, edges));
  }
  return a;
}*/

/*
void write_checkpoint(const std::string checkpointsfolder, const double t, Parameters &prm, const ParticleSet& ps, 
                      const Topology& mesh,
                      const Integrator* integrator){
  prm.t = t;
  prm.n_accepted = integrator->get_accepted();
  prm.n_declined = integrator->get_declined();
  //prm.Nrw = Nrw;
  prm.dump(checkpointsfolder);
  // dump_positions(checkpointsfolder + "/positions.pos", ps.x_rw, ps.Nrw);
  ps.dump_positions(checkpointsfolder + "/positions.pos");
  dump_faces(checkpointsfolder + "/faces.face", mesh.faces);
  dump_edges(checkpointsfolder + "/edges.edge", mesh.edges);
  //dump_colors(checkpointsfolder + "/colors.col", ps.c_rw, ps.Nrw);
  ps.dump_scalar(checkpointsfolder + "/colors.col", "c");
  if (prm.inject){
    dump_vector_field(checkpointsfolder + "/positions_inj.pos", mesh.pos_inj);
    dump_edges(checkpointsfolder + "/edges_inj.edge", mesh.edges_inj);
    dump_list(checkpointsfolder + "/edges_inlet.list", mesh.edges_inlet);
    dump_list(checkpointsfolder + "/nodes_inlet.list", mesh.nodes_inlet);
  }
  if (prm.local_dt){
    //dump_colors(checkpointsfolder + "/tau.dat", ps.tau_rw, ps.Nrw);
    ps.dump_scalar(checkpointsfolder + "/tau.dat", "tau");
  }
}

void load_checkpoint(const Parameters &prm, ParticleSet& ps, 
                     Topology& mesh){
  std::string posfile = prm.restart_folder + "/Checkpoints/positions.pos";
  //load_positions(posfile, pos_init, prm.Nrw);
  ps.load_positions(posfile);
  std::string facefile = prm.restart_folder + "/Checkpoints/faces.face";
  load_faces(facefile, mesh.faces);
  std::string edgefile = prm.restart_folder + "/Checkpoints/edges.edge";
  load_edges(edgefile, mesh.edges);
  std::string colfile = prm.restart_folder + "/Checkpoints/colors.col";
  //load_colors(colfile, ps.c_rw, prm.Nrw);
  ps.load_scalar(colfile, "c");
  if (prm.inject){
    std::string posinjfile = prm.restart_folder + "/Checkpoints/positions_inj.pos";
    load_vector_field(posinjfile, mesh.pos_inj);
    std::string edgesinjfile = prm.restart_folder + "/Checkpoints/edges_inj.edge";
    load_edges(edgesinjfile, mesh.edges_inj);
    std::string edgesinletfile = prm.restart_folder + "/Checkpoints/edges_inlet.list";
    std::string nodesinletfile = prm.restart_folder + "/Checkpoints/nodes_inlet.list";
    load_list(edgesinletfile, mesh.edges_inlet);
    load_list(nodesinletfile, mesh.nodes_inlet);
  }
  if (prm.local_dt){
    std::string taufile = prm.restart_folder + "/Checkpoints/tau.dat";
    //load_colors(taufile, ps.tau_rw, prm.Nrw);
    ps.load_scalar(taufile, "tau");
  }
}*/

std::string get_newfoldername(const std::string rwfolder, Parameters& prm){
  std::ostringstream ss_Dm, ss_dt, ss_Nrw, ss_seed;
  ss_Dm << std::scientific << std::setprecision(7) << prm.Dm;
  ss_dt << std::scientific << std::setprecision(7) << prm.dt;
  ss_Nrw << prm.Nrw;
  ss_seed << prm.seed;
  std::string newfoldername = rwfolder +
                            "/Dm" + ss_Dm.str() + // "_U" + std::to_string(prm.U0) +
                            "_dt" + ss_dt.str() +
                            "_Nrw" + ss_Nrw.str() +
                            "_seed" + ss_seed.str() +
                            prm.tag +
                            "/";
  return newfoldername;
}

#endif