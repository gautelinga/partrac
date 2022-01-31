#ifndef __CHECKPOINTS_HPP
#define __CHECKPOINTS_HPP

#include "H5Cpp.h"
#include "typedefs.hpp"

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

std::string get_newfoldername(const std::string rwfolder, const Parameters& prm){
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