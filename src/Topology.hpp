#ifndef __TOPOLOGY_HPP
#define __TOPOLOGY_HPP

#include "typedefs.hpp"
#include "mesh.hpp"
#include "Initializer.hpp"
#include "stats.hpp"

class Topology {
public:
  Topology(ParticleSet& ps, const Parameters &prm);
  int dim();
  void compute_maps();
  void clear();
  Uint refine();
  Uint coarsen();
  bool filter();
  Uint inject();
  void compute_interior();
  EdgesType edges;
  FacesType faces;
  Edge2FacesType edge2faces;
  Node2EdgesType node2edges;
  std::vector<Vector3d> pos_inj;
  EdgesType edges_inj;
  EdgesListType edges_inlet;
  NodesListType nodes_inlet;
  // Initial curvature?
  InteriorAnglesType interior_ang;
  std::vector<double> mixed_areas;
  std::vector<Vector3d> face_normals;
  void write_checkpoint(const std::string checkpointsfolder, const double t, Parameters &prm) const;
  void load_checkpoint(const std::string checkpointsfolder, const Parameters &prm);
  void dump_hdf5(H5FilePtr h5f, const std::string groupname);
  void load_initial_state(Initializer* init_state);
  void write_statistics(std::ofstream &statfile, const double t, const double ds_max,
                        const bool do_dump_hist, const std::string histfolder, Integrator* integrator);
private:
  ParticleSet& ps;
  double ds_min;
  double ds_max;
  double curv_refine_factor;
  bool cut_if_stuck;
  bool inject_edges;
  bool verbose;
  Uint filter_target;
};

Topology::Topology(ParticleSet& ps, const Parameters &prm) : ps(ps) {
  ds_min = prm.ds_min;
  ds_max = prm.ds_max;
  curv_refine_factor = prm.curv_refine_factor;
  cut_if_stuck = prm.cut_if_stuck;
  inject_edges = prm.inject_edges;
  verbose = prm.verbose;
  filter_target = prm.filter_target;
}

int Topology::dim(){
    if (faces.size() > 0)
        return 2;
    else if (edges.size() > 0)
        return 1;
    return 0;
}

void Topology::compute_maps() {
  // Compute edge2faces map
  compute_edge2faces(edge2faces, faces, edges);
  // Compute node2edges map
  compute_node2edges(node2edges, edges, ps.N());
}

void Topology::clear(){
  edges.clear();
  faces.clear();
}

Uint Topology::refine(){
  return refinement(faces, edges,
                    edge2faces, node2edges,
                    edges_inlet,
                    ps, ds_max,
                    curv_refine_factor,
                    cut_if_stuck);
}

Uint Topology::coarsen(){
  return coarsening(faces, edges,
                    edge2faces, node2edges,
                    edges_inlet, nodes_inlet,
                    ps, ds_min,
                    curv_refine_factor);
}

Uint Topology::inject(){
  return injection(pos_inj,
                   edges_inj,
                   edges_inlet,
                   nodes_inlet,
                   edges, 
                   faces,
                   edge2faces,
                   node2edges,
                   ps,
                   inject_edges,
                   verbose);
}

void Topology::compute_interior(){
  /*
  remove_unused_edges(faces, edges, edges_inlet);
  remove_unused_nodes(edges, nodes_inlet, ps);

  compute_edge2faces(edge2faces, faces, edges);
  compute_node2edges(node2edges, edges, ps.N());
  */

  compute_interior_prop(interior_ang, mixed_areas, face_normals,
                        faces, edges, edge2faces, ps);
  compute_mean_curv(faces, edges, edge2faces, node2edges,
                    ps, interior_ang, mixed_areas, face_normals);
  /**/
}

bool Topology::filter(){
  return filtering(faces, edges, edge2faces, node2edges, ps, filter_target);
}

void Topology::write_checkpoint(const std::string checkpointsfolder, const double t, Parameters &prm) const {
  prm.t = t;
  prm.n_accepted = ps.get_accepted();
  prm.n_declined = ps.get_declined();
  //prm.Nrw = Nrw;
  prm.dump(checkpointsfolder);
  // dump_positions(checkpointsfolder + "/positions.pos", ps.x_rw, ps.Nrw);
  ps.dump_positions(checkpointsfolder + "/positions.pos");
  dump_faces(checkpointsfolder + "/faces.face", faces);
  dump_edges(checkpointsfolder + "/edges.edge", edges);
  //dump_colors(checkpointsfolder + "/colors.col", ps.c_rw, ps.Nrw);
  ps.dump_scalar(checkpointsfolder + "/colors.col", "c");
  if (prm.inject){
    dump_vector_field(checkpointsfolder + "/positions_inj.pos", pos_inj);
    dump_edges(checkpointsfolder + "/edges_inj.edge", edges_inj);
    dump_list(checkpointsfolder + "/edges_inlet.list", edges_inlet);
    dump_list(checkpointsfolder + "/nodes_inlet.list", nodes_inlet);
  }
  if (prm.local_dt){
    //dump_colors(checkpointsfolder + "/tau.dat", ps.tau_rw, ps.Nrw);
    ps.dump_scalar(checkpointsfolder + "/tau.dat", "tau");
  }
}

void Topology::load_checkpoint(const std::string checkpointsfolder, const Parameters &prm){
  std::string posfile = checkpointsfolder + "/positions.pos";
  //load_positions(posfile, pos_init, prm.Nrw);
  ps.load_positions(posfile);
  std::string facefile = checkpointsfolder + "/faces.face";
  load_faces(facefile, faces);
  std::string edgefile = checkpointsfolder + "/edges.edge";
  load_edges(edgefile, edges);
  std::string colfile = checkpointsfolder + "/colors.col";
  //load_colors(colfile, ps.c_rw, prm.Nrw);
  ps.load_scalar(colfile, "c");
  if (prm.inject){
    std::string posinjfile = checkpointsfolder + "/positions_inj.pos";
    load_vector_field(posinjfile, pos_inj);
    std::string edgesinjfile = checkpointsfolder + "/edges_inj.edge";
    load_edges(edgesinjfile, edges_inj);
    std::string edgesinletfile = checkpointsfolder + "/edges_inlet.list";
    std::string nodesinletfile = checkpointsfolder + "/nodes_inlet.list";
    load_list(edgesinletfile, edges_inlet);
    load_list(nodesinletfile, nodes_inlet);
  }
  if (prm.local_dt){
    std::string taufile = checkpointsfolder + "/tau.dat";
    //load_colors(taufile, ps.tau_rw, prm.Nrw);
    ps.load_scalar(taufile, "tau");
  }
}

void Topology::dump_hdf5(H5FilePtr h5f, const std::string groupname){
    if (dim() > 0)
        mesh2hdf(h5f, groupname, ps, faces, edges);
}

void Topology::load_initial_state(Initializer* init_state){
  /*std::vector<Vector3d> pos_init;
  pos_init = initial_positions(prm.init_mode,
                               prm.init_weight,
                               prm.Nrw,
                               {prm.x0, prm.y0, prm.z0},
                               prm.La, prm.Lb,
                               prm.ds_init, t0,
                                 intp, gen, mesh.edges, mesh.faces);*/
  clear();

  edges = init_state->edges;
  faces = init_state->faces;

  if (init_state->clear_initial_edges){
    std::cout << "Clearing initial edges!" << std::endl;
    clear();
  }
  if (init_state->inject){
    pos_inj = init_state->nodes;
    edges_inj = init_state->edges;
    for (Uint i=0; i<init_state->nodes.size(); ++i){
      nodes_inlet.push_back(i);
    }
    for (Uint i=0; i<init_state->edges.size(); ++i){
      edges_inlet.push_back(i);
    }
  }
  ps.add(init_state->nodes, 0);
}

void Topology::write_statistics(std::ofstream &statfile,
                           const double t,
                           const double ds_max,
                           const bool do_dump_hist,
                           const std::string histfolder,
                           Integrator* integrator){
  write_stats(statfile, t, ps, faces, edges, ds_max, do_dump_hist, histfolder,
              integrator->get_accepted(), integrator->get_declined());
}

#endif