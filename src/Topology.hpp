#ifndef __TOPOLOGY_HPP
#define __TOPOLOGY_HPP

#include "typedefs.hpp"
#include "mesh.hpp"
#include "Interpol.hpp"
#include "Initializer.hpp"
#include "stats.hpp"
#include "MPIwrap.hpp"


class Topology {
public:
  Topology(ParticleSet& ps, const Parameters &prm, MPIwrap& mpi);
  int dim();
  void compute_maps();
  void clear();
  Uint refine();
  Uint coarsen();
  bool filter();
  Uint remove_beyond(const int, const double);
  Uint inject();
  void compute_interior();
  void remove_nodes_safe(std::vector<bool>&);
  bool resize(const double);
  void integrate_tau(const double, const double);
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
  void write_checkpoint(const std::string& checkpointsfolder, const double t, Parameters &prm) const;
  void load_checkpoint(const std::string& checkpointsfolder, const Parameters &prm);
  void dump_hdf5(H5::H5File& h5f, const std::string& groupname, std::map<std::string, bool>& output_fields);
  void load_initial_state(std::shared_ptr<Initializer> init_state);
  template<typename T>
  void write_statistics(std::ofstream &statfile, const double t, const double ds_max, //const bool do_dump_hist, const std::string histfolder, 
                        T& integrator);
                        //std::shared_ptr<Integrator> integrator);
private:
  ParticleSet& ps;
  double ds_min;
  double ds_max;
  double curv_refine_factor;
  bool cut_if_stuck;
  bool inject_edges;
  bool verbose;
  Uint filter_target;

  MPIwrap& m_mpi;
};

Topology::Topology(ParticleSet& ps, const Parameters &prm, MPIwrap& mpi) : ps(ps), m_mpi(mpi) {
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

void Topology::integrate_tau(const double dt, const double tau_max){
  if (faces.size() == 0){
    for ( auto & edge : edges ){
      Uint inode = edge.first[0];
      Uint jnode = edge.first[1];
      double ds0 = edge.second;
      double ds = ps.dist(inode, jnode);
      double rho = ds / ds0;
      edge.tau += 0.5*dt*(pow(edge.rho_prev, 2) + pow(rho, 2));
      edge.rho_prev = rho;
    }
    if (tau_max > 0.0){
      std::vector<bool> edge_isactive(edges.size(), true); 
      for (Uint iedge=0; iedge < edges.size(); ++iedge){
        if (edges[iedge].tau > tau_max)
          edge_isactive[iedge] = false;
      }

      // untested!
      FacesType faces_dummy;
      //NodesListType nodes_inlet_dummy;
      EdgesListType edges_inlet_dummy;
      remove_edges(faces_dummy, edges, edge_isactive, edges_inlet_dummy);
      remove_unused_nodes(edges, nodes_inlet, ps);
      compute_node2edges(node2edges, edges, ps.N());
    }
  }
  else {
    for ( auto & face : faces )
    {
      Uint iedge = face.first[0];
      Uint jedge = face.first[1];
      double dA = ps.triangle_area(iedge, jedge, edges);
      double dA0 = face.second;
      double rho = dA/dA0;
      face.tau += 0.5*dt*(pow(face.rho_prev, 2) + pow(rho, 2));
      face.rho_prev = rho;
    }
  }
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

void Topology::remove_nodes_safe(std::vector<bool>& node_isactive){
  remove_nodes_safely(faces, edges,
                      edge2faces, node2edges,
                      edges_inlet, nodes_inlet,
                      node_isactive,
                      ps);
}

bool Topology::filter(){
  return filtering(faces, edges, edge2faces, node2edges, ps, filter_target);
}

Uint Topology::remove_beyond(const int exit_dim, const double Ln){
  std::vector<bool> node_isactive(ps.N(), true);
  Uint count = 0;
  for (Uint i=0; i<ps.N(); ++i){
    Vector3d xi = ps.x(i);
    if (xi[exit_dim] > Ln){
      node_isactive[i] = false;
      ++count;
    }
  }
  remove_nodes_safe(node_isactive);
  return count;
}

bool Topology::resize(const double ds){
  return resizing(edges, node2edges, ps, ds);
}

void Topology::write_checkpoint(const std::string& checkpointsfolder, const double t, Parameters &prm) const {
  prm.t = t;
  //prm.n_accepted = ps.get_accepted();
  //prm.n_declined = ps.get_declined();
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

void Topology::load_checkpoint(const std::string& checkpointsfolder, const Parameters &prm){
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

void Topology::dump_hdf5(H5::H5File& h5f, const std::string& groupname, std::map<std::string, bool>& output_fields){
  /*
  std::cout << "Process " << m_mpi.rank() << ": " << ps.N() << " " << faces.size() << " " << edges.size() << std::endl;
  auto num_nodes_ = m_mpi.gather(ps.N());
  auto num_edges_ = m_mpi.gather(edges.size());
  std::vector<int> cum_nodes_; 
  std::vector<int> cum_edges_;
  cum_nodes_.resize(m_mpi.size());
  cum_edges_.resize(m_mpi.size());
  if (m_mpi.rank() == 0){
    for (int i=0; i<m_mpi.size(); ++i){
      cum_nodes_[i] = (i > 0) ? (cum_nodes_[i-1] + num_nodes_[i-1]) : 0;
      cum_edges_[i] = (i > 0) ? (cum_edges_[i-1] + num_edges_[i-1]) : 0;
    }
  }
  auto inode0 = m_mpi.scatter(cum_nodes_);
  auto iedge0 = m_mpi.scatter(cum_edges_);
  // std::cout << "p" << m_mpi.rank() << " " << inode0 << std::endl;

  EdgesType all_edges = m_mpi.wrapEdges(edges, inode0);
  FacesType all_faces = m_mpi.wrapFaces(faces, iedge0);
  // std::cout << "WORKED: " << m_mpi.rank() << " " << all_edges.size() << std::endl;
  
  for (auto &el : all_edges){
    std::cout << "(" << el.first[0] << ", " << el.first[1] << "): " << el.second << std::endl;
  }
  for (auto &el : all_faces){
    std::cout << "(" << el.first[0] << ", " << el.first[1] << ", " << el.first[2] << "): " << el.second << std::endl;
  }

  m_mpi.barrier();

  Uint all_Nrw = cum_nodes_[m_mpi.size()-1] + num_nodes_[m_mpi.size()-1];
  Interpol* intp_dummy;
  ParticleSet all_ps(intp_dummy, all_Nrw, m_mpi);
  all_ps.reduce(ps, output_fields);

  if (dim() > 0)
    mesh2hdf(h5f, groupname, all_ps, all_faces, all_edges);
  all_ps.dump_hdf5(h5f, groupname, output_fields);*/
  
  if (dim() > 0)
    mesh2hdf(h5f, groupname, ps, faces, edges, output_fields["tau"]);
  ps.dump_hdf5(h5f, groupname, output_fields);
}

void Topology::load_initial_state(std::shared_ptr<Initializer> init_state){
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

template<typename T>
void Topology::write_statistics( std::ofstream &statfile
                               , const double t
                               , const double ds_max
                               , T& integrator){
                           //const bool do_dump_hist,
                           //const std::string histfolder,
                           //std::shared_ptr<Integrator> integrator){
  write_stats(m_mpi, statfile, t, ps, faces, edges, ds_max, //do_dump_hist, histfolder,
              integrator.get_accepted(), integrator.get_declined());
              //integrator->get_accepted(), integrator->get_declined());
}

#endif