#ifndef __TOPOLOGY_HPP
#define __TOPOLOGY_HPP

#include "typedefs.hpp"
#include "mesh.hpp"
#include "Interpol.hpp"
#include "Initializer.hpp"
#include "stats.hpp"


class Topology {
public:
  Topology(ParticleSet& ps, const partrac::Params& prm);
  int dim();
  int dim_settled();
  void check_topology();
  void check_dim();
  void compute_maps();
  void clear();
  Uint refine();
  Uint coarsen(const bool full);
  bool filter();
  Uint remove_beyond(const int, const double);
  Uint inject();
  void compute_interior();
  // Whether compute_interior fills H and n
  bool computes_curvature() const { return curv_refine_factor > 0.; };
  void remove_nodes_safe(std::vector<bool>&);
  bool resize(const double);
  void integrate_tau(const double, const double);
  EdgesType edges;
  FacesType faces;
  Edge2FacesType edge2faces;
  Node2EdgesType node2edges;
  std::vector<Vector3d> pos_inj;
  int dim0 = -1;   // the dimension the mesh settled into, -1 until it has one
  double Dm = 0.;  // not every app with a mesh has an integrator to declare it
  EdgesType edges_inj;
  EdgesListType edges_inlet;
  NodesListType nodes_inlet;
  // Initial curvature?
  InteriorAnglesType interior_ang;
  std::vector<double> mixed_areas;
  std::vector<Vector3d> face_normals;
  void write_checkpoint(const std::string& checkpointsfolder, const double t, partrac::Params& prm) const;
  void load_checkpoint(const std::string& checkpointsfolder, const partrac::Params& prm);
  void dump_hdf5(H5::H5File& h5f, const std::string& groupname, std::map<std::string, bool>& output_fields);
  void load_initial_state(std::shared_ptr<Initializer> init_state, partrac::Params& prm);
  std::vector<StatsColumn> stats_header_columns(const double ds_max){
    return mesh_stats_columns(0., ps, faces, edges, ds_max, 0, 0, dim_settled());
  }
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
  bool injecting;
  bool inject_edges;
  bool verbose;
  Uint filter_target;
};

inline Topology::Topology(ParticleSet& ps, const partrac::Params& prm) : ps(ps) {
  ds_min = prm.get<double>("ds_min");
  ds_max = prm.get<double>("ds_max");
  curv_refine_factor = prm.get<double>("curv_refine_factor");
  cut_if_stuck = prm.get<bool>("cut_if_stuck");
  injecting = prm.has("inject") ? prm.get<bool>("inject") : false;
  inject_edges = prm.get<bool>("inject_edges");
  verbose = prm.get<bool>("verbose");
  filter_target = prm.get<int>("filter_target");
  Dm = prm.has("Dm") ? prm.get<double>("Dm") : 0.;
}

// Nothing left to advect, or a manifold that changed dimension: a strip that
// has lost its last edge, or a sheet its last face, has nothing left to
// stretch, and every row written afterwards would sit under a header that no
// longer describes it. Injection can legitimately start from nothing, so the
// dimension latches on first use.
// An empty set is the caller's policy; a changed dimension is always wrong
inline void Topology::check_topology(){
  if (ps.N() == 0){
    std::cerr << "Error: no particles left. Stopping." << std::endl;
    exit(1);
  }
  check_dim();
}

inline void Topology::check_dim(){
  if (ps.N() == 0) return;   // no nodes, so no manifold to have changed
  const int d = dim();
  if (dim0 < 0){
    if (d > 0){
      dim0 = d;
      if (Dm > 0.)
        std::cerr << "Warning: Dm = " << Dm << " > 0 on a " << d
                  << "-dimensional mesh. Brownian motion and material "
                  << "deformation are normally not compatible." << std::endl;
    }
    return;
  }
  if (d == dim0) return;
  // An injecting run traces a dimension more than its inlet, so it may rise
  if (d > dim0 && injecting){
    dim0 = d;
    return;
  }
  std::cerr << "Error: the mesh changed dimension, from " << dim0 << " to " << d
            << ", with " << ps.N() << " nodes and " << edges.size()
            << " edges left. Stopping." << std::endl;
  exit(1);
}

// The dimension the run settles into. Injection raises what the inlet
// traces by one and the header is written before the first, so the larger
// of the two is what to ask for -- on a restart as well
inline int Topology::dim_settled(){
  int d = dim();
  if (injecting && inject_edges){
    const int inlet = edges_inj.size() > 0 ? 1
                    : (pos_inj.size() > 0 ? 0 : -1);
    if (inlet >= 0)
      d = std::max(d, inlet + 1);
  }
  return d;
}

inline int Topology::dim(){
    if (faces.size() > 0)
        return 2;
    else if (edges.size() > 0)
        return 1;
    return 0;
}

inline void Topology::compute_maps() {
  // Compute edge2faces map
  compute_edge2faces(edge2faces, faces, edges);
  // Compute node2edges map
  compute_node2edges(node2edges, edges, ps.N());
}

inline void Topology::clear(){
  edges.clear();
  faces.clear();
}

inline void Topology::integrate_tau(const double dt, const double tau_max){
  if (faces.size() == 0){
    // Each edge advances its own tau from its own rho_prev and read-only
    // geometry, so the order is immaterial and the result is unchanged
    #pragma omp parallel for
    for (Uint i = 0; i < edges.size(); ++i){
      auto & edge = edges[i];
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

      std::vector<bool> face_isactive(faces.size(), true);   // a strip has none
      std::vector<bool> node_isactive(ps.N(), true);
      remove_inactive(faces, edges, edge2faces, node2edges,
                      edges_inlet, nodes_inlet,
                      face_isactive, edge_isactive, node_isactive, ps);
      check_topology();   // culling every edge would drop the dimension
    }
  }
  else {
    #pragma omp parallel for
    for (Uint i = 0; i < faces.size(); ++i)
    {
      auto & face = faces[i];
      Uint iedge = face.first[0];
      Uint jedge = face.first[1];
      double dA0 = face.second;
      if (!(dA0 > 0.))
        continue;                 // a flat sweep, waiting to be culled
      double dA = ps.triangle_area(iedge, jedge, edges);
      double rho = dA/dA0;
      face.tau += 0.5*dt*(pow(face.rho_prev, 2) + pow(rho, 2));
      face.rho_prev = rho;
    }
    if (tau_max > 0.0){
      std::vector<bool> face_isactive(faces.size(), true);
      for (Uint iface=0; iface < faces.size(); ++iface){
        if (faces[iface].tau > tau_max)
          face_isactive[iface] = false;
      }
      std::vector<bool> edge_isactive(edges.size(), true);
      std::vector<bool> node_isactive(ps.N(), true);
      remove_inactive(faces, edges, edge2faces, node2edges,
                      edges_inlet, nodes_inlet,
                      face_isactive, edge_isactive, node_isactive, ps);
      check_topology();
    }
  }
}

inline Uint Topology::refine(){
  Uint n = refinement(faces, edges,
                      edge2faces, node2edges,
                      edges_inlet, nodes_inlet,
                      pos_inj, edges_inj,
                      ps, ds_max,
                      curv_refine_factor,
                      cut_if_stuck);
  check_topology();
  return n;
}

// Coarsening runs whether or not it was asked for. With it off the threshold
// drops to what is numerically zero: refinement raises a median from the new
// node to the opposite vertex of each face it splits, and that edge's length
// is the face's business, not ds_max's -- on a nearly collinear face it comes
// out at nothing, and the faces on it then have two vertices at one point.
inline Uint Topology::coarsen(const bool full){
  double ds_cut = ds_min;
  if (!full){
    // Scaled from the mesh rather than from ds_max, which a run that does not
    // refine is free to leave enormous
    double ds_longest = 0.;
    for ( const auto & edge : edges )
      ds_longest = std::max(ds_longest,
                            ps.dist(edge.first[0], edge.first[1]));
    ds_cut = 1e-12 * ds_longest;
  }
  Uint n = coarsening(faces, edges,
                      edge2faces, node2edges,
                      edges_inlet, nodes_inlet,
                      ps, ds_cut,
                      curv_refine_factor);
  check_topology();
  return n;
}

inline Uint Topology::inject(){
  Uint n = injection(pos_inj,
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
  // the generation before this one is no longer at the inlet
  std::vector<bool> face_isactive(faces.size(), true);
  std::vector<bool> edge_isactive(edges.size(), true);
  std::vector<bool> node_isactive(ps.N(), true);
  remove_inactive(faces, edges, edge2faces, node2edges,
                  edges_inlet, nodes_inlet,
                  face_isactive, edge_isactive, node_isactive, ps);
  check_topology();
  return n;
}

inline void Topology::compute_interior(){
  // Only the curvature-weighted refinement reads this
  if (!computes_curvature())
    return;
  compute_interior_prop(interior_ang, mixed_areas, face_normals,
                        faces, edges, edge2faces, ps);
  compute_mean_curv(faces, edges, edge2faces, node2edges,
                    ps, interior_ang, mixed_areas, face_normals);
}

inline void Topology::remove_nodes_safe(std::vector<bool>& node_isactive){
  std::vector<bool> face_isactive(faces.size(), true);
  std::vector<bool> edge_isactive(edges.size(), true);
  remove_inactive(faces, edges,
                  edge2faces, node2edges,
                  edges_inlet, nodes_inlet,
                  face_isactive, edge_isactive, node_isactive,
                  ps);
  check_dim();
}

inline bool Topology::filter(){
  bool changed = filtering(faces, edges, edge2faces, node2edges,
                           edges_inlet, nodes_inlet, ps, filter_target);
  check_topology();
  return changed;
}

inline Uint Topology::remove_beyond(const int exit_dim, const double Ln){
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
  check_topology();
  return count;
}

inline bool Topology::resize(const double ds){
  return resizing(edges, node2edges, ps, ds);
}

inline void Topology::write_checkpoint(const std::string& checkpointsfolder, const double t, partrac::Params& prm) const {
  prm.set<double>("t", t);
  // refinement and coarsening change the count, so refresh it here
  prm.set<Uint>("Nrw_current", ps.N());
  prm.dump(checkpointsfolder);
  // dump_positions(checkpointsfolder + "/positions.pos", ps.x_rw, ps.Nrw);
  ps.dump_positions(checkpointsfolder + "/positions.pos");
  dump_faces(checkpointsfolder + "/faces.face", faces);
  dump_edges(checkpointsfolder + "/edges.edge", edges);
  //dump_colors(checkpointsfolder + "/colors.col", ps.c_rw, ps.Nrw);
  ps.dump_scalar(checkpointsfolder + "/colors.col", "c");
  if (prm.get<bool>("inject")){
    dump_vector_field(checkpointsfolder + "/positions_inj.pos", pos_inj);
    dump_edges(checkpointsfolder + "/edges_inj.edge", edges_inj);
    dump_list(checkpointsfolder + "/edges_inlet.list", edges_inlet);
    dump_list(checkpointsfolder + "/nodes_inlet.list", nodes_inlet);
  }
  if (prm.has("local_dt") && prm.get<bool>("local_dt")){
    ps.dump_scalar(checkpointsfolder + "/t_loc.dat", "t_loc");
  }
}

inline void Topology::load_checkpoint(const std::string& checkpointsfolder, const partrac::Params& prm){
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
  if (prm.get<bool>("inject")){
    std::string posinjfile = checkpointsfolder + "/positions_inj.pos";
    load_vector_field(posinjfile, pos_inj);
    std::string edgesinjfile = checkpointsfolder + "/edges_inj.edge";
    load_edges(edgesinjfile, edges_inj);
    std::string edgesinletfile = checkpointsfolder + "/edges_inlet.list";
    std::string nodesinletfile = checkpointsfolder + "/nodes_inlet.list";
    load_list(edgesinletfile, edges_inlet);
    load_list(nodesinletfile, nodes_inlet);
  }
  if (prm.has("local_dt") && prm.get<bool>("local_dt")){
    ps.load_scalar(checkpointsfolder + "/t_loc.dat", "t_loc");
  }
}

inline void Topology::dump_hdf5(H5::H5File& h5f, const std::string& groupname, std::map<std::string, bool>& output_fields){
  
  if (dim() > 0)
    mesh2hdf(h5f, groupname, ps, faces, edges, output_fields["tau"]);
  ps.dump_hdf5(h5f, groupname, output_fields);
}

inline void Topology::load_initial_state(std::shared_ptr<Initializer> init_state, partrac::Params& prm){
  clear();

  edges = init_state->edges;
  faces = init_state->faces;

  if (init_state->clear_initial_edges){
    std::cout << "Clearing initial edges!" << std::endl;
    clear();
  }
  if (init_state->inject){
    // Injection raises what the inlet traces by a dimension, so a sheet inlet
    // would sweep a volume. What clear_initial_edges leaves is the inlet.
    if (faces.size() > 0){
      std::cerr << "Error: an inlet with faces would sweep a volume. Stopping."
                << std::endl;
      exit(1);
    }
    pos_inj = init_state->nodes;
    edges_inj = edges;
    for (Uint i=0; i<pos_inj.size(); ++i){
      nodes_inlet.push_back(i);
    }
    for (Uint i=0; i<edges_inj.size(); ++i){
      edges_inlet.push_back(i);
    }
    // Refinement splits the inlet edge with its template, so a ds_max under its
    // spacing refines the injected curve instead of stalling
    double ds_inj_max = 0.;
    for ( const auto & edge : edges_inj )
      ds_inj_max = std::max(ds_inj_max,
                            (pos_inj[edge.first[0]] - pos_inj[edge.first[1]]).norm());
    if (ds_max > 0. && ds_inj_max > ds_max)
      std::cerr << "Warning: ds_max = " << ds_max << " is below the injection "
                << "template's longest edge " << ds_inj_max << ". The inlet "
                << "will be refined to match, so the injected curve gets finer "
                << "as the run goes on." << std::endl;
  }
  ps.add(init_state->nodes, 0);
  // Nrw is the request; these are what was placed
  prm.set<Uint>("Nrw_init", ps.N());
  prm.set<Uint>("Nrw_current", ps.N());
  check_topology();
}

template<typename T>
void Topology::write_statistics( std::ofstream &statfile
                               , const double t
                               , const double ds_max
                               , T& integrator){
                           //const bool do_dump_hist,
                           //const std::string histfolder,
                           //std::shared_ptr<Integrator> integrator){
  write_stats_row(statfile, mesh_stats_columns(t, ps, faces, edges, ds_max,
                                               integrator.get_accepted(),
                                               integrator.get_declined(),
                                               dim_settled()));
              //integrator->get_accepted(), integrator->get_declined());
}

#endif
