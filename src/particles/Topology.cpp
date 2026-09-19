#include <algorithm>
#include <iostream>
#include <limits>
#include "Error.hpp"
#include "Topology.hpp"
#include "mesh.hpp"
#include "Initializer.hpp"
#include "io.hpp"

Topology::Topology(ParticleSet& ps, const partrac::Params& prm) : ps(ps) {
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

// Stop on an empty set or a changed dimension
void Topology::check_topology(){
  if (ps.N() == 0){
    partrac::fail("no particles left");
  }
  check_dim();
}

void Topology::check_dim(){
  if (ps.N() == 0) return;
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
  // Injection may raise the dimension
  if (d > dim0 && injecting){
    dim0 = d;
    return;
  }
  partrac::fail("the mesh changed dimension, from ", dim0, " to ", d, ", with ", ps.N(),
                " nodes and ", edges.size(), " edges left");
}

// Dimension the run settles into, including injection
int Topology::dim_settled(){
  int d = dim();
  if (injecting && inject_edges){
    const int inlet = edges_inj.size() > 0 ? 1
                    : (pos_inj.size() > 0 ? 0 : -1);
    if (inlet >= 0)
      d = std::max(d, inlet + 1);
  }
  return d;
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
    // Independent per edge
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
      check_topology();   // culling may drop the dimension
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
        continue;                 // degenerate, culled later
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

Uint Topology::refine(){
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

// When not full, collapse only numerically zero edges (degenerate medians)
Uint Topology::coarsen(const bool full){
  // Shortest and longest edge
  double ds_shortest = std::numeric_limits<double>::infinity();
  double ds_longest = 0.;
  #pragma omp parallel for reduction(min:ds_shortest) reduction(max:ds_longest)
  for (Uint i = 0; i < edges.size(); ++i){
    const double ds = ps.dist(edges[i].first[0], edges[i].first[1]);
    ds_shortest = std::min(ds_shortest, ds);
    ds_longest = std::max(ds_longest, ds);
  }
  // Cut scaled from the mesh
  const double ds_cut = full ? ds_min : 1e-12 * ds_longest;

  // Skip if no edge is below the cut
  Uint n = 0;
  if (ds_shortest <= ds_cut * (1. + 1e-9))
    n = coarsening(faces, edges,
                   edge2faces, node2edges,
                   edges_inlet, nodes_inlet,
                   ps, ds_cut,
                   curv_refine_factor);
  check_topology();
  return n;
}

Uint Topology::inject(){
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
  // Cull what is no longer at the inlet
  std::vector<bool> face_isactive(faces.size(), true);
  std::vector<bool> edge_isactive(edges.size(), true);
  std::vector<bool> node_isactive(ps.N(), true);
  remove_inactive(faces, edges, edge2faces, node2edges,
                  edges_inlet, nodes_inlet,
                  face_isactive, edge_isactive, node_isactive, ps);
  check_topology();
  return n;
}

void Topology::compute_interior(){
  // Only for curvature-weighted refinement
  if (!computes_curvature())
    return;
  compute_interior_prop(interior_ang, mixed_areas, face_normals,
                        faces, edges, edge2faces, ps);
  compute_mean_curv(faces, edges, edge2faces, node2edges,
                    ps, interior_ang, mixed_areas, face_normals);
}

void Topology::remove_nodes_safe(std::vector<bool>& node_isactive){
  std::vector<bool> face_isactive(faces.size(), true);
  std::vector<bool> edge_isactive(edges.size(), true);
  remove_inactive(faces, edges,
                  edge2faces, node2edges,
                  edges_inlet, nodes_inlet,
                  face_isactive, edge_isactive, node_isactive,
                  ps);
  check_dim();
}

bool Topology::filter(){
  bool changed = filtering(faces, edges, edge2faces, node2edges,
                           edges_inlet, nodes_inlet, ps, filter_target);
  check_topology();
  return changed;
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
  check_topology();
  return count;
}

bool Topology::resize(const double ds){
  return resizing(edges, node2edges, ps, ds);
}

bool Topology::resize_doublings(const double ds){
  return resizing_doublings(edges, edge_doublings(), ps, ds);
}

void Topology::write_checkpoint(const std::string& checkpointsfolder, const double t, partrac::Params& prm) const {
  prm.set<double>("t", t);
  prm.set<Uint>("Nrw_current", ps.N());
  prm.dump(checkpointsfolder);
  // dump_positions(checkpointsfolder + "/positions.pos", ps.x_rw, ps.Nrw);
  ps.dump_positions(checkpointsfolder + "/positions.pos");
  dump_faces(checkpointsfolder + "/faces.face", faces);
  dump_edges(checkpointsfolder + "/edges.edge", edges);
  if (records_doublings){
    std::vector<Uint> d(doublings);
    d.resize(edges.size(), 0);
    dump_list(checkpointsfolder + "/doublings.list", d);
  }
  //dump_colors(checkpointsfolder + "/colors.col", ps.c_rw, ps.Nrw);
  ps.dump_scalar(checkpointsfolder + "/colors.col", "c");
  if (prm.get<bool>("inject")){
    dump_vector_field(checkpointsfolder + "/positions_inj.pos", pos_inj);
    dump_edges(checkpointsfolder + "/edges_inj.edge", edges_inj);
    dump_list(checkpointsfolder + "/edges_inlet.list", edges_inlet);
    dump_list(checkpointsfolder + "/nodes_inlet.list", nodes_inlet);
  }
  if (records_t_loc || (prm.has("local_dt") && prm.get<bool>("local_dt"))){
    ps.dump_scalar(checkpointsfolder + "/t_loc.dat", "t_loc");
  }
  // Carried fields and ids
  if (ps.carries() == TransportElement::Vector){
    ps.dump_vector(checkpointsfolder + "/rhohat.vec", "rhohat");
    ps.dump_scalar(checkpointsfolder + "/w.dat", "w");
    ps.dump_scalar(checkpointsfolder + "/S.dat", "S");
  }
  if (ps.carries() == TransportElement::Tensor)
    ps.dump_tensor(checkpointsfolder + "/F.ten", "F");
  if (ps.records_generation())
    ps.dump_scalar(checkpointsfolder + "/generation.dat", "generation");
  ps.dump_ids(checkpointsfolder + "/id.list");
}

void Topology::load_checkpoint(const std::string& checkpointsfolder, const partrac::Params& prm){
  std::string posfile = checkpointsfolder + "/positions.pos";
  //load_positions(posfile, pos_init, prm.Nrw);
  ps.load_positions(posfile);
  std::string facefile = checkpointsfolder + "/faces.face";
  load_faces(facefile, faces);
  std::string edgefile = checkpointsfolder + "/edges.edge";
  load_edges(edgefile, edges);
  if (records_doublings){
    load_list(checkpointsfolder + "/doublings.list", doublings);
    doublings.resize(edges.size(), 0);
  }
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
  if (records_t_loc || (prm.has("local_dt") && prm.get<bool>("local_dt"))){
    ps.load_scalar(checkpointsfolder + "/t_loc.dat", "t_loc");
  }
  if (ps.carries() == TransportElement::Vector){
    ps.load_vector(checkpointsfolder + "/rhohat.vec", "rhohat");
    ps.load_scalar(checkpointsfolder + "/w.dat", "w");
    ps.load_scalar(checkpointsfolder + "/S.dat", "S");
  }
  if (ps.carries() == TransportElement::Tensor)
    ps.load_tensor(checkpointsfolder + "/F.ten", "F");
  if (ps.records_generation())
    ps.load_scalar(checkpointsfolder + "/generation.dat", "generation");
  ps.load_ids(checkpointsfolder + "/id.list");
}

bool Topology::sort_by_cell(){
  const std::vector<Uint> old2new = ps.sort_by_cell();
  if (old2new.empty())
    return false;
  renumber(old2new);
  return true;
}

void Topology::shuffle(std::mt19937& gen){
  const std::vector<Uint> old2new = ps.shuffle(gen);
  if (!old2new.empty())
    renumber(old2new);
}

// Remap node indices
void Topology::renumber(const std::vector<Uint>& old2new){
  for (auto & edge : edges)
    for (Uint j = 0; j < 2; ++j)
      edge.first[j] = old2new[edge.first[j]];
  for (auto & inode : nodes_inlet)
    inode = old2new[inode];
  compute_node2edges(node2edges, edges, ps.N());
}

void Topology::dump_hdf5(H5::H5File& h5f, const std::string& groupname, std::map<std::string, bool>& output_fields){
  
  if (dim() > 0)
    mesh2hdf(h5f, groupname, ps, faces, edges, output_fields["tau"], records_doublings ? &edge_doublings() : nullptr);
  ps.dump_hdf5(h5f, groupname, output_fields);
}

void Topology::load_initial_state(std::shared_ptr<Initializer> init_state, partrac::Params& prm){
  clear();

  edges = init_state->edges;
  faces = init_state->faces;

  if (init_state->clear_initial_edges){
    std::cout << "Clearing initial edges!" << std::endl;
    clear();
  }
  if (init_state->inject){
    // A sheet inlet would sweep a volume
    if (faces.size() > 0){
      partrac::fail("an inlet with faces would sweep a volume");
    }
    pos_inj = init_state->nodes;
    edges_inj = edges;
    for (Uint i=0; i<pos_inj.size(); ++i){
      nodes_inlet.push_back(i);
    }
    for (Uint i=0; i<edges_inj.size(); ++i){
      edges_inlet.push_back(i);
    }
    // Inlet is refined along with its template
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
  // Actual counts; Nrw is the request
  prm.set<Uint>("Nrw_init", ps.N());
  prm.set<Uint>("Nrw_current", ps.N());
  check_topology();
}
