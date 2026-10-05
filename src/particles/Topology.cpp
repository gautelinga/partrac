#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include "Error.hpp"
#include "Topology.hpp"
#include "mesh.hpp"
#include "Initializer.hpp"
#include "io.hpp"

TopologyOptions cloud_options(const partrac::Params& prm){
  TopologyOptions o;
  o.ds_min = prm.get<double>("ds_min");
  o.ds_max = prm.get<double>("ds_max");
  o.injecting = prm.has("inject") ? prm.get<bool>("inject") : false;
  o.verbose = prm.get<bool>("verbose");
  o.Dm = prm.has("Dm") ? prm.get<double>("Dm") : 0.;
  return o;
}

TopologyOptions mesh_options(const partrac::Params& prm){
  TopologyOptions o = cloud_options(prm);
  o.curv_refine_factor = prm.get<double>("curv_refine_factor");
  o.cut_if_stuck = prm.get<bool>("cut_if_stuck");
  o.inject_edges = prm.get<bool>("inject_edges");
  o.filter_target = prm.get<int>("filter_target");
  return o;
}

Topology::Topology(ParticleSet& ps, const TopologyOptions& opts) : ps(ps), opts(opts) {}

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
      if (opts.Dm > 0.)
        std::cerr << "Warning: Dm = " << opts.Dm << " > 0 on a " << d
                  << "-dimensional mesh. Brownian motion and material "
                  << "deformation are normally not compatible." << std::endl;
    }
    return;
  }
  if (d == dim0) return;
  // Injection may raise the dimension
  if (d > dim0 && opts.injecting){
    dim0 = d;
    return;
  }
  partrac::fail("the mesh changed dimension, from ", dim0, " to ", d, ", with ", ps.N(),
                " nodes and ", edges.size(), " edges left");
}

// Dimension the run settles into, including injection
int Topology::dim_settled(){
  int d = dim();
  if (opts.injecting && opts.inject_edges){
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
      remove_inactive(*this, face_isactive, edge_isactive, node_isactive, ps);
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
      remove_inactive(*this, face_isactive, edge_isactive, node_isactive, ps);
      check_topology();
    }
  }
}

Uint Topology::refine(){
  Uint n = refinement(*this, ps, opts.ds_max, opts.curv_refine_factor,
                      opts.cut_if_stuck ? StuckEdge::Cut : StuckEdge::Stop);
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
  const double ds_cut = full ? opts.ds_min : 1e-12 * ds_longest;

  // Skip if no edge is below the cut
  Uint n = 0;
  if (ds_shortest <= ds_cut * (1. + 1e-9))
    n = coarsening(*this, ps, ds_cut, opts.curv_refine_factor);
  check_topology();
  return n;
}

Uint Topology::inject(){
  Uint n = injection(*this, ps, opts.inject_edges, opts.verbose);
  // Cull what is no longer at the inlet
  std::vector<bool> face_isactive(faces.size(), true);
  std::vector<bool> edge_isactive(edges.size(), true);
  std::vector<bool> node_isactive(ps.N(), true);
  remove_inactive(*this, face_isactive, edge_isactive, node_isactive, ps);
  check_topology();
  return n;
}

void Topology::compute_interior(){
  // Only for curvature-weighted refinement
  if (!computes_curvature())
    return;
  compute_interior_prop(interior_ang, mixed_areas, face_normals, *this, ps);
  compute_mean_curv(*this, ps, interior_ang, mixed_areas, face_normals);
}

void Topology::remove_nodes_safe(std::vector<bool>& node_isactive){
  std::vector<bool> face_isactive(faces.size(), true);
  std::vector<bool> edge_isactive(edges.size(), true);
  remove_inactive(*this, face_isactive, edge_isactive, node_isactive, ps);
  check_dim();
}

bool Topology::filter(){
  bool changed = filtering(*this, ps, opts.filter_target);
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

// Node pairs, stored lengths, tau and rho_prev
static void write_edges(H5::H5File& h5f, const EdgesType& edges, const std::string& suffix){
  const Uint n = edges.size();
  std::vector<Uint> nodes(2*n);
  std::vector<double> dl0(n), tau(n), rho_prev(n);
  for (Uint i = 0; i < n; ++i){
    nodes[2*i] = edges[i].first[0];
    nodes[2*i+1] = edges[i].first[1];
    dl0[i] = edges[i].second;
    tau[i] = edges[i].tau;
    rho_prev[i] = edges[i].rho_prev;
  }
  ulong2hdf5(h5f, "edges" + suffix, nodes, n, 2);
  scalar2hdf5(h5f, "dl0" + suffix, dl0, n);
  scalar2hdf5(h5f, "edge_tau" + suffix, tau, n);
  scalar2hdf5(h5f, "edge_rho_prev" + suffix, rho_prev, n);
}

static void read_edges(const H5::H5File& h5f, EdgesType& edges, const std::string& suffix, const Uint n_nodes){
  const Uint n = hdf5_rows(h5f, "edges" + suffix);
  std::vector<Uint> nodes;
  std::vector<double> dl0, tau, rho_prev;
  hdf52ulongs(h5f, "edges" + suffix, n, 2, nodes);
  hdf52doubles(h5f, "dl0" + suffix, n, 1, dl0);
  hdf52doubles(h5f, "edge_tau" + suffix, n, 1, tau);
  hdf52doubles(h5f, "edge_rho_prev" + suffix, n, 1, rho_prev);
  for (Uint i = 0; i < n; ++i){
    if (nodes[2*i] >= n_nodes || nodes[2*i+1] >= n_nodes)
      partrac::fail(h5f.getFileName(), ": edge ", i, " of 'edges", suffix, "' joins nodes ", nodes[2*i],
                    " and ", nodes[2*i+1], ", of ", n_nodes);
    edges.push_back({{nodes[2*i], nodes[2*i+1]}, dl0[i], tau[i], rho_prev[i]});
  }
}

// Edge triples, stored areas, tau and rho_prev
static void write_faces(H5::H5File& h5f, const FacesType& faces){
  const Uint n = faces.size();
  std::vector<Uint> tri(3*n);
  std::vector<double> dA0(n), tau(n), rho_prev(n);
  for (Uint i = 0; i < n; ++i){
    for (Uint j = 0; j < 3; ++j)
      tri[3*i+j] = faces[i].first[j];
    dA0[i] = faces[i].second;
    tau[i] = faces[i].tau;
    rho_prev[i] = faces[i].rho_prev;
  }
  ulong2hdf5(h5f, "face_edges", tri, n, 3);
  scalar2hdf5(h5f, "dA0", dA0, n);
  scalar2hdf5(h5f, "face_tau", tau, n);
  scalar2hdf5(h5f, "face_rho_prev", rho_prev, n);
}

static void read_faces(const H5::H5File& h5f, FacesType& faces, const Uint n_edges){
  const Uint n = hdf5_rows(h5f, "face_edges");
  std::vector<Uint> tri;
  std::vector<double> dA0, tau, rho_prev;
  hdf52ulongs(h5f, "face_edges", n, 3, tri);
  hdf52doubles(h5f, "dA0", n, 1, dA0);
  hdf52doubles(h5f, "face_tau", n, 1, tau);
  hdf52doubles(h5f, "face_rho_prev", n, 1, rho_prev);
  for (Uint i = 0; i < n; ++i){
    for (Uint j = 0; j < 3; ++j)
      if (tri[3*i+j] >= n_edges)
        partrac::fail(h5f.getFileName(), ": face ", i, " names edge ", tri[3*i+j], ", of ", n_edges);
    faces.push_back({{tri[3*i], tri[3*i+1], tri[3*i+2]}, dA0[i], tau[i], rho_prev[i]});
  }
}

// Indices below n
static void read_list(const H5::H5File& h5f, const std::string& name, std::vector<Uint>& li,
                      const Uint n, const std::string& what){
  hdf52ulongs(h5f, name, hdf5_rows(h5f, name), 1, li);
  for (Uint i = 0; i < li.size(); ++i)
    if (li[i] >= n)
      partrac::fail(h5f.getFileName(), ": entry ", i, " of '", name, "' names ", what, " ", li[i], ", of ", n);
}

// t in a params file, or NaN
static double params_t(const std::string& path){
  std::ifstream f(path);
  std::string line;
  while (std::getline(f, line))
    if (line.rfind("t=", 0) == 0){
      char* end = nullptr;
      const double t = std::strtod(line.c_str() + 2, &end);
      if (end != line.c_str() + 2) return t;
    }
  return std::numeric_limits<double>::quiet_NaN();
}

bool Topology::has_t_loc(const partrac::Params& prm) const {
  return records_t_loc || (prm.has("local_dt") && prm.get<bool>("local_dt"));
}

// Written beside, then moved over the last: checkpoint.h5 first, params.dat last
void Topology::write_checkpoint(const std::string& checkpointsfolder, const double t, partrac::Params& prm) const {
  prm.set<double>("t", t);
  prm.set<Uint>("Nrw_current", ps.N());
  prm.dump_tmp(checkpointsfolder);
  const std::string path = checkpointsfolder + "/checkpoint.h5";
  try {
    H5::H5File h5f(path + ".tmp", H5F_ACC_TRUNC);
    H5::Attribute at = h5f.createAttribute("t", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(H5S_SCALAR));
    at.write(H5::PredType::NATIVE_DOUBLE, &t);
    at.close();
    ps.write_checkpoint(h5f, has_t_loc(prm));
    write_edges(h5f, edges, "");
    write_faces(h5f, faces);
    if (records_doublings){
      std::vector<Uint> d(doublings);
      d.resize(edges.size(), 0);
      ulong2hdf5(h5f, "doublings", d, d.size());
    }
    if (prm.get<bool>("inject")){
      vector2hdf5(h5f, "positions_inj", pos_inj, pos_inj.size());
      write_edges(h5f, edges_inj, "_inj");
      ulong2hdf5(h5f, "edges_inlet", edges_inlet, edges_inlet.size());
      ulong2hdf5(h5f, "nodes_inlet", nodes_inlet, nodes_inlet.size());
    }
    h5f.close();
  } catch (const H5::Exception&){
    partrac::fail("cannot write the checkpoint ", path, ".tmp");
  }
  if (std::rename((path + ".tmp").c_str(), path.c_str()) != 0)
    partrac::fail("cannot move ", path, ".tmp to ", path, ": ", std::strerror(errno));
  prm.commit_dump(checkpointsfolder);
}

// checkpoint.h5, or the text files of older runs
void Topology::load_checkpoint(const std::string& checkpointsfolder, const partrac::Params& prm){
  const std::string path = checkpointsfolder + "/checkpoint.h5";
  if (!std::ifstream(path)){
    if (!std::ifstream(checkpointsfolder + "/positions.pos"))
      partrac::fail(checkpointsfolder, ": no checkpoint.h5 or positions.pos");
    load_text_checkpoint(checkpointsfolder, prm);
    return;
  }
  H5::Exception::dontPrint();
  H5::H5File h5f;
  double t = 0.;
  try {
    h5f.openFile(path, H5F_ACC_RDONLY);
    if (H5Aexists(h5f.getId(), "t") <= 0)
      partrac::fail(path, ": no attribute 't'");
    h5f.openAttribute("t").read(H5::PredType::NATIVE_DOUBLE, &t);
  } catch (const H5::Exception&){
    partrac::fail("cannot read the checkpoint ", path);
  }
  // One checkpoint
  if (t != prm.get<double>("t"))
    partrac::fail(path, " is at t = ", t, ", its params.dat at t = ", prm.get<double>("t"),
                  params_t(checkpointsfolder + "/params.dat.tmp") == t
                  ? "; params.dat.tmp beside it is the params.dat written with it" : "");
  ps.read_checkpoint(h5f, has_t_loc(prm));
  read_edges(h5f, edges, "", ps.N());
  read_faces(h5f, faces, edges.size());
  if (records_doublings)
    hdf52ulongs(h5f, "doublings", edges.size(), 1, doublings);
  if (prm.get<bool>("inject")){
    const Uint n = hdf5_rows(h5f, "positions_inj");
    hdf52vectors(h5f, "positions_inj", n, pos_inj);
    read_edges(h5f, edges_inj, "_inj", n);
    read_list(h5f, "edges_inlet", edges_inlet, edges.size(), "edge");
    read_list(h5f, "nodes_inlet", nodes_inlet, ps.N(), "node");
  }
}

void Topology::load_text_checkpoint(const std::string& checkpointsfolder, const partrac::Params& prm){
  std::string posfile = checkpointsfolder + "/positions.pos";
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
  if (has_t_loc(prm)){
    ps.load_scalar(checkpointsfolder + "/t_loc.dat", "t_loc");
  }
  if (ps.carries() == TransportElement::Vector){
    ps.load_vector(checkpointsfolder + "/rhohat.vec", "rhohat");
    ps.load_scalar(checkpointsfolder + "/w.dat", "w");
    ps.load_scalar(checkpointsfolder + "/S.dat", "S");
  }
  // F factored, or whole in a checkpoint from before the factors
  if (ps.carries() == TransportElement::Tensor){
    if (std::ifstream(checkpointsfolder + "/Q.ten")){
      ps.load_tensor(checkpointsfolder + "/Q.ten", "Q");
      ps.load_vector(checkpointsfolder + "/logstretch.vec", "logstretch");
      ps.load_vector(checkpointsfolder + "/U.vec", "U");
    }
    else
      ps.load_tensor(checkpointsfolder + "/F.ten", "F");
  }
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

void Topology::dump_hdf5(H5::H5File& h5f, const std::string& groupname, const OutputFields& output_fields){
  
  if (dim() > 0)
    mesh2hdf(h5f, groupname, ps, faces, edges, output_fields.tau, records_doublings ? &edge_doublings() : nullptr);
  ps.dump_hdf5(h5f, groupname, output_fields);
}

void Topology::load_initial_state(const InitialState& init_state, partrac::Params& prm){
  clear();

  edges = init_state.edges;
  faces = init_state.faces;

  if (prm.get<bool>("clear_initial_edges")){
    std::cout << "Clearing initial edges!" << std::endl;
    clear();
  }
  if (prm.get<bool>("inject")){
    // A sheet inlet would sweep a volume
    if (faces.size() > 0){
      partrac::fail("an inlet with faces would sweep a volume");
    }
    pos_inj = init_state.nodes;
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
    if (opts.ds_max > 0. && ds_inj_max > opts.ds_max)
      std::cerr << "Warning: ds_max = " << opts.ds_max << " is below the injection "
                << "template's longest edge " << ds_inj_max << ". The inlet "
                << "will be refined to match, so the injected curve gets finer "
                << "as the run goes on." << std::endl;
  }
  ps.add(init_state.nodes, 0);
  // Actual counts; Nrw is the request
  prm.set<Uint>("Nrw_init", ps.N());
  prm.set<Uint>("Nrw_current", ps.N());
  check_topology();
}
