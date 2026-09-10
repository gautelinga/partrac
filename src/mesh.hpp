#ifndef __MESH_HPP
#define __MESH_HPP

#include <map>
#include "typedefs.hpp"
#include "utils.hpp"
#include "ParticleSet.hpp"

// using namespace std;
// Declarations
void compute_node2edges(Node2EdgesType&, const EdgesType&, const Uint);
void compute_edge2faces(Edge2FacesType&, const FacesType&, const EdgesType&);
void remove_faces(FacesType&, const std::vector<bool>&);
void remove_edges(FacesType&, EdgesType&, const std::vector<bool>&, EdgesListType&);
void remove_inactive(FacesType&, EdgesType&, Edge2FacesType&, Node2EdgesType&,
                     EdgesListType&, NodesListType&, std::vector<bool>&,
                     std::vector<bool>&, std::vector<bool>&, ParticleSet&);

// Definitions

void mesh2hdf( H5::H5File& h5f, const std::string& groupname
             , const ParticleSet& ps
             , const FacesType& faces
             , const EdgesType& edges
             , const bool output_tau
             ){
  // This function is here because it contains ps. Consider stripping that.
  // Faces
  if (faces.size() > 0){
    hsize_t faces_dims[2];
    faces_dims[0] = faces.size();
    faces_dims[1] = 3;
    H5::DataSpace faces_dspace(2, faces_dims);
    std::vector<Uint> faces_arr(faces_dims[0]*faces_dims[1]);

    std::vector<double> dA(faces_dims[0]);
    std::vector<double> dA0(faces_dims[0]);
    std::vector<double> tau_(faces_dims[0]);
    
    for (Uint iface=0; iface < faces_dims[0]; ++iface){
      std::set<Uint> unique_nodes;
      for (Uint i=0; i<3; ++i){
        Uint iedge = faces[iface].first[i];
        for (Uint j=0; j<2; ++j){
          Uint inode = edges[iedge].first[j];
          unique_nodes.insert(inode);
        }
      }
      Uint count = 0;
      for ( auto & inode : unique_nodes ){
        faces_arr[iface*faces_dims[1] + count] = inode;
        ++count;
      }
      dA[iface] = ps.triangle_area(iface, faces, edges);
      dA0[iface] = faces[iface].second;
      if (output_tau){
        tau_[iface] = faces[iface].tau;
      }
    }
    H5::DataSet faces_dset = h5f.createDataSet(groupname + "/faces",
                                           H5::PredType::NATIVE_ULONG,
                                           faces_dspace);
    faces_dset.write(faces_arr.data(), H5::PredType::NATIVE_ULONG);

    scalar2hdf5(h5f, groupname + "/dA", dA, faces_dims[0]);
    scalar2hdf5(h5f, groupname + "/dA0", dA0, faces_dims[0]);

    if (output_tau)
      scalar2hdf5(h5f, groupname + "/tau", tau_, faces_dims[0]);
  }
  // Edges
  else if (edges.size() > 0){
    hsize_t edges_dims[2];
    edges_dims[0] = edges.size();
    edges_dims[1] = 2;
    H5::DataSpace edges_dspace(2, edges_dims);
    std::vector<Uint> edges_arr(edges_dims[0]*edges_dims[1]);

    std::vector<double> dl(edges_dims[0]);
    std::vector<double> dl0(edges_dims[0]);
    std::vector<double> tau_(edges_dims[0]);
    for (Uint iedge=0; iedge < edges_dims[0]; ++iedge){
      for (Uint j=0; j<2; ++j){
        edges_arr[iedge*edges_dims[1] + j] = edges[iedge].first[j];
      }
      dl[iedge] = ps.dist(edges[iedge].first[0], edges[iedge].first[1]);
      dl0[iedge] = edges[iedge].second;
      if (output_tau)
        tau_[iedge] = edges[iedge].tau;
    }

    H5::DataSet edges_dset = h5f.createDataSet(groupname + "/edges",
                                             H5::PredType::NATIVE_ULONG,
                                             edges_dspace);
    edges_dset.write(edges_arr.data(), H5::PredType::NATIVE_ULONG);

    scalar2hdf5(h5f, groupname + "/dl", dl, edges_dims[0]);
    scalar2hdf5(h5f, groupname + "/dl0", dl0, edges_dims[0]);
    if (output_tau)
      scalar2hdf5(h5f, groupname + "/tau", tau_, edges_dims[0]);
  }
}

std::array<Uint, 2> sort_edges(Uint inode, Uint kedge, Uint ledge,
                          EdgesType &edges){
  if (inode == edges[kedge].first[0] || inode == edges[kedge].first[1]){
    return {kedge, ledge};
  }
  return {ledge, kedge};
}

template<typename Row>
Uint get_common_entry(const Row& iedges, const Row& jedges){
  std::vector<Uint> out;
  std::set_intersection(iedges.begin(), iedges.end(), jedges.begin(), jedges.end(),
                        std::back_inserter(out));
  assert(out.size() == 1);
  return *out.begin();
}

Uint get_common_entry(Uint kedge, Uint ledge,
                     EdgesType &edges){
  for (Uint i=0; i<2; ++i){
    Uint knode = edges[kedge].first[i];
    for (Uint j=0; j<2; ++j){
      if (knode == edges[ledge].first[j]){
        return knode;
      }
    }
  }
  std::cout << "Error: No common entry!" << std::endl;
  exit(1);
  return -1;
}

std::array<Uint, 3> get_close_entities(Uint iedge, Uint jedge, Uint kedge, Uint ledge,
                                       EdgesType &edges){
  Uint inode = edges[iedge].first[0];
  Uint knode;
  std::array<Uint, 2> mnedges;
  if (iedge == jedge){
    knode = get_common_entry(kedge, ledge, edges);
    mnedges = sort_edges(inode, kedge, ledge, edges);
  }
  else if (iedge == kedge){
    knode = get_common_entry(ledge, jedge, edges);
    mnedges = sort_edges(inode, ledge, jedge, edges);
  }
  else if (iedge == ledge){
    knode = get_common_entry(jedge, kedge, edges);
    mnedges = sort_edges(inode, jedge, kedge, edges);
  }
  else {
    std::cout << "Error: Found no close entities." << std::endl;
    exit(1);
  }
  // std::cout << mnedges[0] << " " << mnedges[1] << std::endl;
  return {knode, mnedges[0], mnedges[1]};
}


Uint sheet_refinement(FacesType &faces,
                      EdgesType &edges,
                      Edge2FacesType &edge2faces,
                      Node2EdgesType &node2edges,                    
                      EdgesListType &edges_inlet,                
                      NodesListType &nodes_inlet,
                      ParticleSet& ps,
                      const double ds_max,
                      const double curv_refine_factor,
                      const bool cut_if_stuck,
                      const bool check_if_inside=true){
  bool changed;
  Uint n_add = 0;
  std::set<Uint> edges_to_remove;

  do {
    changed = false;

    std::vector<double> ds_ratio_;
    for ( auto & edge : edges ){
      Uint inode = edge.first[0];
      Uint jnode = edge.first[1];
      double ds = ps.dist(inode, jnode);
      // std::cout << inode << " " << jnode << " " << ds << std::endl;
      //double kappa = 0.5*(abs(H_rw[inode]) + abs(H_rw[jnode]));
      //double ds_max_loc = ds_max/(1.0 + curv_refine_factor*kappa);
      double ds_max_loc = ds_max;
      ds_ratio_.push_back(ds/ds_max_loc);
    }

    // Don't refine inlet edges
    for ( auto & inode : edges_inlet ){
      ds_ratio_[inode] = 0.0;
    }

    // Only the edges over threshold are split, so only they need ordering:
    // sorting every edge to find the few percent that qualify was half of
    // refinement. A stable sort of the subset keeps exactly the order the
    // full one gave, so the splits happen in the same order as before.
    std::vector<size_t> ids_;
    for (size_t iedge = 0; iedge < ds_ratio_.size(); ++iedge){
      if (ds_ratio_[iedge] >= 1.0)
        ids_.push_back(iedge);
    }
    std::stable_sort(ids_.begin(), ids_.end(),
                     [&ds_ratio_](size_t i1, size_t i2){ return ds_ratio_[i1] > ds_ratio_[i2]; });

    for ( auto & iedge : ids_ ){

      if (!ps.has_space()){
        // No more points can fit
        break;
      }
      double ds_ratio = ds_ratio_[iedge];
      if (ds_ratio < 1.0)
        break;

      Uint inode = edges[iedge].first[0];
      Uint jnode = edges[iedge].first[1];
      double ds0 = edges[iedge].second;

      Uint new_inode = ps.N();
      // Add point
      bool added = ps.insert_node_between(inode, jnode, check_if_inside);
      if (added){
        changed = true;
        ++n_add;
        //
        edges[iedge].first[1] = new_inode;
        Uint new_iedge = edges.size();
        edges.push_back({{new_inode, jnode}, ds0/2});
        edge2faces.push_back({});
        // Append new node to node2edges list
        node2edges.push_back({iedge, new_iedge});
        // Modify existing entry
        std::replace(node2edges[jnode].begin(), node2edges[jnode].end(), iedge, new_iedge);

        for (auto itface = edge2faces[iedge].begin();
             itface != edge2faces[iedge].end(); ++itface){
          Uint jedge = faces[*itface].first[0];
          Uint kedge = faces[*itface].first[1];
          Uint ledge = faces[*itface].first[2];
          double dA0 = faces[*itface].second;
          double tau = faces[*itface].tau;
          double rho_prev = faces[*itface].rho_prev;
          Uint new_jedge = edges.size();
          std::array<Uint, 3> close_entities = get_close_entities(iedge, jedge, kedge, ledge, edges);
          Uint knode = close_entities[0];
          Uint medge = close_entities[1];
          Uint nedge = close_entities[2];
          edges.push_back({{new_inode, knode}, ds0/2});  // ds0/2 - or what else?

          edge2faces.push_back({});

          Uint new_iface = faces.size();
          //faces[*itface] = {{iedge, new_jedge, medge}, dA0/2};
          faces[*itface].first = {iedge, new_jedge, medge};
          faces[*itface].second = dA0/2;
          //faces.push_back({{new_iedge, nedge, new_jedge}, dA0/2});
          faces.push_back({{new_iedge, nedge, new_jedge}, dA0/2, tau, rho_prev});

          {
            auto & row = edge2faces[nedge];
            row.erase(std::remove(row.begin(), row.end(), *itface), row.end());
          }
          edge2faces[nedge].push_back(new_iface);
          edge2faces[new_iedge].push_back(new_iface);
          edge2faces[new_jedge].push_back(*itface);
          edge2faces[new_jedge].push_back(new_iface);
          //
          node2edges[new_inode].push_back(new_jedge);
          node2edges[knode].push_back(new_jedge);
        }
      }
      else {
        //std::cout << "Here we should remove this edge." << std::endl;
        //exit(1);
        if (cut_if_stuck)
          edges_to_remove.insert(iedge);
      }
    }
  } while (changed);
  if (edges_to_remove.size() > 0){
    //if (!cut_if_stuck){
    //  std::cout << "Edge is stuck! Turn on 'cut_if_stuck' to continue in such cases." << std::endl;
    //  exit(1);
    //}
    std::vector<bool> face_isactive(faces.size(), true);
    std::vector<bool> edge_isactive(edges.size(), true);
    std::vector<bool> node_isactive(ps.N(), true);
    for (auto & jedge : edges_to_remove )
      edge_isactive[jedge] = false;
    remove_inactive(faces, edges, edge2faces, node2edges,
                    edges_inlet, nodes_inlet,
                    face_isactive, edge_isactive, node_isactive, ps);
  }
  return n_add;
}

Uint strip_refinement(FacesType &faces,
                      EdgesType &edges,
                      Edge2FacesType &edge2faces,
                      Node2EdgesType &node2edges,
                      EdgesListType &edges_inlet,
                      NodesListType &nodes_inlet,
                      ParticleSet& ps,
                      const double ds_max,
                      const double curv_refine_factor,
                      const bool cut_if_stuck){
  Uint n_add = 0;
  Uint iedge = 0;
  std::set<Uint> edges_to_remove;
  // Don't refine inlet edges
  std::vector<bool> is_inlet(edges.size(), false);
  for ( auto & jedge : edges_inlet )
    is_inlet[jedge] = true;
  while (iedge < edges.size()){
    if (iedge < is_inlet.size() && is_inlet[iedge]){
      ++iedge;
      continue;
    }
    Uint inode = edges[iedge].first[0];
    Uint jnode = edges[iedge].first[1];
    double ds0 = edges[iedge].second;
    double ds = ps.dist(inode, jnode);
    double tau = edges[iedge].tau;
    double rho_prev = edges[iedge].rho_prev;
    //double kappa = sqrt(abs(H_rw[inode]*H_rw[jnode]));
    //double kappa = 0.5*(abs(H_rw[inode]) + abs(H_rw[jnode]));
    double ds_max_loc = ds_max/(1.0);  // + curv_refine_factor*kappa);
    if (ds > ds_max_loc && ps.has_space()){
      Uint new_inode = ps.N();
      /*bool added = append_new_node(inode, jnode, x_rw, u_rw,
                                   rho_rw, p_rw, c_rw, tau_rw, H_rw, n_rw,
                                   a_rw, Nrw, do_output_all, intp, int_order);*/
      bool added = ps.insert_node_between(inode, jnode);
      if (added){
        //std::cout << "before: " << edges[iedge].first[0] << " " << edges[iedge].first[1] << std::endl;
        edges[iedge].first[1] = new_inode;
        edges[iedge].second = ds0/2;
        Uint new_iedge = edges.size();

        //edges.push_back({{new_inode, jnode}, ds0/2});
        edges.push_back({{new_inode, jnode}, ds0/2, tau, rho_prev});
        //std::cout << "after:  " << edges[iedge].first[0] << " " << edges[iedge].first[1] << std::endl;
        //std::cout << "...and: " << edges[new_iedge].first[0] << " " << edges[new_iedge].first[1] << std::endl;
        edge2faces.push_back({});

        node2edges.push_back({iedge, new_iedge});

        std::replace(node2edges[jnode].begin(), node2edges[jnode].end(), iedge, new_iedge);

        //compute_node2edges(node2edges, edges, new_inode);

        ++n_add;
      }
      else {
        // exit(1);
        edges_to_remove.insert(iedge);
        ++iedge;
      }
    }
    else {
      ++iedge;
    }
  }
  if (edges_to_remove.size() > 0){
    if (!cut_if_stuck){
      std::cout << "Edge is stuck! Turn on 'cut_if_stuck' to continue in such cases." << std::endl;
      exit(1);
    }
    std::vector<bool> face_isactive(faces.size(), true);   // a strip has none
    std::vector<bool> edge_isactive(edges.size(), true);
    std::vector<bool> node_isactive(ps.N(), true);
    for (std::set<Uint>::const_iterator sit = edges_to_remove.begin();
         sit != edges_to_remove.end(); ++sit){
      edge_isactive[*sit] = false;
    }
    remove_inactive(faces, edges, edge2faces, node2edges,
                    edges_inlet, nodes_inlet,
                    face_isactive, edge_isactive, node_isactive, ps);
  }
  return n_add;
}

Uint refinement(FacesType &faces,
                EdgesType &edges,
                Edge2FacesType &edge2faces,
                Node2EdgesType &node2edges,
                EdgesListType &edges_inlet,
                NodesListType &nodes_inlet,
                ParticleSet& ps, const double ds_max,
                const double curv_refine_factor,
                const bool cut_if_stuck){
  Uint n_add = 0;
  if (faces.size() > 0){
    n_add = sheet_refinement(faces, edges, edge2faces, node2edges,
                             edges_inlet, nodes_inlet,
                             ps, ds_max, curv_refine_factor, cut_if_stuck);
  }
  else {
    n_add = strip_refinement(faces, edges, edge2faces, node2edges,
                             edges_inlet, nodes_inlet,
                             ps, ds_max, curv_refine_factor, cut_if_stuck);
  }
  return n_add;
}

void compute_edge2faces(Edge2FacesType &edge2faces,
                        const FacesType &faces,
                        const EdgesType &edges){
  // Rows keep their capacity across rebuilds
  edge2faces.resize(edges.size());
  for (auto & row : edge2faces)
    row.clear();

  for (Uint iface=0; iface < faces.size(); ++iface){
    for (Uint i=0; i < 3; ++i){
      edge2faces[faces[iface].first[i]].push_back(iface);
    }
  }
}

void compute_node2edges(Node2EdgesType &node2edges,
                        const EdgesType &edges,
                        const Uint Nrw){
  node2edges.resize(Nrw);
  for (auto & row : node2edges)
    row.clear();

  for (Uint iedge=0; iedge < edges.size(); ++iedge){
    for (Uint i=0; i < 2; ++i){
      node2edges[edges[iedge].first[i]].push_back(iedge);
    }
  }
}

// GL: Template this
template<typename T>
void print(const T vec){
  for ( auto & elem : vec ){
    std::cout << elem << " ";
  }
  std::cout << std::endl;
}

// Nodes joined to inode by active edges, sorted, and the edge joining each
void get_conodes(std::vector<Uint> &conodes,
                 std::vector<Uint> &coedges,
                 const Uint inode,
                 const Node2EdgesType &node2edges,
                 const EdgesType &edges,
                 const std::vector<bool>& edge_isactive){
  conodes.clear();
  coedges.clear();
  for ( auto & iedge : node2edges[inode] ){
    if (edge_isactive[iedge]){
      for ( auto & jnode : edges[iedge].first ){
        if (jnode != inode){
          conodes.push_back(jnode);
          coedges.push_back(iedge);
        }
      }
    }
  }
  // Sort by node; a doubled edge keeps the later one
  for (Uint a = 1; a < conodes.size(); ++a){
    const Uint n = conodes[a], e = coedges[a];
    Uint b = a;
    for (; b > 0 && conodes[b-1] > n; --b){
      conodes[b] = conodes[b-1];
      coedges[b] = coedges[b-1];
    }
    conodes[b] = n;
    coedges[b] = e;
  }
  Uint m = 0;
  for (Uint a = 0; a < conodes.size(); ++a){
    if (m > 0 && conodes[m-1] == conodes[a])
      coedges[m-1] = coedges[a];
    else {
      conodes[m] = conodes[a];
      coedges[m] = coedges[a];
      ++m;
    }
  }
  conodes.resize(m);
  coedges.resize(m);
}

// The edge joining node to the node get_conodes was called for
Uint coedge_to(const std::vector<Uint> &conodes,
               const std::vector<Uint> &coedges,
               const Uint node){
  for (Uint i = 0; i < conodes.size(); ++i){
    if (conodes[i] == node)
      return coedges[i];
  }
  assert(false);
  return 0;
}

void get_cofaces(std::vector<Uint> &jfaces,
                 const Uint iedge,
                 const Edge2FacesType &edge2faces,
                 const FacesType &faces){
  jfaces.clear();
  for ( auto & iface : edge2faces[iedge] ){
    for ( auto & kedge : faces[iface].first ){
      if (kedge != iedge){
        for (auto & kface : edge2faces[kedge] ){
          if (kface != iface){
            jfaces.push_back(kface);
          }
        }
      }
    }
  }
  std::sort(jfaces.begin(), jfaces.end());
  jfaces.erase(std::unique(jfaces.begin(), jfaces.end()), jfaces.end());
}

bool is_border_node(const Uint inode,
                    const Node2EdgesType &node2edges,
                    const Edge2FacesType &edge2faces){
  for ( auto & iedge : node2edges[inode] ){
    if (edge2faces[iedge].size() == 1){
      return true;
    }
  }
  return false;
}

// Turning angle of the rim at a node
double rim_turn(const Uint inode,
                const ParticleSet& ps,
                const EdgesType &edges,
                const Edge2FacesType &edge2faces,
                const Node2EdgesType &node2edges){
  std::vector<Uint> nbrs;
  for ( auto & iedge : node2edges[inode] ){
    if (edge2faces[iedge].size() < 2)
      nbrs.push_back(edges[iedge].first[0] == inode ? edges[iedge].first[1]
                                                    : edges[iedge].first[0]);
  }
  if (nbrs.size() != 2)
    return 0.;
  Vector3d a = ps.x(inode) - ps.x(nbrs[0]);
  Vector3d b = ps.x(nbrs[1]) - ps.x(inode);
  if (a.norm() <= 0. || b.norm() <= 0.)
    return 0.;
  return acos(std::max(-1., std::min(1., a.dot(b)/(a.norm()*b.norm()))));
}

bool get_new_pos(Vector3d &x,
                 const ParticleSet& ps,
                 const Uint iedge,
                 const EdgesType &edges,
                 const Edge2FacesType &edge2faces,
                 const Node2EdgesType &node2edges){
  Uint inode = edges[iedge].first[0];
  Uint jnode = edges[iedge].first[1];

  bool inode_is_border = is_border_node(inode, node2edges, edge2faces);
  bool jnode_is_border = is_border_node(jnode, node2edges, edge2faces);

  bool both_are_border = inode_is_border && jnode_is_border;
  bool none_are_border = !inode_is_border && !jnode_is_border;

  if (both_are_border && edge2faces[iedge].size() != 1){
    return false;
  }
  else if (none_are_border){
    x = 0.5*(ps.x(inode) + ps.x(jnode));
  }
  else if (both_are_border){
    // Rim nodes merge at their midpoint, unless one is a corner: keep the corner
    const double sharp = 0.25;   // radians; smooth rim curvature is far below
    const double turn_i = rim_turn(inode, ps, edges, edge2faces, node2edges);
    const double turn_j = rim_turn(jnode, ps, edges, edge2faces, node2edges);
    if (std::max(turn_i, turn_j) < sharp)
      x = 0.5*(ps.x(inode) + ps.x(jnode));
    else
      x = turn_i >= turn_j ? ps.x(inode) : ps.x(jnode);
  }
  else if (inode_is_border) {
    x = ps.x(inode);
  }
  else {
    assert(jnode_is_border);
    x = ps.x(jnode);
  }
  return true;
}

// Cross product of a face's first two edges: twice the area by its norm, the
// normal by its direction. With moved, both ends of iedge are read at x.
Vector3d face_cross(const Uint jface,
                    const Uint iedge,
                    const Vector3d &x,
                    const bool moved,
                    const ParticleSet& ps,
                    const FacesType &faces,
                    const EdgesType &edges){
  const Uint inode = edges[iedge].first[0];
  const Uint jnode = edges[iedge].first[1];
  auto pos = [&](const Uint n){
    return (moved && (n == inode || n == jnode)) ? x : ps.x(n);
  };
  const Uint jedge = faces[jface].first[0];
  const Uint kedge = faces[jface].first[1];
  const Vector3d drj = pos(edges[jedge].first[1]) - pos(edges[jedge].first[0]);
  const Vector3d drk = pos(edges[kedge].first[1]) - pos(edges[kedge].first[0]);
  return drj.cross(drk);
}

// Check that no surviving face is turned over or flattened onto its own
// opposite edge when both nodes move to x, and take the areas the move gives
bool normals_are_ok(const Uint iedge,
                    const Vector3d &x,
                    const ParticleSet& ps,
                    const std::vector<Uint> &jfaces,
                    const FacesType &faces,
                    const EdgesType &edges,
                    const std::vector<Vector3d> &cross_old,
                    std::vector<double> &dAs_new){
  dAs_new.resize(jfaces.size());
  for (Uint k = 0; k < jfaces.size(); ++k){
    const Vector3d c = face_cross(jfaces[k], iedge, x, true, ps, faces, edges);
    const double s2 = cross_old[k].squaredNorm();
    if (s2 > 0. && cross_old[k].dot(c) <= 1e-10*s2)
      return false;   // a face with no area has none to lose
    dAs_new[k] = c.norm()/2;
  }
  return true;
}

void get_incident_faces(std::vector<Uint> &kfaces,
                        const Uint iedge,
                        const EdgesType &edges,
                        const Edge2FacesType &edge2faces,
                        const Node2EdgesType &node2edges){
  Uint inode = edges[iedge].first[0];
  Uint jnode = edges[iedge].first[1];
  kfaces.clear();
  for ( auto & kedge : node2edges[inode] ){
    for ( auto & kface : edge2faces[kedge] )
      kfaces.push_back(kface);
  }
  for ( auto & kedge : node2edges[jnode] ){
    for ( auto & kface : edge2faces[kedge] )
      kfaces.push_back(kface);
  }
  std::sort(kfaces.begin(), kfaces.end());
  kfaces.erase(std::unique(kfaces.begin(), kfaces.end()), kfaces.end());
  for ( auto & kface : edge2faces[iedge] )
    kfaces.erase(std::remove(kfaces.begin(), kfaces.end(), kface), kfaces.end());
}

// Scratch for collapse_edge, owned by the sweep
struct CollapseBuffers {
  std::vector<Uint> inodes, jnodes, icoedges, jcoedges, joint_nodes;
  std::vector<Uint> jfaces, kfaces, iedges_vec, ifaces;
  std::vector<std::pair<Uint, Uint>> replace_edges;
  std::vector<Vector3d> cross_old;
  std::vector<double> dAs_old, dAs_new, dA0_new, v_new;
};

bool collapse_edge(const Uint iedge,
                   FacesType &faces,
                   EdgesType &edges,
                   Edge2FacesType &edge2faces,
                   Node2EdgesType &node2edges,
                   std::vector<bool> &face_isactive,
                   std::vector<bool> &edge_isactive,
                   std::vector<bool> &node_isactive,
                   ParticleSet& ps,
                   CollapseBuffers &buf){

  // This function will collapse the edge 'iedge' and thus remove it.
  // The 1-2 facets next to it will be removed.
  // The two nodes it connects to will be replaced by one of them.
  // The dicts node2edges and edge2faces will be updated.

  // Should only be called if edge is active (collapsable).
  assert(edge_isactive[iedge]);

  // Get nodes of that edge
  Uint inode = std::min(edges[iedge].first[0], edges[iedge].first[1]);
  Uint jnode = std::max(edges[iedge].first[0], edges[iedge].first[1]);
  // double ds0 = edges[iedge].second;
  
  // Nodes that are connected to the respective nodes, and the edges that
  // connect them, sorted by node
  auto &inodes = buf.inodes, &jnodes = buf.jnodes;
  auto &icoedges = buf.icoedges, &jcoedges = buf.jcoedges;

  get_conodes(inodes, icoedges, inode, node2edges, edges, edge_isactive);
  get_conodes(jnodes, jcoedges, jnode, node2edges, edges, edge_isactive);
  auto &joint_nodes = buf.joint_nodes;
  joint_nodes.clear();
  
  // Nodes that are connected to both nodes (should be <= 2)
  set_intersection(inodes.begin(), inodes.end(),
                   jnodes.begin(), jnodes.end(),
                   back_inserter(joint_nodes));

  //std::cout << "iedge=" << iedge << std::endl;
  //print(inodes);
  //print(jnodes);
  //print(joint_nodes);

  // Check if it has too many/too few common joint nodes.
  std::size_t num_faces = edge2faces[iedge].size();

  // std::cout << "num_faces=" << num_faces << std::endl;

  if (joint_nodes.size() > num_faces ){
    // Wrong number of joint nodes
    //std::cout << "Wrong number of joint nodes: " << joint_nodes.size() << std::endl;
    return false;
  }

  //if ()

  Vector3d x;
  bool new_node_is_good = get_new_pos(x, ps, iedge, edges, edge2faces, node2edges);
  //bool new_node_is_good = ps.get_new_pos ...?
  if (!new_node_is_good){
    // Corner edge
    // std::cout << "Corner node!" << std::endl;
    return false;
  }
  //std::cout << "New pos: " << x << " " << y << " " << z << std::endl;

  auto &jfaces = buf.jfaces;
  get_cofaces(jfaces, iedge, edge2faces, faces);

  //print(jfaces);

  double dA_res = 0.;
  double dA0_res = 0.;
  for ( auto & iface : edge2faces[iedge] ){
    //dA_res += area(*faceit, x_rw, faces, edges);
    dA_res += ps.triangle_area(iface, faces, edges);
    dA0_res += faces[iface].second;
    //faces[*faceit].second = 0.;
  }

  auto &kfaces = buf.kfaces;
  get_incident_faces(kfaces, iedge, edges, edge2faces, node2edges);
  // Areas of the incident faces now, and their normals for the check below
  auto &cross_old = buf.cross_old;
  auto &dAs_old = buf.dAs_old;
  cross_old.resize(kfaces.size());
  dAs_old.resize(kfaces.size());
  for (Uint k = 0; k < kfaces.size(); ++k){
    cross_old[k] = face_cross(kfaces[k], iedge, x, false, ps, faces, edges);
    dAs_old[k] = cross_old[k].norm()/2;
  }

  // No incident faces to receive the removed mass
  if (kfaces.empty())
    return false;

  // Carry tau across only if every face in the patch has one. All equal:
  // nothing to do. Fresh material (tau = 0) next to mixed material: decline.
  bool share_tau = false;
  {
    const double tau_first = faces[kfaces.front()].tau;
    bool all_equal = true, any_zero = false;
    auto look = [&](const Uint f){
      if (faces[f].tau != tau_first) all_equal = false;
      if (faces[f].tau <= 0.) any_zero = true;
    };
    for ( auto & iface : edge2faces[iedge] ) look(iface);
    for ( auto & kface : kfaces ) look(kface);
    if (!all_equal){
      if (any_zero)
        return false;
      share_tau = true;
    }
  }
  double v_res = 0.;   // variance content of the removed faces
  if (share_tau){
    for ( auto & iface : edge2faces[iedge] )
      v_res += faces[iface].second/sqrt(faces[iface].tau);
  }

  // A collapse may not flip a normal; this also takes the areas after the move
  auto &dAs_new = buf.dAs_new;
  if (!normals_are_ok(iedge, x, ps, kfaces, faces, edges, cross_old, dAs_new))
    return false;

  // assert(jfaces.size()==4 || jfaces.size()==2);

  ps.replace_nodes(x, inode, jnode);

  for (auto & jedge : node2edges[jnode] ){
    for (Uint j=0; j<2; ++j){
      if ( edges[jedge].first[j] == jnode ){
        edges[jedge].first[j] = inode;
      }
    }
  }

  auto &iedges_vec = buf.iedges_vec;
  iedges_vec.clear();
  set_symmetric_difference(node2edges[inode].begin(),
                           node2edges[inode].end(),
                           node2edges[jnode].begin(),
                           node2edges[jnode].end(),
                           back_inserter(iedges_vec));
  // iedges_vec is sorted and unique, being a symmetric difference

  // Doubled edges at the joint nodes (at most two)
  auto &replace_edges = buf.replace_edges;
  replace_edges.clear();
  for (auto & knode : joint_nodes ){
    const Uint ie = coedge_to(inodes, icoedges, knode);
    const Uint je = coedge_to(jnodes, jcoedges, knode);
    Uint kedge_min = std::min(ie, je);
    Uint kedge_max = std::max(ie, je);
    replace_edges.push_back({kedge_max, kedge_min});

    assert(edge_isactive[kedge_max]);
    edge_isactive[kedge_max] = false;
    iedges_vec.erase(std::remove(iedges_vec.begin(), iedges_vec.end(), kedge_max),
                     iedges_vec.end());
  }
  auto replaced = [&](const Uint e){
    for (auto & mit : replace_edges)
      if (mit.first == e) return mit.second;
    return e;
  };

  //assert(edge_isactive[iedge]);
  edge_isactive[iedge] = false;

  for (auto & jface : jfaces ){
    for (Uint j=0; j<3; ++j)
      faces[jface].first[j] = replaced(faces[jface].first[j]);
    //std::cout << "faces[" << *jfaceit << "].first="
    //     << faces[*jfaceit].first[0] << " "
    //     << faces[*jfaceit].first[1] << " "
    //     << faces[*jfaceit].first[2] << std::endl;
  }
  // Distribute the removed mass over the incident faces that survive.
  // The removed face's reference area is dA0_res, and its actual area is dA_res.
  // The ratio r_res = dA0_res/dA_res is used to scale the contribution to the
  // surviving faces based on how much their area changed.
  // The mass is renormalized to conserve the total reference area of the patch.
  const double r_res = dA_res > 0. ? dA0_res/dA_res : 0.; // density of removed face
  double dA0_patch = dA0_res;
  double dA0_est = 0.;
  auto &dA0_new = buf.dA0_new;
  dA0_new.resize(kfaces.size());
  for (Uint k=0; k < kfaces.size(); ++k){
    const double dA0_k = faces[kfaces[k]].second;
    const double w = dAs_new[k] - dAs_old[k];
    dA0_patch += dA0_k;
    if (w > 0.)
      dA0_new[k] = dA0_k + r_res*w; // gets mass proportional to its area gain, with the density of the removed face
    else if (w < 0. && dAs_old[k] > 0.)
      dA0_new[k] = dA0_k*dAs_new[k]/dAs_old[k]; // loses mass proportional to its area loss, with its own density
    else
      dA0_new[k] = dA0_k;
    dA0_est += dA0_new[k];
  }
  // Share out the variance content dA0/sqrt(tau) the same way, from the old
  // dA0 and tau. tau follows from the two shares.
  const double q_res = dA_res > 0. ? v_res/dA_res : 0.;
  double v_patch = v_res;
  double v_est = 0.;
  auto &v_new = buf.v_new;
  v_new.resize(kfaces.size());
  if (share_tau){
    for (Uint k=0; k < kfaces.size(); ++k){
      const double v_k = faces[kfaces[k]].second/sqrt(faces[kfaces[k]].tau);
      const double w = dAs_new[k] - dAs_old[k];
      v_patch += v_k;
      if (w > 0.)
        v_new[k] = v_k + q_res*w;
      else if (w < 0. && dAs_old[k] > 0.)
        v_new[k] = v_k*dAs_new[k]/dAs_old[k];
      else
        v_new[k] = v_k;
      v_est += v_new[k];
    }
    assert (v_est > 0.);
  }

  assert (dA0_est > 0.); // Ensure that the estimated reference area of the patch is positive before renormalizing
  for (Uint k=0; k < kfaces.size(); ++k)
    faces[kfaces[k]].second = dA0_patch*dA0_new[k]/dA0_est;
  if (share_tau){
    for (Uint k=0; k < kfaces.size(); ++k){
      const double v = v_patch*v_new[k]/v_est;
      faces[kfaces[k]].tau = pow(faces[kfaces[k]].second/v, 2);
    }
  }
  // Recompute rho_prev: dA0 and the geometry have changed
  for (Uint k=0; k < kfaces.size(); ++k)
    faces[kfaces[k]].rho_prev = dAs_new[k]/faces[kfaces[k]].second;

  // Deactivate faces adjacent to iedge
  for (auto & jface : edge2faces[iedge] ){
    assert(face_isactive[jface]);
    face_isactive[jface] = false;
  }

  //print(node2edges[std::min(inode, jnode)]);
  node2edges[inode].assign(iedges_vec.begin(), iedges_vec.end());
  //print(node2edges[std::min(inode, jnode)]);
  node2edges[jnode].clear();
  assert(node_isactive[jnode]);
  node_isactive[jnode] = false;
  
  for (auto & knode : joint_nodes){
    // Drop the doubled edges from the joint nodes' lists
    auto & row = node2edges[knode];
    for (auto & mit : replace_edges )
      row.erase(std::remove(row.begin(), row.end(), mit.first), row.end());
    std::sort(row.begin(), row.end());
    row.erase(std::unique(row.begin(), row.end()), row.end());
  }

  for (auto & mit : replace_edges ){
    Uint kedge_max = mit.first;
    Uint kedge_min = mit.second;
    auto &ifaces = buf.ifaces;
    ifaces.clear();
    set_symmetric_difference(edge2faces[kedge_min].begin(),
                             edge2faces[kedge_min].end(),
                             edge2faces[kedge_max].begin(),
                             edge2faces[kedge_max].end(),
                             back_inserter(ifaces));
    //print(edge2faces[kedge_min]);
    //print(edge2faces[kedge_max]);
    edge2faces[kedge_min].assign(ifaces.begin(), ifaces.end());
    //print(edge2faces[kedge_min]);
    edge2faces[kedge_max].clear();
  }

  edge2faces[iedge].clear();

  return true;
}

void remove_faces(FacesType &faces, const std::vector<bool>& face_isactive){
  assert(faces.size() == face_isactive.size());
  // Compact in one pass
  Uint n = 0;
  for (Uint i=0; i < faces.size(); ++i){
    if (face_isactive[i]){
      if (n != i)
        faces[n] = faces[i];
      ++n;
    }
  }
  faces.erase(faces.begin() + n, faces.end());
}

void remove_edges(FacesType &faces, EdgesType &edges,
                  const std::vector<bool> &edge_isactive,
                  EdgesListType &edges_inlet){
  assert(edges.size() == edge_isactive.size());
  // Map old edge indices to new; a removed edge keeps its old index
  std::vector<Uint> new_index(edges.size());
  Uint n = 0;
  for (Uint i=0; i < edges.size(); ++i){
    new_index[i] = i;
    if (edge_isactive[i]){
      if (n != i)
        edges[n] = edges[i];
      new_index[i] = n++;
    }
  }
  edges.erase(edges.begin() + n, edges.end());
  for (auto & face : faces){
    for (Uint j=0; j<3; ++j)
      face.first[j] = new_index[face.first[j]];
  }
  for (auto & iedge : edges_inlet)
    iedge = new_index[iedge];
}

void remove_nodes(EdgesType& edges,
                  NodesListType& nodes_inlet,
                  ParticleSet& ps,
                  const std::vector<bool> &node_isactive
                  ){
  assert(ps.N() == node_isactive.size());
  std::vector<Uint> used_nodes;
  for (Uint i=0; i<ps.N(); ++i){
    if (node_isactive[i])
      used_nodes.push_back(i);
  }
  std::vector<Uint> new_index(ps.N());
  for (Uint i=0; i < ps.N(); ++i)
    new_index[i] = i;
  for (Uint i=0; i<used_nodes.size(); ++i){
    new_index[used_nodes[i]] = i;
    ps.copy_node(i, used_nodes[i]);
  }
  ps.set_N(used_nodes.size());

  for (auto & edge : edges){
    for (Uint j=0; j<2; ++j)
      edge.first[j] = new_index[edge.first[j]];
  }
  for (auto & inode : nodes_inlet)
    inode = new_index[inode];
}

// Remove everything the flags mark dead, and everything that leaves dangling.
// A caller sets the flags it knows and passes all-true for the rest. What
// propagates depends on the dimension: a cloud has no edges, a strip no faces.
void remove_inactive(FacesType &faces, EdgesType &edges,
                     Edge2FacesType &edge2faces, Node2EdgesType &node2edges,
                     EdgesListType &edges_inlet, NodesListType &nodes_inlet,
                     std::vector<bool> &face_isactive,
                     std::vector<bool> &edge_isactive,
                     std::vector<bool> &node_isactive,
                     ParticleSet& ps){
  assert(faces.size() == face_isactive.size());
  assert(edges.size() == edge_isactive.size());
  assert(ps.N() == node_isactive.size());
  const bool has_faces = faces.size() > 0;
  const bool has_edges = edges.size() > 0;

  // The inlet is what the next generation stitches to, so it stays until it is
  // no longer the inlet
  for ( auto & iedge : edges_inlet )
    edge_isactive[iedge] = true;
  for ( auto & inode : nodes_inlet )
    node_isactive[inode] = true;

  // A dead node takes its edges, a dead edge its faces
  for (Uint inode=0; inode < ps.N(); ++inode){
    if (!node_isactive[inode]){
      for ( auto & iedge : node2edges[inode] )
        edge_isactive[iedge] = false;
    }
  }
  if (has_faces){
    for (Uint iedge=0; iedge < edges.size(); ++iedge){
      if (!edge_isactive[iedge]){
        for ( auto & iface : edge2faces[iedge] )
          face_isactive[iface] = false;
      }
    }
  }

  // A face the sweep gave no area is kept while it is at the inlet, and goes
  // as soon as it is not
  if (has_faces){
    std::vector<bool> at_inlet(faces.size(), false);
    for ( auto & iedge : edges_inlet ){
      for ( auto & iface : edge2faces[iedge] )
        at_inlet[iface] = true;
    }
    for (Uint iface=0; iface < faces.size(); ++iface){
      if (faces[iface].second <= 0. && !at_inlet[iface])
        face_isactive[iface] = false;
    }
  }

  // and whatever is left with nothing above it goes too, inlets excepted
  if (has_faces){
    std::vector<bool> is_used(edges.size(), false);
    for (Uint iface=0; iface < faces.size(); ++iface){
      if (face_isactive[iface]){
        for (Uint k=0; k < 3; ++k)
          is_used[faces[iface].first[k]] = true;
      }
    }
    for ( auto & iedge : edges_inlet )
      is_used[iedge] = true;
    for (Uint iedge=0; iedge < edges.size(); ++iedge){
      if (!is_used[iedge])
        edge_isactive[iedge] = false;
    }
  }
  if (has_edges){
    std::vector<bool> is_used(ps.N(), false);
    for (Uint iedge=0; iedge < edges.size(); ++iedge){
      if (edge_isactive[iedge]){
        for (Uint k=0; k < 2; ++k)
          is_used[edges[iedge].first[k]] = true;
      }
    }
    for ( auto & inode : nodes_inlet )
      is_used[inode] = true;
    for (Uint inode=0; inode < ps.N(); ++inode){
      if (!is_used[inode])
        node_isactive[inode] = false;
    }
  }

  remove_faces(faces, face_isactive);
  remove_edges(faces, edges, edge_isactive, edges_inlet);
  remove_nodes(edges, nodes_inlet, ps, node_isactive);

  if (has_faces)
    compute_edge2faces(edge2faces, faces, edges);
  compute_node2edges(node2edges, edges, ps.N());
}

void check_geometry(FacesType& faces, EdgesType& edges, std::vector<bool>& face_isactive, std::vector<bool>& edge_isactive){
  std::vector<std::vector<Uint>> e2f(edges.size());

  for ( Uint iface=0; iface < faces.size(); ++iface ){
    //std::cout << iface << " ";
    if (face_isactive[iface]){
      for ( auto & iedge : faces[iface].first ){
        //std::cout << iedge << " ";
        e2f[iedge].push_back(iface);
      }
    }
    //std::cout << std::endl;
  }

  for ( Uint iedge=0; iedge < e2f.size(); ++iedge ){
    //std::cout << iedge << " " << e2f[iedge].size() << std::endl;
    if (e2f[iedge].size() > 2) {
      std::cout << "Non-manifold!" << std::endl;
      break;
    }
    if (!edge_isactive[iedge] && e2f[iedge].size() > 0){
      std::cout << "Active face points to inactive edge!" << std::endl;
    }
  }
}

Uint sheet_coarsening(FacesType &faces,
                      EdgesType &edges,
                      Edge2FacesType &edge2faces,
                      Node2EdgesType &node2edges,
                      EdgesListType& edges_inlet,
                      NodesListType& nodes_inlet,
                      ParticleSet& ps,
                      const double ds_min,
                      const double curv_refine_factor){

  std::vector<bool> face_isactive(faces.size(), true);
  std::vector<bool> edge_isactive(edges.size(), true);
  std::vector<bool> node_isactive(ps.N(), true);

  // Do not allow inlet edges or adjacent edges to be collapsed
  std::vector<bool> edge_allow_collapse(edges.size(), true);
  for ( auto & inode : nodes_inlet ){
    for ( auto & jedge : node2edges[inode] ){
      edge_allow_collapse[jedge] = false;
    }
  }
  // This one is for debugging topology:
  // check_geometry(faces, edges, face_isactive, edge_isactive);

  bool changed;
  Uint n_coll = 0;
  Uint iedge;
  CollapseBuffers buf;

  // One linear sweep per pass, as strip_coarsening does. A collapse marks
  // entities inactive rather than erasing them, so every index stays valid
  // and the sweep carries on past it; an edge it shortened earlier in the
  // sweep is picked up by the next pass. Restarting the sweep from zero after
  // each collapse instead cost a full pass per collapse -- 3e10 distance
  // evaluations for 1e5 collapses, 70% of the run.
  do {
    changed = false;
    for (iedge = 0; iedge < edges.size(); ++iedge){
      if (!edge_isactive[iedge] || !edge_allow_collapse[iedge])
        continue;
      Uint inode = edges[iedge].first[0];
      Uint jnode = edges[iedge].first[1];
      // double ds0 = edges[iedge].second;
      double ds = ps.dist(inode, jnode);
      //double kappa = sqrt(abs(H_rw[inode]*H_rw[jnode]));
      //double ds_min_loc = ds_min/(1.0 + curv_refine_factor*kappa);
      double ds_min_loc = ds_min;
      if (ds < ds_min_loc){
        bool coll = collapse_edge(iedge, faces, edges,
                                  edge2faces,
                                  node2edges,
                                  face_isactive,
                                  edge_isactive,
                                  node_isactive,
                                  ps, buf);

        // Debugging topology:
        // check_geometry(faces, edges, face_isactive, edge_isactive);
        if (coll){
          changed = true;
          ++n_coll;
        }
      }
    }
  } while(changed);

  remove_inactive(faces, edges, edge2faces, node2edges,
                  edges_inlet, nodes_inlet,
                  face_isactive, edge_isactive, node_isactive, ps);
  return n_coll;
}

// Adjacent edges normally agree on tau; refuse to merge the rare exception
static bool tau_compatible(const double a, const double b){
  if (a == b) return true;                 // including a run not tracking tau
  const double lo = std::min(a, b), hi = std::max(a, b);
  return lo > 0. && hi <= 3.*lo;
}

Uint strip_coarsening(FacesType &faces,
                      EdgesType &edges,
                      Edge2FacesType &edge2faces,
                      Node2EdgesType &node2edges,
                      EdgesListType& edges_inlet,
                      NodesListType& nodes_inlet,
                      ParticleSet& ps,
                      const double ds_min,
                      const double curv_refine_factor){
  bool changed;
  Uint n_coll = 0;
  Uint iedge;

  std::vector<bool> face_isactive(faces.size(), true);   // a strip has none
  std::vector<bool> edge_isactive(edges.size(), true);
  std::vector<bool> node_isactive(ps.N(), true);

  std::vector<bool> edge_allow_collapse(edges.size(), true);
  for ( auto & inode : nodes_inlet ){
    for ( auto & jedge : node2edges[inode] ){
      edge_allow_collapse[jedge] = false;
    }
  }

  do {
    changed = false;
    iedge = 0;
    while (iedge < edges.size()){
      if (edge_isactive[iedge] && edge_allow_collapse[iedge]){
        Uint inode = edges[iedge].first[0];
        Uint jnode = edges[iedge].first[1];
        double ds0 = edges[iedge].second;
        double ds = ps.dist(inode, jnode);
        // double kappa = sqrt(abs(H_rw[inode]*H_rw[jnode]));
        // double ds_min_loc = ds_min/(1.0 + curv_refine_factor*kappa);
        double ds_min_loc = ds_min;
        Uint new_inode = std::min(inode, jnode);
        Uint old_inode = std::max(inode, jnode);
        // Collapse interior edges only
        const bool has_j = node2edges[old_inode].size() > 1;
        const bool has_k = node2edges[new_inode].size() > 1;

        if (ds < ds_min_loc && edge_isactive[iedge] && has_j && has_k){
          std::vector<Uint> jedges(node2edges[old_inode].begin(),
                              node2edges[old_inode].end());
          std::vector<Uint> kedges(node2edges[new_inode].begin(),
                              node2edges[new_inode].end());
          Uint jedge = get_other(jedges[0], jedges[1], iedge);
          Uint kedge = get_other(kedges[0], kedges[1], iedge);

          // before anything is rewired, so there is nothing to undo
          if (!tau_compatible(edges[jedge].tau, edges[iedge].tau) ||
              !tau_compatible(edges[kedge].tau, edges[iedge].tau)){
            ++iedge;
            continue;
          }

          edge_isactive[iedge] = false;
          node_isactive[old_inode] = false;

          std::replace(edges[jedge].first.begin(), edges[jedge].first.end(),
                  old_inode, new_inode);

          for (const Uint nedge : {jedge, kedge}){
            if (edges[nedge].tau == edges[iedge].tau) continue;
            const double wn = edges[nedge].second, wi = ds0/2;
            // update tau to preserve concentration variance
            const double inv = (wn/sqrt(edges[nedge].tau)
                                + wi/sqrt(edges[iedge].tau)) / (wn + wi);
            edges[nedge].tau = 1./(inv*inv);
          }
          // each neighbour takes half the reference length
          edges[jedge].second += ds0/2;
          edges[kedge].second += ds0/2;

          // collapse_nodes decides where to place the surviving node
          ps.collapse_nodes(inode, jnode, node2edges);

          // ds0 and the geometry have changed, so the elongation must be recomputed
          for (const Uint nedge : {jedge, kedge})
            edges[nedge].rho_prev = ps.dist(edges[nedge].first[0],
                                            edges[nedge].first[1])
                                    / edges[nedge].second;

          node2edges[old_inode].clear();
          std::replace(node2edges[new_inode].begin(),
                  node2edges[new_inode].end(), iedge, jedge);

          changed = true;
          ++n_coll;
        }
      }
      ++iedge;
    }
  } while(changed);

  remove_inactive(faces, edges, edge2faces, node2edges,
                  edges_inlet, nodes_inlet,
                  face_isactive, edge_isactive, node_isactive, ps);
  return n_coll;
}

Uint coarsening(FacesType &faces,
                EdgesType &edges,
                Edge2FacesType &edge2faces,
                Node2EdgesType &node2edges,
                EdgesListType& edges_inlet,
                NodesListType& nodes_inlet,
                ParticleSet& ps,
                const double ds_min,
                const double curv_refine_factor){
  if (faces.size() > 0){ // TODO: better requirement for injection?
    return sheet_coarsening(faces, edges,
                            edge2faces, node2edges,
                            edges_inlet, nodes_inlet,
                            ps, ds_min, curv_refine_factor);
  }
  else if (edges_inlet.size() == 0){ // Omits first steps of injection which might have no faces
    return strip_coarsening(faces, edges,
                            edge2faces, node2edges,
                            edges_inlet, nodes_inlet,
                            ps, ds_min, curv_refine_factor);
  }
  return 0;
}

bool strip_filtering(FacesType &faces,
                     EdgesType &edges,
                     Edge2FacesType &edge2faces,
                     Node2EdgesType &node2edges,
                     EdgesListType &edges_inlet,
                     NodesListType &nodes_inlet,
                     ParticleSet& ps,
                     const Uint filter_target){
  // std::cout << "Edges size: " << edges.size() << std::endl;
  // std::cout << "Filter target: " << filter_target << std::endl;

  if (edges.size() <= filter_target)
    return false;
  while (edges.size() > filter_target){
    int index = rand() % edges.size();
    edges.erase(edges.begin() + index);
  }
  // the edges are already gone; this collects the nodes they leave behind.
  // Not compatible with injection: the schema refuses inject beside filter
  std::vector<bool> face_isactive(faces.size(), true);   // a strip has none
  std::vector<bool> edge_isactive(edges.size(), true);
  std::vector<bool> node_isactive(ps.N(), true);
  remove_inactive(faces, edges, edge2faces, node2edges,
                  edges_inlet, nodes_inlet,
                  face_isactive, edge_isactive, node_isactive, ps);
  return true;
}

bool sheet_filtering(FacesType &faces,
                     EdgesType &edges,
                     Edge2FacesType &edge2faces,
                     Node2EdgesType &node2edges,
                     EdgesListType &edges_inlet,
                     NodesListType &nodes_inlet,
                     ParticleSet& ps,
                     const Uint filter_target){
  std::cout << "SHEET FILTERING NOT TESTED" << std::endl;
  exit(1);
  if (faces.size() <= filter_target)
    return false;
  std::vector<Uint> ids(faces.size());
  iota(ids.begin(), ids.end(), 0);

  // FIXME: Should be made thread safe
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(ids.begin(), ids.end(), g);

  std::vector<bool> face_isactive(faces.size(), false);
  for (Uint iface=0; iface < filter_target; ++iface){
    face_isactive[iface] = true;
  }

  std::vector<bool> edge_isactive(edges.size(), true);
  std::vector<bool> node_isactive(ps.N(), true);
  remove_inactive(faces, edges, edge2faces, node2edges,
                  edges_inlet, nodes_inlet,
                  face_isactive, edge_isactive, node_isactive, ps);
  return true;
}

bool filtering(FacesType &faces,
               EdgesType &edges,
               Edge2FacesType &edge2faces,
               Node2EdgesType &node2edges,
               EdgesListType &edges_inlet,
               NodesListType &nodes_inlet,
               ParticleSet& ps,
               const Uint filter_target){
  if (faces.size() > 0){
    return sheet_filtering(faces, edges, edge2faces, node2edges,
                           edges_inlet, nodes_inlet, ps, filter_target);
  }
  else{
    return strip_filtering(faces, edges, edge2faces, node2edges,
                           edges_inlet, nodes_inlet, ps, filter_target);
  }
}

bool resizing(EdgesType &edges,
              Node2EdgesType &node2edges,
              ParticleSet& ps,
              const double ds){
  bool resized = false;
  for ( auto & edge : edges ){
    Uint inode = edge.first[0];
    Uint jnode = edge.first[1];

    Vector3d xi = ps.x(inode);
    Vector3d xj = ps.x(jnode);
    Vector3d dx = xj - xi;  // such that xj = xi + dx

    double rescale_factor = ds / dx.norm();
    if (rescale_factor < 1.0){
      resized = true;
      // reference and geometry scale together, so ds/ds0 is unchanged
      edge.second *= rescale_factor;
      ps.set_x(jnode, xi + dx * rescale_factor); // such that xj = xi + dx
    }
  }
  return resized;
}

std::array<double, 3> mixed_area_contrib(const double ang0,
                                    const double ang1,
                                    const double ang2,
                                    const double s01,
                                    const double s02,
                                    const double s12){
  double a0, a1, a2;
  if (ang0 <= M_PI_2 && ang1 <= M_PI_2 && ang2 <= M_PI_2){
    double da0 = s12 / tan(ang0);
    double da1 = s02 / tan(ang1);
    double da2 = s01 / tan(ang2);
    a0 = (da1 + da2) / 8;
    a1 = (da2 + da0) / 8;
    a2 = (da0 + da1) / 8;
  }
  else {
    double face_area = sqrt(s01*s02)*sin(ang0)/2;
    a0 = face_area/4;
    a1 = face_area/4;
    a2 = face_area/4;
    if (ang0 > M_PI_2){
      a0 *= 2;
    }
    else if (ang1 > M_PI_2){
      a1 *= 2;
    }
    else {
      a2 *= 2;
    }
  }
  return {a0, a1, a2};
}

double get_angle(const Vector3d &a, const Vector3d &b){
  return acos(a.dot(b)/(a.norm()*b.norm()));
}

void compute_interior_prop(InteriorAnglesType &interior_ang,
                           std::vector<double> &mixed_areas,
                           std::vector<Vector3d> &face_normals,
                           const FacesType &faces,
                           const EdgesType &edges,
                           const Edge2FacesType &edge2faces,
                           ParticleSet& ps){
  interior_ang.clear();
  mixed_areas.clear();
  face_normals.clear();
  mixed_areas.assign(ps.N(), 0.);
  // face2nodes.clear();
  std::set<Uint> not_visited;
  for (Uint iface=0; iface < faces.size(); ++iface){
    not_visited.insert(iface);
    std::set<Uint> inodes_set;
    for (Uint i=0; i < 3; ++i){
      for (Uint j=0; j < 2; ++j){
        inodes_set.insert(edges[faces[iface].first[i]].first[j]);
      }
    }
    assert(inodes_set.size()==3);
    std::vector<Uint> inodes(inodes_set.begin(), inodes_set.end());
    std::array<Vector3d, 3> v;
    for (Uint i=0; i<3; ++i){
      Uint inode = inodes[i];
      v[i] = ps.x(inode);
    }
    double ang0 = get_angle(v[1]-v[0], v[2]-v[0]);
    double ang1 = get_angle(v[2]-v[1], v[0]-v[1]);
    double ang2 = get_angle(v[0]-v[2], v[1]-v[2]);
    interior_ang.push_back({{inodes[0], ang0},
                            {inodes[1], ang1},
                            {inodes[2], ang2}});
    double s01 = (v[1]-v[0]).squaredNorm();
    double s02 = (v[2]-v[0]).squaredNorm();
    double s12 = (v[2]-v[1]).squaredNorm();

    std::array<double, 3> a_loc = mixed_area_contrib(ang0, ang1, ang2, s01, s02, s12);
    for (Uint i=0; i<3; ++i){
      mixed_areas[inodes[i]] += a_loc[i];
    }
    //Vector3d n_loc = get_normal(iface, faces, edges, x_rw);
    Vector3d n_loc = ps.facet_normal(iface, faces, edges);
    face_normals.push_back(n_loc);
  }
  std::set<Uint> to_visit;
  while (not_visited.size() > 0){
    to_visit.insert(*not_visited.begin());
    not_visited.erase(not_visited.begin());
    while (to_visit.size() > 0){
      Uint iface = *to_visit.begin();
      to_visit.erase(to_visit.begin());
      not_visited.erase(iface);
      for (Uint i=0; i<3; ++i){
        Uint iedge = faces[iface].first[i];
        for (auto faceit = edge2faces[iedge].begin();
             faceit != edge2faces[iedge].end(); ++faceit){
          Uint jface = *faceit;
          if (iface != jface && contains(not_visited, jface)){
            to_visit.insert(jface);
            if (face_normals[iface].dot(face_normals[jface]) < 0){
              face_normals[jface] *= -1;
            }
          }
        }
      }
    }
  }
}

void compute_sheet_curv(const FacesType &faces,
                        const EdgesType &edges,
                        const Edge2FacesType &edge2faces,
                        const Node2EdgesType &node2edges,
                        ParticleSet& ps,
                        const InteriorAnglesType &interior_ang,
                        const std::vector<double> &mixed_areas,
                        const std::vector<Vector3d> &face_normals
                        ){
  std::vector<double> edge_w(edges.size(), 0.);
  for (Uint iedge=0; iedge < edges.size(); ++iedge){
    for (auto faceit = edge2faces[iedge].begin();
         faceit != edge2faces[iedge].end(); ++faceit){
      std::vector<Uint> other_edges;
      for (Uint i=0; i<3; ++i){
        Uint jedge = faces[*faceit].first[i];
        if (iedge != jedge)
          other_edges.push_back(jedge);
      }
      // assert(other_edges.size()==2);
      // GL: Hack to avoid crashing. Curvature calculations need improvement if they are going to be used!
      if (other_edges.size() == 2){
        Uint inode = get_intersection(edges[other_edges[0]].first,
                                      edges[other_edges[1]].first);
        // assert(contains(interior_ang[*faceit], inode));
        if (contains(interior_ang[*faceit], inode)){
          std::map<Uint, double> angles = interior_ang[*faceit];
          edge_w[iedge] += 1.0/tan(angles[inode]);
        }
      }
    }
  }
  ps.set_normals(interior_ang, face_normals);
  ps.compute_curvature(edges, node2edges, edge_w, mixed_areas);
}

void compute_strip_curv(const EdgesType &edges,
                        const Node2EdgesType &node2edges,
                        ParticleSet& ps){
  ps.compute_strip_curvature(edges, node2edges);
}

void compute_mean_curv(const FacesType &faces,
                       const EdgesType &edges,
                       const Edge2FacesType &edge2faces,
                       const Node2EdgesType &node2edges,
                       ParticleSet& ps,
                       const InteriorAnglesType &interior_ang,
                       const std::vector<double> &mixed_areas,
                       const std::vector<Vector3d> &face_normals
                       ){
  if (faces.size() > 1){
    compute_sheet_curv(faces, edges,
                       edge2faces, node2edges, ps,
                       interior_ang,
                       mixed_areas, face_normals);
  }
  else {
    compute_strip_curv(edges, node2edges, ps);
  }
}

// Lay down a generation of the inlet. With inject_edges the generations are
// stitched to each other as well, which raises what they trace by a dimension:
// a cloud draws lines, a curve sweeps a surface.
// The stitching is one quad per inlet edge, split by the diagonal N_j -- O_j+1.
// A node the flow has left on the inlet is reused rather than injected again;
// its quad is then a triangle, closed by the old edge or by the new one.
bool injection(const std::vector<Vector3d> &pos_inj,
               const EdgesType &edges_inj,
               EdgesListType &edges_inlet,
               NodesListType &nodes_inlet,
               EdgesType &edges,
               FacesType &faces,
               Edge2FacesType& edge2faces,
               Node2EdgesType& node2edges,
               ParticleSet &ps,
               const bool inject_edges,
               const bool verbose
               ){
  const Uint n_inj = pos_inj.size();
  if (n_inj == 0 || !ps.has_space(n_inj))
    return true;
  assert(nodes_inlet.size() == n_inj);

  // Only the stitching can go degenerate
  std::vector<bool> reused(n_inj, false);
  if (inject_edges){
    double tol = 0.;
    for ( auto & edge : edges_inj )
      tol += edge.second;
    if (!edges_inj.empty())
      tol *= 1e-10/edges_inj.size();
    for (Uint i=0; i < n_inj; ++i)
      reused[i] = (ps.x(nodes_inlet[i]) - pos_inj[i]).norm() <= tol;
  }
  std::vector<Uint> fresh;
  for (Uint i=0; i < n_inj; ++i){
    if (!reused[i])
      fresh.push_back(i);
  }

  const Uint irw0 = ps.N();
  ps.add(pos_inj, fresh, irw0);
  for (Uint k=0; k < fresh.size(); ++k)
    node2edges.push_back({});
  if (verbose)
    std::cout << "Added " << fresh.size() << " nodes." << std::endl;

  // What carries each template node this generation
  NodesListType node_new(n_inj);
  for (Uint i=0, k=0; i < n_inj; ++i)
    node_new[i] = reused[i] ? nodes_inlet[i] : irw0 + k++;

  {
    const Uint iedge0 = edges.size();
    const Uint iface0 = faces.size();
    auto add_edge = [&](const Uint inode, const Uint jnode, const double ds0){
      const Uint iedge = edges.size();
      edges.push_back({{inode, jnode}, ds0});
      edge2faces.push_back({});
      node2edges[inode].push_back(iedge);
      node2edges[jnode].push_back(iedge);
      return iedge;
    };
    auto add_face = [&](const Uint iedge, const Uint jedge, const Uint kedge,
                        const double dA0){
      const Uint iface = faces.size();
      faces.push_back({{iedge, jedge, kedge}, dA0});
      edge2faces[iedge].push_back(iface);
      edge2faces[jedge].push_back(iface);
      edge2faces[kedge].push_back(iface);
    };
    // Twice the area of a corner
    auto wedge = [&](const Uint i, const Uint j, const Uint k){
      return (ps.x(j) - ps.x(i)).cross(ps.x(k) - ps.x(i)).norm();
    };

    // The new curve. An edge both of whose ends stayed put is the old one.
    EdgesListType edge_new(edges_inj.size());
    for (Uint j=0; j < edges_inj.size(); ++j){
      const Uint a = edges_inj[j].first[0];
      const Uint b = edges_inj[j].first[1];
      edge_new[j] = (reused[a] && reused[b])
                  ? edges_inlet[j]
                  : add_edge(node_new[a], node_new[b], edges_inj[j].second);
    }
    // The rungs, one per node that was actually injected
    std::vector<Uint> rung(n_inj, 0);
    for (Uint i=0; inject_edges && i < n_inj; ++i){
      if (!reused[i])
        rung[i] = add_edge(node_new[i], nodes_inlet[i],
                           ps.dist(node_new[i], nodes_inlet[i]));
    }
    // The quads. A corner the flow has flattened is laid down like any other
    // and culled once it is no longer at the inlet
    for (Uint j=0; inject_edges && j < edges_inj.size(); ++j){
      const Uint a = edges_inj[j].first[0];
      const Uint b = edges_inj[j].first[1];
      if (reused[a] && reused[b])
        continue;                     // no material arrived here at all
      const Uint Oa = nodes_inlet[a], Ob = nodes_inlet[b];
      const Uint Na = node_new[a], Nb = node_new[b];
      const Uint old_edge = edges_inlet[j];
      Uint diag;
      if (reused[a])
        diag = old_edge;              // the diagonal is already there
      else if (reused[b])
        diag = edge_new[j];
      else
        diag = add_edge(Na, Ob, ps.dist(Na, Ob));
      if (!reused[a])
        add_face(diag, old_edge, rung[a], wedge(Na, Oa, Ob)/2);
      if (!reused[b])
        add_face(diag, edge_new[j], rung[b], wedge(Na, Nb, Ob)/2);
    }
    if (verbose){
      std::cout << "Added " << edges.size()-iedge0 << " edges." << std::endl;
      std::cout << "Added " << faces.size()-iface0 << " faces." << std::endl;
    }
    edges_inlet = edge_new;
  }
  nodes_inlet = node_new;
  return true;
}


#endif
