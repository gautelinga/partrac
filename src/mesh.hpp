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
void remove_unused_edges(FacesType&, EdgesType&, EdgesListType&);
void remove_unused_nodes(EdgesType&, NodesListType&, ParticleSet&);

// Definitions
/*void add_particles(std::vector<Vector3d> &pos_init, Interpol *intp,
                   std::vector<Vector3d> &x_rw, std::vector<Vector3d> &u_rw,
                   std::vector<double> &c_rw, std::vector<double> &tau_rw,
                   std::vector<double> &rho_rw, std::vector<double> &p_rw,
                   std::vector<Vector3d> &a_rw,
                   const std::string &restart_folder, const int int_order,
                   const Uint irw0) {
                     // TO BECOME REDUNTANT
  Uint Nrw = pos_init.size();
  for (Uint irw=irw0; irw < irw0+Nrw; ++irw){
    // Assign initial position
    x_rw[irw] = pos_init[irw-irw0]; // could be done more efficiently

    if (restart_folder == ""){
      c_rw[irw] = double(irw)/(Nrw-1);
    }
    tau_rw[irw] = 0.;  // anything else?

    intp->probe(x_rw[irw]);
    u_rw[irw] = intp->get_u();

    rho_rw[irw] = intp->get_rho();
    p_rw[irw] = intp->get_p();

    // Second-order terms
    if (int_order >= 2){
      a_rw[irw] = intp->get_Ju() + intp->get_a();
    }
  }
}*/

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
      for (std::set<Uint>::const_iterator sit=unique_nodes.begin();
           sit != unique_nodes.end(); ++sit){
        faces_arr[iface*faces_dims[1] + count] = *sit;
        ++count;
      }
      dA[iface] = ps.triangle_area(iface, faces, edges);
      dA0[iface] = faces[iface].second;
    }
    H5::DataSet faces_dset = h5f.createDataSet(groupname + "/faces",
                                           H5::PredType::NATIVE_ULONG,
                                           faces_dspace);
    faces_dset.write(faces_arr.data(), H5::PredType::NATIVE_ULONG);

    scalar2hdf5(h5f, groupname + "/dA", dA, faces_dims[0]);
    scalar2hdf5(h5f, groupname + "/dA0", dA0, faces_dims[0]);
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

Uint get_common_entry(EdgesListType& iedges, EdgesListType& jedges){
  EdgesListType out;
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
  exit(0);
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
    exit(0);
  }
  // std::cout << mnedges[0] << " " << mnedges[1] << std::endl;
  return {knode, mnedges[0], mnedges[1]};
}

bool append_new_node(const Uint inode, const Uint jnode,
                     std::vector<Vector3d>& x_rw,
                     std::vector<Vector3d>& u_rw,
                     std::vector<double>& rho_rw, std::vector<double>& p_rw, std::vector<double>& c_rw, std::vector<double>& tau_rw,
                     std::vector<double>& H_rw, std::vector<Vector3d>& n_rw,
                     std::vector<Vector3d>& a_rw,
                     Uint& Nrw, const bool do_output_all,
                     Interpol *intp,
                     const int int_order){
                       // To become obsolete...
  Vector3d x_rw_new = 0.5*(x_rw[inode]+x_rw[jnode]);
  int refinement_insertion_levels = 10;
  intp->probe(x_rw_new);
  if (!intp->inside_domain()){
    Vector3d dx_rw_new = x_rw[inode]-x_rw[jnode];
    double dx0 = dx_rw_new.norm();
    Vector3d n0 = u_rw[inode]+u_rw[jnode];
    n0 /= -n0.norm();
    intp->probe(x_rw_new + dx0*n0);
    //exit(0);
    if (!intp->inside_domain()){
      std::cout << "Insertion failed! Information:" << std::endl;
      std::cout << n0 << std::endl;
      std::cout << dx0 << std::endl;
      std::cout << x_rw_new << std::endl;
      std::cout << x_rw_new + dx0*n0 << std::endl;
      // exit(0);

      return false;
    }
    double ddx = dx0/2;
    double dx1 = dx0;
    for (int i=2; i<(2+refinement_insertion_levels); ++i){
      intp->probe(x_rw_new + ddx*n0);
      if (intp->inside_domain()){
        ddx -= dx0/pow(2, i);
        dx1 = ddx;
      }
      else {
        ddx += dx0/pow(2, i);
      }
    }
    x_rw_new += dx1*n0;
    intp->probe(x_rw_new);
  }

  x_rw[Nrw] = x_rw_new;

  c_rw[Nrw] = 0.5*(c_rw[inode]+c_rw[jnode]);
  tau_rw[Nrw] = 0.5*(tau_rw[inode]+tau_rw[jnode]);
  H_rw[Nrw] = 0.5*(H_rw[inode]+H_rw[jnode]);

  n_rw[Nrw] = 0.5*(n_rw[inode]+n_rw[jnode]);
  n_rw[Nrw] /= n_rw[Nrw].norm();

  // u_rw[Nrw] = intp->get_u();  // For some reason this goes wrong?
  u_rw[Nrw] = 0.5*(u_rw[inode]+u_rw[jnode]);
  if (do_output_all){
    // rho_rw[Nrw] = intp->get_rho();
    // p_rw[Nrw] = intp->get_p();
    rho_rw[Nrw] = 0.5*(rho_rw[inode]+rho_rw[jnode]);
    p_rw[Nrw] = 0.5*(p_rw[inode]+p_rw[jnode]);
  }
  // Second-order terms
  if (int_order >= 2){
    // a_rw[Nrw] = intp->get_Ju() + intp->get_a();
    a_rw[Nrw] = 0.5*(a_rw[inode] + a_rw[jnode]);
  }
  ++Nrw;
  return true;
}

Uint sheet_refinement(FacesType &faces,
                      EdgesType &edges,
                      Edge2FacesType &edge2faces,
                      Node2EdgesType &node2edges,                    
                      EdgesListType &edges_inlet,                
                      ParticleSet& ps,
                      const double ds_max,
                      const double curv_refine_factor,
                      const bool cut_if_stuck){
  bool changed;
  Uint n_add = 0;
  std::set<Uint> edges_to_remove;

  do {
    changed = false;

    std::vector<double> ds_ratio_;
    for (EdgesType::const_iterator edgeit = edges.begin();
         edgeit != edges.end(); ++edgeit){
      Uint inode = edgeit->first[0];
      Uint jnode = edgeit->first[1];
      double ds = ps.dist(inode, jnode);
      // std::cout << inode << " " << jnode << " " << ds << std::endl;
      //double kappa = 0.5*(abs(H_rw[inode]) + abs(H_rw[jnode]));
      //double ds_max_loc = ds_max/(1.0 + curv_refine_factor*kappa);
      double ds_max_loc = ds_max;
      ds_ratio_.push_back(ds/ds_max_loc);
    }
    std::vector<size_t> ids_ = argsort_descending(ds_ratio_);

    // Don't refine inlet edges
    for (EdgesListType::const_iterator nodeit = edges_inlet.begin();
         nodeit != edges_inlet.end(); ++nodeit){
      ds_ratio_[*nodeit] = 0.0;
    }

    for (std::vector<size_t>::iterator itedge = ids_.begin();
         itedge != ids_.end(); ++itedge){

      Uint iedge = *itedge;
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
      bool added = ps.insert_node_between(inode, jnode);
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

        for (FacesListType::iterator itface = edge2faces[iedge].begin();
             itface != edge2faces[iedge].end(); ++itface){
          Uint jedge = faces[*itface].first[0];
          Uint kedge = faces[*itface].first[1];
          Uint ledge = faces[*itface].first[2];
          double dA0 = faces[*itface].second;
          Uint new_jedge = edges.size();
          std::array<Uint, 3> close_entities = get_close_entities(iedge, jedge, kedge, ledge, edges);
          Uint knode = close_entities[0];
          Uint medge = close_entities[1];
          Uint nedge = close_entities[2];
          edges.push_back({{new_inode, knode}, ds0/2});  // ds0/2 - or what else?
          edge2faces.push_back({});

          Uint new_iface = faces.size();
          faces[*itface] = {{iedge, new_jedge, medge}, dA0/2};
          faces.push_back({{new_iedge, nedge, new_jedge}, dA0/2});

          edge2faces[nedge].remove(*itface);
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
        std::cout << "Here we should remove this edge." << std::endl;
        //exit(0);
        edges_to_remove.insert(iedge);
      }
    }
  } while (changed);
  if (edges_to_remove.size() > 0){
    if (!cut_if_stuck){
      std::cout << "Edge is stuck! Turn on 'cut_if_stuck' to continue in such cases." << std::endl;
      exit(0);
    }
    std::vector<bool> edge_isactive(edges.size(), true);
    std::vector<bool> face_isactive(faces.size(), true);
    for (std::set<Uint>::const_iterator sit = edges_to_remove.begin();
         sit != edges_to_remove.end(); ++sit){
      edge_isactive[*sit] = false;
      for (FacesListType::const_iterator jfaceit=edge2faces[*sit].begin();
           jfaceit != edge2faces[*sit].end(); ++jfaceit){
          face_isactive[*jfaceit] = false;
      }
    }
    NodesListType nodes_inlet_dummy;
    EdgesListType edges_inlet_dummy;

    remove_faces(faces, face_isactive);
    remove_edges(faces, edges, edge_isactive, edges_inlet_dummy);
    remove_unused_nodes(edges, nodes_inlet_dummy, ps);
    compute_edge2faces(edge2faces, faces, edges);
    compute_node2edges(node2edges, edges, ps.N());
  }
  return n_add;
}

Uint strip_refinement(EdgesType &edges,
                      Node2EdgesType &node2edges,
                      ParticleSet& ps,
                      const double ds_max,
                      const double curv_refine_factor,
                      const bool cut_if_stuck){
  Uint n_add = 0;
  Uint iedge = 0;
  std::set<Uint> edges_to_remove;
  while (iedge < edges.size()){
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

        node2edges.push_back({iedge, new_iedge});

        std::replace(node2edges[jnode].begin(), node2edges[jnode].end(), iedge, new_iedge);

        //compute_node2edges(node2edges, edges, new_inode);

        ++n_add;
      }
      else {
        // exit(0);
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
      exit(0);
    }
    std::vector<bool> edge_isactive(edges.size(), true);
    for (std::set<Uint>::const_iterator sit = edges_to_remove.begin();
         sit != edges_to_remove.end(); ++sit){
      edge_isactive[*sit] = false;
    }
    FacesType faces_dummy;
    NodesListType nodes_inlet_dummy;
    EdgesListType edges_inlet_dummy;
    remove_edges(faces_dummy, edges, edge_isactive, edges_inlet_dummy);
    remove_unused_nodes(edges, nodes_inlet_dummy, ps);
    compute_node2edges(node2edges, edges, ps.N());
  }
  return n_add;
}

Uint refinement(FacesType &faces,
                EdgesType &edges,
                Edge2FacesType &edge2faces,
                Node2EdgesType &node2edges,
                EdgesListType &edges_inlet,
                ParticleSet& ps, const double ds_max,
                const double curv_refine_factor,
                const bool cut_if_stuck){
  Uint n_add = 0;
  if (faces.size() > 0){
    n_add = sheet_refinement(faces, edges, edge2faces, node2edges, edges_inlet,
                             ps, ds_max, curv_refine_factor, cut_if_stuck);
  }
  else {
    n_add = strip_refinement(edges,
                             node2edges,
                             ps, ds_max,
                             curv_refine_factor,
                             cut_if_stuck);
  }
  return n_add;
}

void compute_edge2faces(Edge2FacesType &edge2faces,
                        const FacesType &faces,
                        const EdgesType &edges){
  edge2faces.clear();
  for (EdgesType::const_iterator edgeit = edges.begin();
       edgeit != edges.end(); ++edgeit){
    edge2faces.push_back({});
  }
  for (Uint iface=0; iface < faces.size(); ++iface){
    for (Uint i=0; i < 3; ++i){
      edge2faces[faces[iface].first[i]].push_back(iface);
    }
  }
}

void compute_node2edges(Node2EdgesType &node2edges,
                        const EdgesType &edges,
                        const Uint Nrw){
  node2edges.clear();
  for (Uint irw=0; irw<Nrw; ++irw){
    node2edges.push_back({});
  }
  for (Uint iedge=0; iedge < edges.size(); ++iedge){
    for (Uint i=0; i < 2; ++i){
      node2edges[edges[iedge].first[i]].push_back(iedge);
    }
  }
}

void print(const std::vector<Uint> vec){
  for (std::vector<Uint>::const_iterator vit = vec.begin();
       vit != vec.end(); ++vit){
    std::cout << *vit << " ";
  }
  std::cout << std::endl;
}

void print(const std::set<Uint> vec){
  for (std::set<Uint>::const_iterator vit = vec.begin();
       vit != vec.end(); ++vit){
    std::cout << *vit << " ";
  }
  std::cout << std::endl;
}

void print(const std::list<Uint> vec){
  for (std::list<Uint>::const_iterator vit = vec.begin();
       vit != vec.end(); ++vit){
    std::cout << *vit << " ";
  }
  std::cout << std::endl;
}

void get_conodes(std::set<Uint> &inodes,
                 std::map<Uint, Uint> &icoedges,
                 const Uint inode,
                 const Node2EdgesType &node2edges,
                 const EdgesType &edges){
  //std::cout << "Getting conodes of " << inode << std::endl;
  inodes.clear();
  icoedges.clear();

  //print(inodes);

  for (EdgesListType::const_iterator edgeit=node2edges[inode].begin();
       edgeit != node2edges[inode].end(); ++edgeit){
    for (Uint j=0; j < 2; ++j){
      Uint jnode = edges[*edgeit].first[j];
      if (jnode != inode){
        inodes.insert(jnode);
        icoedges[jnode] = *edgeit;
      }
    }
  }
  //print(inodes);
}

std::set<Uint> get_cofaces(const Uint iedge,
                      const Edge2FacesType &edge2faces,
                      const FacesType &faces){
  FacesListType ifaces = edge2faces[iedge];
  std::set<Uint> jfaces;
  for (FacesListType::const_iterator ifaceit=edge2faces[iedge].begin();
       ifaceit != edge2faces[iedge].end(); ++ifaceit){
    for (Uint k=0; k<3; ++k){
      Uint kedge = faces[*ifaceit].first[k];
      if (kedge != iedge){
        for (FacesListType::const_iterator kfaceit=edge2faces[kedge].begin();
             kfaceit != edge2faces[kedge].end(); ++kfaceit){
          if (*kfaceit != *ifaceit){
            jfaces.insert(*kfaceit);
          }
        }
      }
    }
  }
  return jfaces;
}

bool is_border_node(const Uint inode,
                    const Node2EdgesType &node2edges,
                    const Edge2FacesType &edge2faces){
  for (EdgesListType::const_iterator edgeit=node2edges[inode].begin();
       edgeit != node2edges[inode].end(); ++edgeit){
    if (edge2faces[*edgeit].size() == 1){
      return true;
    }
  }
  return false;
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
  else if (both_are_border || none_are_border){
    x = 0.5*(ps.x(inode) + ps.x(jnode));
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

bool normals_are_ok(const Uint iedge,
                    const Vector3d &x,
                    std::vector<Vector3d>& x_rw,
                    const std::set<Uint> &jfaces,
                    const FacesType &faces,
                    const EdgesType &edges){
  Uint inode = edges[iedge].first[0];
  Uint jnode = edges[iedge].first[1];

  for (std::set<Uint>::const_iterator jfaceit=jfaces.begin();
       jfaceit != jfaces.end(); ++jfaceit){
    Uint jedge = faces[*jfaceit].first[0];
    Uint kedge = faces[*jfaceit].first[1];

    Vector3d n_before = get_normal(jedge, kedge, edges, x_rw, {}, x);
    Vector3d n_after = get_normal(jedge, kedge, edges, x_rw,
                                  {inode, jnode}, x);
    if (n_before.dot(n_after) <= 0.0)
      return false;
  }
  return true;
}

std::vector<Uint> get_incident_faces(const Uint iedge,
                                const EdgesType &edges,
                                const Edge2FacesType &edge2faces,
                                const Node2EdgesType &node2edges){
  Uint inode = edges[iedge].first[0];
  Uint jnode = edges[iedge].first[1];
  std::set<Uint> kfaces_set;
  for (EdgesListType::const_iterator edgeit=node2edges[inode].begin();
       edgeit != node2edges[inode].end(); ++edgeit){
    for (FacesListType::const_iterator faceit=edge2faces[*edgeit].begin();
         faceit != edge2faces[*edgeit].end(); ++faceit){
      kfaces_set.insert(*faceit);
    }
  }
  for (EdgesListType::const_iterator edgeit=node2edges[jnode].begin();
       edgeit != node2edges[jnode].end(); ++edgeit){
    for (FacesListType::const_iterator faceit=edge2faces[*edgeit].begin();
         faceit != edge2faces[*edgeit].end(); ++faceit){
      kfaces_set.insert(*faceit);
    }
  }
  for (FacesListType::const_iterator faceit=edge2faces[iedge].begin();
       faceit != edge2faces[iedge].end(); ++faceit){
    kfaces_set.erase(*faceit);
  }
  std::vector<Uint> kfaces(kfaces_set.begin(), kfaces_set.end());
  return kfaces;
}

bool collapse_edge(const Uint iedge,
                   FacesType &faces,
                   EdgesType &edges,
                   Edge2FacesType &edge2faces,
                   Node2EdgesType &node2edges,
                   std::vector<bool> &face_isactive,
                   std::vector<bool> &edge_isactive,
                   std::vector<bool> &node_isactive,
                   ParticleSet& ps){

  assert(edge_isactive[iedge]);

  Uint inode = edges[iedge].first[0];
  Uint jnode = edges[iedge].first[1];
  // double ds0 = edges[iedge].second;
  std::set<Uint> inodes;
  std::set<Uint> jnodes;
  std::map<Uint, Uint> icoedges;
  std::map<Uint, Uint> jcoedges;
  get_conodes(inodes, icoedges, inode, node2edges, edges);
  get_conodes(jnodes, jcoedges, jnode, node2edges, edges);
  std::vector<Uint> joint_nodes;
  set_intersection(inodes.begin(), inodes.end(),
                   jnodes.begin(), jnodes.end(),
                   back_inserter(joint_nodes));

  //std::cout << "iedge=" << iedge << std::endl;
  //print(inodes);
  //print(jnodes);
  //print(joint_nodes);

  // Check if it has too many/too few common joint nodes.
  if (joint_nodes.size() > 2){
    // Wrong number of joint nodes
    // std::cout << "Wrong number of joint nodes." << std::endl;
    return false;
  }
  Vector3d x;
  bool new_node_is_good = get_new_pos(x, ps, iedge, edges, edge2faces, node2edges);
  //bool new_node_is_good = ps.get_new_pos ...?
  if (!new_node_is_good){
    // Corner edge
    // std::cout << "Corner node!" << std::endl;
    return false;
  }
  //std::cout << "New pos: " << x << " " << y << " " << z << std::endl;

  std::set<Uint> jfaces = get_cofaces(iedge, edge2faces, faces);

  //print(jfaces);

  double dA_res = 0.;
  double dA0_res = 0.;
  for (FacesListType::const_iterator faceit=edge2faces[iedge].begin();
       faceit != edge2faces[iedge].end(); ++faceit){
    //dA_res += area(*faceit, x_rw, faces, edges);
    dA_res += ps.triangle_area(*faceit, faces, edges);
    dA0_res += faces[*faceit].second;
    //faces[*faceit].second = 0.;
  }

  std::vector<Uint> kfaces = get_incident_faces(iedge, edges, edge2faces, node2edges);
  std::vector<double> dAs_old = ps.triangle_areas(kfaces, faces, edges);

  // assert(jfaces.size()==4 || jfaces.size()==2);

  /*
  if (!normals_are_ok(iedge, x, x_rw,
                      jfaces, faces, edges)){
    // Normals are not ok
    // std::cout << "Normals are not OK!" << std::endl;
    return false;
  }
  */
  /*
  intp->probe(x);  //

  Uint irws[2] = {inode, jnode};
  for (Uint i=0; i<2; ++i){
    Uint irw = irws[i];
    x_rw[irw] = x;
    u_rw[irw] = intp->get_u();

    if (do_output_all){
      rho_rw[irw] = intp->get_rho();
      p_rw[irw] = intp->get_p();
    }
    // Second-order terms
    if (int_order >= 2){
      a_rw[irw] = intp->get_Ju() + intp->get_a();
    }
    c_rw[irw] = 0.5*(c_rw[inode]+c_rw[jnode]);
    tau_rw[irw] = 0.5*(tau_rw[inode]+tau_rw[jnode]);
    H_rw[irw] = 0.5*(H_rw[inode]+H_rw[jnode]);
    n_rw[irw] = 0.5*(n_rw[inode]+n_rw[jnode]);
  }
  */
  ps.replace_nodes(x, inode, jnode);

  for (EdgesListType::const_iterator edgeit=node2edges[std::max(inode, jnode)].begin();
       edgeit != node2edges[std::max(inode, jnode)].end(); ++edgeit){
    for (Uint j=0; j<2; ++j){
      if (edges[*edgeit].first[j] == std::max(inode, jnode)){
        edges[*edgeit].first[j] = std::min(inode, jnode);
      }
    }
  }

  std::vector<Uint> iedges_vec;
  set_symmetric_difference(node2edges[inode].begin(),
                           node2edges[inode].end(),
                           node2edges[jnode].begin(),
                           node2edges[jnode].end(),
                           back_inserter(iedges_vec));
  std::set<Uint> iedges_set(iedges_vec.begin(), iedges_vec.end());
  //print(node2edges[inode]);
  //print(node2edges[jnode]);
  //print(iedges_set);

  std::map<Uint, Uint> replace_edges;
  for (std::vector<Uint>::const_iterator joint_node_it=joint_nodes.begin();
       joint_node_it != joint_nodes.end(); ++joint_node_it){
    Uint kedge_min = std::min(icoedges[*joint_node_it], jcoedges[*joint_node_it]);
    Uint kedge_max = std::max(icoedges[*joint_node_it], jcoedges[*joint_node_it]);
    replace_edges[kedge_max] = kedge_min;
    assert(edge_isactive[kedge_max]);
    edge_isactive[kedge_max] = false;
    if (contains(iedges_set, kedge_max))
      iedges_set.erase(kedge_max);
    //std::cout << "kedge_max=" << kedge_max << ", kedge_min=" << kedge_min << std::endl;
  }
  //print(iedges_set);

  assert(edge_isactive[iedge]);
  edge_isactive[iedge] = false;

  for (std::set<Uint>::const_iterator jfaceit=jfaces.begin();
       jfaceit != jfaces.end(); ++jfaceit){
    //std::cout << "faces[" << *jfaceit << "].first="
    //     << faces[*jfaceit].first[0] << " "
    //     << faces[*jfaceit].first[1] << " "
    //     << faces[*jfaceit].first[2] << std::endl;
    for (Uint j=0; j<3; ++j){
      if (contains(replace_edges, faces[*jfaceit].first[j])){
        faces[*jfaceit].first[j] = replace_edges[faces[*jfaceit].first[j]];
      }
    }
    //std::cout << "faces[" << *jfaceit << "].first="
    //     << faces[*jfaceit].first[0] << " "
    //     << faces[*jfaceit].first[1] << " "
    //     << faces[*jfaceit].first[2] << std::endl;
  }
  std::vector<double> dAs_new = ps.triangle_areas(kfaces, faces, edges);
  double wsum = 0.;
  std::vector<double> w_;
  for (Uint k=0; k < kfaces.size(); ++k){
    //Uint kface = kfaces[k];
    // std::cout << faces[kface].second << "+=" << dA0_res << "/" << dA_res << "*(" << dAs_new[k] << "-" << dAs_old[k] << ")" << std::endl;
    // std::cout << dA0_res << " " << dA_res << std::endl;
    // if (dA_res > 0.0){
    //     faces[kface].second += dA0_res/dA_res*(dAs_new[k]-dAs_old[k]);
    // faces[kface].second += dA0_res/kfaces.size();
    double w = dAs_new[k]-dAs_old[k]; // May give negative mass
    w = std::max(w, 0.0);
    //faces[kface].second += dA0_res*w;
    w_.push_back(w);
    wsum += w;
  }
  if (wsum <= 0.0){
    wsum = 0.0;
    for (Uint k=0; k < kfaces.size(); ++k){
      double w = dAs_new[k];
      w_[k] = w;
      wsum += w;
    }
  }
  for (Uint k=0; k < kfaces.size(); ++k){
    Uint kface = kfaces[k];
    faces[kface].second += dA0_res*w_[k]/wsum;
  }

  for (FacesListType::const_iterator jfaceit=edge2faces[iedge].begin();
       jfaceit != edge2faces[iedge].end(); ++jfaceit){
    assert(face_isactive[*jfaceit]);
    face_isactive[*jfaceit] = false;
  }

  //print(node2edges[std::min(inode, jnode)]);
  node2edges[std::min(inode, jnode)].assign(iedges_set.begin(), iedges_set.end());
  //print(node2edges[std::min(inode, jnode)]);
  node2edges[std::max(inode, jnode)].clear();
  assert(node_isactive[std::max(inode, jnode)]);
  node_isactive[std::max(inode, jnode)] = false;
  for (std::vector<Uint>::const_iterator joint_node_it=joint_nodes.begin();
       joint_node_it != joint_nodes.end(); ++joint_node_it){
    std::set<Uint> tmp_set(node2edges[*joint_node_it].begin(),
                      node2edges[*joint_node_it].end());
    //print(node2edges[*joint_node_it]);
    for (std::map<Uint, Uint>::const_iterator mit=replace_edges.begin();
         mit != replace_edges.end(); ++mit){
      tmp_set.erase(mit->first);
    }
    node2edges[*joint_node_it].assign(tmp_set.begin(), tmp_set.end());
    //print(node2edges[*joint_node_it]);
  }

  for (std::map<Uint, Uint>::const_iterator mit=replace_edges.begin();
       mit != replace_edges.end(); ++mit){
    Uint kedge_max = mit->first;
    Uint kedge_min = mit->second;
    std::vector<Uint> ifaces;
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
  for (int i=int(face_isactive.size())-1; i >= 0; --i){
    if (!face_isactive[i])
      faces.erase(faces.begin()+i);
  }
}

void remove_edges(FacesType &faces, EdgesType &edges,
                  const std::vector<bool> &edge_isactive,
                  EdgesListType &edges_inlet){
  assert(edges.size() == edge_isactive.size());
  std::vector<Uint> used_edges;
  for (Uint i=0; i<edges.size(); ++i){
    used_edges.push_back(i);
  }

  for (int i=edge_isactive.size()-1; i >= 0; --i){
    if (!edge_isactive[i]){
      edges.erase(edges.begin()+i);
      used_edges.erase(used_edges.begin()+i);
    }
  }
  std::map<Uint, Uint> replace_edges;
  for (Uint i=0; i<used_edges.size(); ++i){
    replace_edges[used_edges[i]] = i;
  }
  for (Uint iface=0; iface < faces.size(); ++iface){
    for (Uint j=0; j<3; ++j){
      Uint iedge = faces[iface].first[j];
      if (contains(replace_edges, iedge)){
        faces[iface].first[j] = replace_edges[iedge];
      }
    }
  }
  //for (Uint i=0; i < edges_inlet.size(); ++i){
  for (EdgesListType::iterator edgeit=edges_inlet.begin();
       edgeit != edges_inlet.end(); ++edgeit){
    Uint iedge = *edgeit;
    if (contains(replace_edges, iedge)){
      *edgeit = replace_edges[iedge];
    }
  }
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
  std::map<Uint, Uint> replace_nodes;
  for (Uint i=0; i<used_nodes.size(); ++i){
    replace_nodes[used_nodes[i]] = i;
    Uint j = used_nodes[i];
    ps.copy_node(i, j);
  }
  ps.set_N(used_nodes.size());

  for (Uint iedge=0; iedge < edges.size(); ++iedge){
    for (Uint j=0; j<2; ++j){
      Uint inode = edges[iedge].first[j];
      if (contains(replace_nodes, inode)){
        edges[iedge].first[j] = replace_nodes[inode];
      }
    }
  }
  //for (Uint i=0; i<nodes_inlet.size(); ++i){
  for (NodesListType::iterator nodeit=nodes_inlet.begin();
       nodeit != nodes_inlet.end(); ++nodeit){
    //Uint inode = nodes_inlet[i];
    Uint inode = *nodeit;
    if (contains(replace_nodes, inode)){
      //nodes_inlet[inode] = replace_nodes[inode];
      *nodeit = replace_nodes[inode];
    }
  }
}

void remove_unused_edges(FacesType &faces, EdgesType &edges, EdgesListType &edges_inlet){
  std::vector<bool> edge_isactive(edges.size(), false);
  for (Uint i=0; i < faces.size(); ++i){
    for (Uint k=0; k < 3; ++k){
      edge_isactive[faces[i].first[k]] = true;
    }
  }
  remove_edges(faces, edges, edge_isactive, edges_inlet);
}

void remove_unused_nodes(EdgesType &edges,
                         NodesListType &nodes_inlet,
                         ParticleSet& ps){
  std::vector<bool> node_isactive(ps.N(), false);
  for (Uint i=0; i < edges.size(); ++i){
    for (Uint k=0; k < 2; ++k){
      node_isactive[edges[i].first[k]] = true;
    }
  }
  remove_nodes(edges, nodes_inlet, ps, node_isactive);
}

/*void assert_equal(const Edge2FacesType &a,
                  const Edge2FacesType &b){
  assert(a.size() == b.size());
  for (Uint i=0; i < a.size(); ++i){
    assert(a[i].size() == b[i].size());
  }
  exit(0);
}*/

void remove_nodes_safely(FacesType &faces, EdgesType &edges,
                         Edge2FacesType &edge2faces, Node2EdgesType &node2edges,
                         EdgesListType &edges_inlet, NodesListType &nodes_inlet,
                         const std::vector<bool> &node_isactive,
                         ParticleSet& ps) {
  assert(ps.N() == node_isactive.size());
  std::vector<bool> edge_isactive(edges.size(), true);
  for (Uint i=0; i<ps.N(); ++i){
    if (!node_isactive[i]){
      for (EdgesListType::const_iterator edgeit = node2edges[i].begin();
           edgeit != node2edges[i].end(); ++edgeit){
        edge_isactive[*edgeit] = false;
      }
    }
  }
  if (faces.size() > 0){
    std::cout << "WARNING: remove_nodes_safely is NOT TESTED for sheets!" << std::endl;
    exit(0);
    // Remove node from sheet (cut hole)
    std::vector<bool> face_isactive(faces.size(), true);
    for (Uint iedge=0; iedge<edges.size(); ++iedge){
      if (!edge_isactive[iedge]){
        for (FacesListType::const_iterator faceit = edge2faces[iedge].begin();
             faceit != edge2faces[iedge].end(); ++faceit){
          face_isactive[*faceit] = false;
        }
      }
    }
    remove_faces(faces, face_isactive);
    remove_unused_edges(faces, edges, edges_inlet);
    remove_unused_nodes(edges, nodes_inlet, ps);

    compute_edge2faces(edge2faces, faces, edges);
    compute_node2edges(node2edges, edges, ps.N());
  }
  else {
    // Remove node from strip (cut hole)
    EdgesListType edges_inlet_dummy;
    remove_edges(faces, edges, edge_isactive, edges_inlet_dummy);
    remove_unused_nodes(edges, nodes_inlet, ps);
    compute_node2edges(node2edges, edges, ps.N());
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
  bool changed;
  Uint n_coll = 0;
  Uint iedge;

  std::vector<bool> face_isactive(faces.size(), true);
  std::vector<bool> edge_isactive(edges.size(), true);
  std::vector<bool> node_isactive(ps.N(), true);

  // Do not allow inlet edges or adjacent edges to be collapsed
  std::vector<bool> edge_allow_collapse(edges.size(), true);
  for (NodesListType::const_iterator nodeit=nodes_inlet.begin();
       nodeit != nodes_inlet.end(); ++nodeit){
    for (EdgesListType::const_iterator edgeit = node2edges[*nodeit].begin();
         edgeit != node2edges[*nodeit].end(); ++edgeit){
      edge_allow_collapse[*edgeit] = false;
    }
  }
  //
  do {
    changed = false;
    iedge = 0;
    while (iedge < edges.size()){
      Uint inode = edges[iedge].first[0];
      Uint jnode = edges[iedge].first[1];
      // double ds0 = edges[iedge].second;
      double ds = ps.dist(inode, jnode);
      //double kappa = sqrt(abs(H_rw[inode]*H_rw[jnode]));
      //double ds_min_loc = ds_min/(1.0 + curv_refine_factor*kappa);
      double ds_min_loc = ds_min;
      if (ds < ds_min_loc && edge_isactive[iedge] && edge_allow_collapse[iedge]){
        // std::cout << "TRYING to collapse edge " << iedge
        //      << ". ds = " << ds << " and ds_min = " << ds_min << std::endl;
        bool coll = collapse_edge(iedge, faces, edges,
                                  edge2faces,
                                  node2edges,
                                  face_isactive,
                                  edge_isactive,
                                  node_isactive,
                                  ps);
        // what's wrong?
        //compute_edge2faces(edge2faces, faces, edges);
        //Edge2FacesType edge2faces_alt;
        //compute_edge2faces(edge2faces_alt, faces, edges);
        //assert_equal(edge2faces_alt, edge2faces);
        //compute_node2edges(node2edges, edges, Nrw);
        if (coll){
          changed = true;
          ++n_coll;
          break;
        }
      }
      ++iedge;
    }
  } while(changed);

  remove_faces(faces, face_isactive);
  //remove_edges(faces, edges, edge_isactive);
  //remove_nodes(edges, x_rw, y_rw, z_rw, Nrw, node_isactive);
  remove_unused_edges(faces, edges, edges_inlet);
  remove_unused_nodes(edges, nodes_inlet, ps);

  compute_edge2faces(edge2faces, faces, edges);
  compute_node2edges(node2edges, edges, ps.N());
  return n_coll;
}

Uint strip_coarsening(EdgesType &edges,
                      Node2EdgesType &node2edges,
                      NodesListType &nodes_inlet,
                      ParticleSet& ps,
                      const double ds_min,
                      const double curv_refine_factor){
  bool changed;
  Uint n_coll = 0;
  Uint iedge;

  std::vector<bool> edge_isactive(edges.size(), true);
  std::vector<bool> node_isactive(ps.N(), true);

  std::vector<bool> edge_allow_collapse(edges.size(), true);
  for (NodesListType::const_iterator nodeit=nodes_inlet.begin();
       nodeit != nodes_inlet.end(); ++nodeit){
    for (EdgesListType::const_iterator edgeit=node2edges[*nodeit].begin();
         edgeit != node2edges[*nodeit].end(); ++edgeit){
      edge_allow_collapse[*edgeit] = false;
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
        if (ds < ds_min_loc && edge_isactive[iedge]){
          edge_isactive[iedge] = false;

          Uint new_inode = std::min(inode, jnode);
          Uint old_inode = std::max(inode, jnode);

          node_isactive[old_inode] = false;

          std::vector<Uint> jedges(node2edges[old_inode].begin(),
                              node2edges[old_inode].end());
          Uint jedge = get_other(jedges[0], jedges[1], iedge);
          std::replace(edges[jedge].first.begin(), edges[jedge].first.end(),
                  old_inode, new_inode);
          std::vector<Uint> kedges(node2edges[new_inode].begin(),
                              node2edges[new_inode].end());
          Uint kedge = get_other(kedges[0], kedges[1], iedge);

          edges[jedge].second += ds0/2;
          edges[kedge].second += ds0/2;

          // Do something about tau, rho_prev?

          node2edges[old_inode].clear();
          std::replace(node2edges[new_inode].begin(),
                  node2edges[new_inode].end(), iedge, jedge);

          ps.collapse_nodes(inode, jnode, node2edges);

          changed = true;
          ++n_coll;
          //--iedge;
        }
      }
      ++iedge;
    }
  } while(changed);

  FacesType faces;
  EdgesListType edges_inlet_dummy;
  remove_edges(faces, edges, edge_isactive, edges_inlet_dummy);
  remove_unused_nodes(edges, nodes_inlet, ps);
  compute_node2edges(node2edges, edges, ps.N());
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
    return strip_coarsening(edges, node2edges,
                            nodes_inlet,
                            ps, ds_min, curv_refine_factor);
  }
  return 0;
}

bool strip_filtering(EdgesType &edges,
                     Node2EdgesType &node2edges,
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
  NodesListType nodes_inlet_dummy; // Not compatible with injection (yet)
  remove_unused_nodes(edges, nodes_inlet_dummy, ps);
  compute_node2edges(node2edges, edges, ps.N());

  return true;
}

bool sheet_filtering(FacesType &faces,
                     EdgesType &edges,
                     Edge2FacesType &edge2faces,
                     Node2EdgesType &node2edges,
                     ParticleSet& ps,
                     const Uint filter_target){
  std::cout << "SHEET FILTERING NOT TESTED" << std::endl;
  exit(0);
  if (faces.size() <= filter_target)
    return false;
  std::vector<Uint> ids(faces.size());
  iota(ids.begin(), ids.end(), 0);
  random_shuffle(ids.begin(), ids.end());

  std::vector<bool> face_isactive(faces.size(), false);
  for (Uint iface=0; iface < filter_target; ++iface){
    face_isactive[iface] = true;
  }

  EdgesListType edges_inlet_dummy; // Not compatible with injection (yet)
  NodesListType nodes_inlet_dummy; // Not compatible with injection (yet)

  remove_faces(faces, face_isactive);
  remove_unused_edges(faces, edges, edges_inlet_dummy);
  remove_unused_nodes(edges, nodes_inlet_dummy, ps);

  compute_edge2faces(edge2faces, faces, edges);
  compute_node2edges(node2edges, edges, ps.N());

  return true;
}

bool filtering(FacesType &faces,
               EdgesType &edges,
               Edge2FacesType &edge2faces,
               Node2EdgesType &node2edges,
               ParticleSet& ps,
               const Uint filter_target){
  if (faces.size() > 0){
    return sheet_filtering(faces, edges, edge2faces, node2edges, ps, filter_target);
  }
  else{
    return strip_filtering(edges, node2edges, ps, filter_target);
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
    double ds0 = edge.second;

    Vector3d xi = ps.x(inode);
    Vector3d xj = ps.x(jnode);
    Vector3d dx = xj - xi;  // such that xj = xi + dx

    double rescale_factor = ds / dx.norm();
    if (rescale_factor < 1.0){
      resized = true;
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
        for (FacesListType::const_iterator faceit = edge2faces[iedge].begin();
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
    for (FacesListType::const_iterator faceit = edge2faces[iedge].begin();
         faceit != edge2faces[iedge].end(); ++faceit){
      std::vector<Uint> other_edges;
      for (Uint i=0; i<3; ++i){
        Uint jedge = faces[*faceit].first[i];
        if (iedge != jedge)
          other_edges.push_back(jedge);
      }
      assert(other_edges.size()==2);
      Uint inode = get_intersection(edges[other_edges[0]].first,
                                    edges[other_edges[1]].first);
      assert(contains(interior_ang[*faceit], inode));
      std::map<Uint, double> angles = interior_ang[*faceit];
      edge_w[iedge] += 1.0/tan(angles[inode]);
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
  if (ps.has_space(pos_inj.size())){
    Uint irw0 = ps.N();
    Uint iedge0 = edges.size();
    Uint iface0 = faces.size();
    /*
    add_particles(pos_inj, intp,
                  x_rw, u_rw, c_rw, tau_rw, rho_rw, p_rw, a_rw,
                  U0, prm.restart_folder, prm.int_order, irw0);
    Nrw += pos_inj.size();
    */
    ps.add(pos_inj, irw0);
    for (Uint irw=0; irw < pos_inj.size(); ++irw){
      node2edges.push_back({});
    }
    if (verbose){
      std::cout << "Added " << pos_inj.size() << " nodes." << std::endl;
      /*
      std::cout << "at ..." << std::endl;
      for (Uint irw=irw0; irw<Nrw; ++irw){
        std::cout << x_rw[irw][0] << " " << x_rw[irw][1] << " " << x_rw[irw][2] << std::endl;
      }
      */
    }
    if (inject_edges){
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
        double ds0 = ps.dist(inode, jnode);
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
        double ds0 = ps.dist(inode, jnode);
        edges.push_back({{inode, jnode}, ds0});
        edge2faces.push_back({});
        node2edges[inode].push_back(iedge);
        node2edges[jnode].push_back(iedge);
        //long 
        double dA01 = ps.triangle_area(iedge, kedge, edges);
        Uint medge = get_common_entry(node2edges[knode], node2edges[inode]);
        faces.push_back({{iedge, kedge, medge}, dA01});
        edge2faces[iedge].push_back(iface);
        edge2faces[kedge].push_back(iface);
        edge2faces[medge].push_back(iface);
        ++iface;
        //long 
        double dA02 = ps.triangle_area(iedge, ledge, edges);
        Uint nedge = get_common_entry(node2edges[lnode], node2edges[jnode]);
        faces.push_back({{iedge, ledge, nedge}, dA02});
        edge2faces[iedge].push_back(iface);
        edge2faces[ledge].push_back(iface);
        edge2faces[nedge].push_back(iface);
        ++iface;
        ++iedge;
      }
      assert(iedge == edges.size());
      if (verbose){
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
  return true;
}


#endif
