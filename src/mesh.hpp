#include "utils.hpp"

#ifndef __MESH_HPP
#define __MESH_HPP

using namespace std;

void print_mesh(const FacesType& faces, const EdgesType& edges,
                const Edge2FacesType& edge2faces,
                Vector3d* x_rw, const Uint Nrw){
  cout << "==================" << endl;
  cout << "nodes" << endl;
  for (Uint irw=0; irw < Nrw; ++irw){
    cout << irw << " "
         << x_rw[irw][0] << " "
         << x_rw[irw][1] << " "
         << x_rw[irw][2] << endl;
  }

  cout << "edges" << endl;
  for (size_t ie=0; ie < edges.size(); ++ie){
    cout << ie << ": " << edges[ie].first[0] << " " << edges[ie].first[1] << " " << edges[ie].second << endl;
  }

  cout << "faces" << endl;
  for (size_t ifa=0; ifa < faces.size(); ++ifa){
    cout << ifa << ": " << faces[ifa].first[0] << " " << faces[ifa].first[1] << " " << faces[ifa].first[2] << " " << faces[ifa].second << endl;
  }

  cout << "edge2faces" << endl;
  for (size_t ie2f = 0; ie2f < edge2faces.size(); ++ie2f){
    cout << ie2f << ": ";
    for (FacesListType::const_iterator vit=edge2faces[ie2f].begin();
         vit != edge2faces[ie2f].end(); ++vit){
      cout << *vit << " ";
    }
    cout << endl;
  }
  cout << "==================" << endl;
}

void dump_mesh(const FacesType& faces, const EdgesType& edges,
               const Edge2FacesType& edge2faces,
               Vector3d* x_rw, const Uint Nrw){
  ofstream nodesf("mesh.node");
  for (Uint irw=0; irw < Nrw; ++irw){
    nodesf << x_rw[irw][0] << " "
           << x_rw[irw][1] << " "
           << x_rw[irw][2] << endl;
  }
  nodesf.close();

  ofstream edgesf("mesh.edge");
  for (size_t ie=0; ie < edges.size(); ++ie){
    edgesf << edges[ie].first[0] << " " << edges[ie].first[1] << " " << edges[ie].second << endl;
  }
  nodesf.close();

  ofstream facesf("mesh.face");
  for (size_t ifa=0; ifa < faces.size(); ++ifa){
    facesf << faces[ifa].first[0] << " " << faces[ifa].first[1] << " "
           << faces[ifa].first[2] << " " << faces[ifa].second << endl;
  }
  facesf.close();

  ofstream e2ff("mesh.e2f");
  for (size_t ie2f = 0; ie2f < edge2faces.size(); ++ie2f){
    for (FacesListType::const_iterator vit=edge2faces[ie2f].begin();
         vit != edge2faces[ie2f].end(); ++vit){
      e2ff << *vit << " ";
    }
    e2ff << endl;
  }
  e2ff.close();
}

void posdata2txt(ofstream &pos_out,
                 Vector3d* x_rw,
                 Vector3d* u_rw,
                 const Uint Nrw){
  for (Uint irw=0; irw < Nrw; ++irw){
    pos_out << irw << " "
            << x_rw[irw][0] << " "
            << x_rw[irw][1] << " "
            << x_rw[irw][2] << " "
            << u_rw[irw][0] << " "
            << u_rw[irw][1] << " "
            << u_rw[irw][2] << std::endl;
  }
}

void tensor2hdf5(H5File* h5f, const string dsetname,
                 double* axx_rw, double* axy_rw, double* axz_rw,
                 double* ayx_rw, double* ayy_rw, double* ayz_rw,
                 double* azx_rw, double* azy_rw, double* azz_rw,
                 const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 3*3;
  DataSpace dspace(2, dims);
  double* data = new double[Nrw*3*3];
  for (Uint irw=0; irw < Nrw; ++irw){
    data[irw*3*3+0] = axx_rw[irw];
    data[irw*3*3+1] = axy_rw[irw];
    data[irw*3*3+2] = axz_rw[irw];
    data[irw*3*3+3] = ayx_rw[irw];
    data[irw*3*3+4] = ayy_rw[irw];
    data[irw*3*3+5] = ayz_rw[irw];
    data[irw*3*3+6] = azx_rw[irw];
    data[irw*3*3+7] = azy_rw[irw];
    data[irw*3*3+8] = azz_rw[irw];
  }
  DataSet dset = h5f->createDataSet(dsetname,
                                    PredType::NATIVE_DOUBLE,
                                    dspace);
  dset.write(data, PredType::NATIVE_DOUBLE);
}

void vector2hdf5(H5File* h5f, const string dsetname,
                 double* ax_rw, double* ay_rw, double* az_rw,
                 const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 3;
  DataSpace dspace(2, dims);
  double* data = new double[Nrw*3];
  for (Uint irw=0; irw < Nrw; ++irw){
    data[irw*3+0] = ax_rw[irw];
    data[irw*3+1] = ay_rw[irw];
    data[irw*3+2] = az_rw[irw];
  }
  DataSet dset = h5f->createDataSet(dsetname,
                                    PredType::NATIVE_DOUBLE,
                                    dspace);
  dset.write(data, PredType::NATIVE_DOUBLE);
}

void vector2hdf5(H5File* h5f, const string dsetname,
                 Vector3d* a_rw, const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 3;
  DataSpace dspace(2, dims);
  double* data = new double[Nrw*3];
  for (Uint irw=0; irw < Nrw; ++irw){
    for (Uint d=0; d<3; ++d){
      data[irw*3+d] = a_rw[irw][d];
    }
  }
  DataSet dset = h5f->createDataSet(dsetname,
                                    PredType::NATIVE_DOUBLE,
                                    dspace);
  dset.write(data, PredType::NATIVE_DOUBLE);
}


void scalar2hdf5(H5File* h5f, const string dsetname, double* c_rw,
                 const Uint Nrw){
  hsize_t dims[2];
  dims[0] = Nrw;
  dims[1] = 1;
  DataSpace dspace(2, dims);
  double* data = new double[Nrw];
  for (Uint irw=0; irw < Nrw; ++irw){
    data[irw] = c_rw[irw];
  }
  DataSet dset = h5f->createDataSet(dsetname,
                                    PredType::NATIVE_DOUBLE,
                                    dspace);
  dset.write(data, PredType::NATIVE_DOUBLE);
}

void mesh2hdf(H5File* h5f, const string groupname,
              Vector3d* x_rw,
              const FacesType& faces, const EdgesType& edges){
  // Faces
  if (faces.size() > 0){
    hsize_t faces_dims[2];
    faces_dims[0] = faces.size();
    faces_dims[1] = 3;
    DataSpace faces_dspace(2, faces_dims);
    Uint* faces_arr = new Uint[faces_dims[0]*faces_dims[1]];

    double* dA = new double[faces_dims[0]];
    double* dA0 = new double[faces_dims[0]];
    for (Uint iface=0; iface < faces_dims[0]; ++iface){
      set<Uint> unique_nodes;
      for (Uint i=0; i<3; ++i){
        Uint iedge = faces[iface].first[i];
        for (Uint j=0; j<2; ++j){
          Uint inode = edges[iedge].first[j];
          unique_nodes.insert(inode);
        }
      }
      Uint count = 0;
      for (set<Uint>::const_iterator sit=unique_nodes.begin();
           sit != unique_nodes.end(); ++sit){
        faces_arr[iface*faces_dims[1] + count] = *sit;
        ++count;
      }
      dA[iface] = area(iface, x_rw, faces, edges);
      dA0[iface] = faces[iface].second;
    }
    DataSet faces_dset = h5f->createDataSet(groupname + "/faces",
                                             PredType::NATIVE_ULONG,
                                             faces_dspace);
    faces_dset.write(faces_arr, PredType::NATIVE_ULONG);

    scalar2hdf5(h5f, groupname + "/dA", dA, faces_dims[0]);
    scalar2hdf5(h5f, groupname + "/dA0", dA0, faces_dims[0]);
  }
  // Edges
  else if (edges.size() > 0){
    hsize_t edges_dims[2];
    edges_dims[0] = edges.size();
    edges_dims[1] = 2;
    DataSpace edges_dspace(2, edges_dims);
    Uint* edges_arr = new Uint[edges_dims[0]*edges_dims[1]];

    double* dl = new double[edges_dims[0]];
    double* dl0 = new double[edges_dims[0]];
    for (Uint iedge=0; iedge < edges_dims[0]; ++iedge){
      for (Uint j=0; j<2; ++j){
        edges_arr[iedge*edges_dims[1] + j] = edges[iedge].first[j];
      }
      dl[iedge] = dist(edges[iedge].first[0], edges[iedge].first[1], x_rw);
      dl0[iedge] = edges[iedge].second;
    }

    DataSet edges_dset = h5f->createDataSet(groupname + "/edges",
                                             PredType::NATIVE_ULONG,
                                             edges_dspace);
    edges_dset.write(edges_arr, PredType::NATIVE_ULONG);

    scalar2hdf5(h5f, groupname + "/dl", dl, edges_dims[0]);
    scalar2hdf5(h5f, groupname + "/dl0", dl0, edges_dims[0]);
  }
}

array<Uint, 2> sort_edges(Uint inode, Uint kedge, Uint ledge,
                          EdgesType &edges){
  if (inode == edges[kedge].first[0] || inode == edges[kedge].first[1]){
    return {kedge, ledge};
  }
  return {ledge, kedge};
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
  cout << "Error: No common entry!" << endl;
  exit(0);
  return -1;
}

array<Uint, 3> get_close_entities(Uint iedge, Uint jedge, Uint kedge, Uint ledge,
                                 EdgesType &edges){
  Uint inode = edges[iedge].first[0];
  Uint knode;
  array<Uint, 2> mnedges;
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
    cout << "Error: Found no close entities." << endl;
    exit(0);
  }
  // cout << mnedges[0] << " " << mnedges[1] << endl;
  return {knode, mnedges[0], mnedges[1]};
}

void append_new_node(const Uint inode, const Uint jnode,
                     Vector3d* x_rw,
                     Vector3d* u_rw,
                     double* rho_rw, double* p_rw, double* c_rw,
                     double* H_rw, Vector3d* n_rw,
                     Vector3d* a_rw,
                     Uint& Nrw, const bool do_output_all,
                     Interpol *intp,
                     const double U0, const int int_order){
  x_rw[Nrw] = 0.5*(x_rw[inode]+x_rw[jnode]);

  c_rw[Nrw] = 0.5*(c_rw[inode]+c_rw[jnode]);
  H_rw[Nrw] = 0.5*(H_rw[inode]+H_rw[jnode]);

  n_rw[Nrw] = 0.5*(n_rw[inode]+n_rw[jnode]);
  n_rw[Nrw] /= n_rw[Nrw].norm();

  intp->probe(x_rw[Nrw]);

  u_rw[Nrw] = U0*intp->get_u();
  if (do_output_all){
    rho_rw[Nrw] = intp->get_rho();
    p_rw[Nrw] = intp->get_p();
  }
  // Second-order terms
  if (int_order >= 2){
    a_rw[Nrw] = U0*U0*intp->get_Ju() + U0*intp->get_a();
  }
  ++Nrw;
}

Uint sheet_refinement(FacesType &faces,
                      EdgesType &edges,
                      Edge2FacesType &edge2faces,
                      Node2EdgesType &node2edges,
                      Vector3d* x_rw,
                      Vector3d* u_rw,
                      double* rho_rw, double* p_rw, double* c_rw,
                      double* H_rw, Vector3d* n_rw,
                      Vector3d* a_rw,
                      Uint &Nrw, const Uint Nrw_max, const double ds_max,
                      const bool do_output_all,
                      Interpol *intp,
                      const double U0,
                      const int int_order,
                      const double curv_refine_factor){
  bool changed;
  Uint n_add = 0;
  do {
    changed = false;

    vector<double> ds_ratio_;
    for (EdgesType::const_iterator edgeit = edges.begin();
         edgeit != edges.end(); ++edgeit){
      Uint inode = edgeit->first[0];
      Uint jnode = edgeit->first[1];
      double ds = dist(inode, jnode, x_rw);
      // cout << inode << " " << jnode << " " << ds << endl;
      double kappa = 0.5*(abs(H_rw[inode]) + abs(H_rw[jnode]));
      double ds_max_loc = ds_max/(1.0 + curv_refine_factor*kappa);

      ds_ratio_.push_back(ds/ds_max_loc);
    }
    vector<size_t> ids_ = argsort_descending(ds_ratio_);

    for (vector<size_t>::iterator itedge = ids_.begin();
         itedge != ids_.end(); ++itedge){

      Uint iedge = *itedge;
      if (Nrw >= Nrw_max){
        // No more points can fit
        break;
      }
      double ds_ratio = ds_ratio_[iedge];
      if (ds_ratio < 1.0)
        break;
      changed = true;

      Uint inode = edges[iedge].first[0];
      Uint jnode = edges[iedge].first[1];
      double ds0 = edges[iedge].second;

      Uint new_inode = Nrw;
      // Add point
      append_new_node(inode, jnode, x_rw, u_rw,
                      rho_rw, p_rw, c_rw, H_rw, n_rw,
                      a_rw, Nrw, do_output_all, intp, U0, int_order);
      ++n_add;
      //
      edges[iedge].first[1] = new_inode;
      Uint new_iedge = edges.size();
      edges.push_back({{new_inode, jnode}, ds0/2});
      edge2faces.push_back({});
      // Append new node to node2edges list
      node2edges.push_back({iedge, new_iedge});
      // Modify existing entry
      replace(node2edges[jnode].begin(), node2edges[jnode].end(), iedge, new_iedge);

      for (FacesListType::iterator itface = edge2faces[iedge].begin();
           itface != edge2faces[iedge].end(); ++itface){
        Uint jedge = faces[*itface].first[0];
        Uint kedge = faces[*itface].first[1];
        Uint ledge = faces[*itface].first[2];
        double dA0 = faces[*itface].second;
        Uint new_jedge = edges.size();
        array<Uint, 3> close_entities = get_close_entities(iedge, jedge, kedge, ledge, edges);
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
  } while (changed);
  return n_add;
}

Uint strip_refinement(EdgesType &edges,
                      Node2EdgesType &node2edges,
                      Vector3d* x_rw,
                      Vector3d* u_rw,
                      double* rho_rw, double* p_rw, double* c_rw,
                      double* H_rw, Vector3d* n_rw,
                      Vector3d* a_rw,
                      Uint &Nrw, const Uint Nrw_max, const double ds_max,
                      const bool do_output_all,
                      Interpol *intp,
                      const double U0, const int int_order,
                      const double curv_refine_factor){
  Uint n_add = 0;
  for (Uint iedge=0; iedge < edges.size(); ++iedge){
    Uint inode = edges[iedge].first[0];
    Uint jnode = edges[iedge].first[1];
    double ds0 = edges[iedge].second;
    double ds = dist(inode, jnode, x_rw);
    //double kappa = sqrt(abs(H_rw[inode]*H_rw[jnode]));
    double kappa = 0.5*(abs(H_rw[inode]) + abs(H_rw[jnode]));
    double ds_max_loc = ds_max/(1.0 + curv_refine_factor*kappa);
    if (ds > ds_max_loc && Nrw < Nrw_max){
      edges[iedge].first[1] = Nrw;
      edges[iedge].second = ds0/2;
      Uint new_iedge = edges.size();
      edges.push_back({{Nrw, jnode}, ds0/2});

      node2edges.push_back({iedge, new_iedge});
      replace(node2edges[jnode].begin(), node2edges[jnode].end(), iedge, new_iedge);

      append_new_node(inode, jnode, x_rw, u_rw,
                      rho_rw, p_rw, c_rw, H_rw, n_rw,
                      a_rw, Nrw, do_output_all, intp, U0, int_order);
      ++n_add;
    }
  }
  return n_add;
}

Uint refinement(FacesType &faces,
                EdgesType &edges,
                Edge2FacesType &edge2faces,
                Node2EdgesType &node2edges,
                Vector3d* x_rw,
                Vector3d* u_rw,
                double* rho_rw, double* p_rw, double* c_rw,
                double* H_rw, Vector3d* n_rw,
                Vector3d* a_rw,
                Uint &Nrw, const Uint Nrw_max, const double ds_max,
                const bool do_output_all,
                Interpol *intp,
                const double U0, const int int_order,
                const double curv_refine_factor){
  Uint n_add = 0;
  if (faces.size() > 0){
    n_add = sheet_refinement(faces, edges, edge2faces, node2edges,
                             x_rw,
                             u_rw,
                             rho_rw, p_rw, c_rw, H_rw, n_rw,
                             a_rw,
                             Nrw, Nrw_max, ds_max,
                             do_output_all, intp,
                             U0, int_order, curv_refine_factor);
  }
  else {
    n_add = strip_refinement(edges,
                             node2edges,
                             x_rw,
                             u_rw,
                             rho_rw, p_rw, c_rw, H_rw, n_rw,
                             a_rw,
                             Nrw, Nrw_max, ds_max,
                             do_output_all, intp,
                             U0, int_order,
                             curv_refine_factor);
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

void print(const vector<Uint> vec){
  for (vector<Uint>::const_iterator vit = vec.begin();
       vit != vec.end(); ++vit){
    cout << *vit << " ";
  }
  cout << endl;
}

void print(const set<Uint> vec){
  for (set<Uint>::const_iterator vit = vec.begin();
       vit != vec.end(); ++vit){
    cout << *vit << " ";
  }
  cout << endl;
}

void print(const list<Uint> vec){
  for (list<Uint>::const_iterator vit = vec.begin();
       vit != vec.end(); ++vit){
    cout << *vit << " ";
  }
  cout << endl;
}

void get_conodes(set<Uint> &inodes,
                 map<Uint, Uint> &icoedges,
                 const Uint inode,
                 const Node2EdgesType &node2edges,
                 const EdgesType &edges){
  //cout << "Getting conodes of " << inode << endl;
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

set<Uint> get_cofaces(const Uint iedge,
                      const Edge2FacesType &edge2faces,
                      const FacesType &faces){
  FacesListType ifaces = edge2faces[iedge];
  set <Uint> jfaces;
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

bool get_new_pos(Vector3d &xi,
                 Vector3d* x_rw,
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
    xi = 0.5*(x_rw[inode] + x_rw[jnode]);
  }
  else if (inode_is_border) {
    xi = x_rw[inode];
  }
  else {
    assert(jnode_is_border);
    xi = x_rw[jnode];
  }
  return true;
}

bool normals_are_ok(const Uint iedge,
                    const Vector3d &x,
                    Vector3d* x_rw,
                    const set<Uint> &jfaces,
                    const FacesType &faces,
                    const EdgesType &edges){
  Uint inode = edges[iedge].first[0];
  Uint jnode = edges[iedge].first[1];

  for (set<Uint>::const_iterator jfaceit=jfaces.begin();
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

vector<Uint> get_incident_faces(const Uint iedge,
                                const EdgesType &edges,
                                const Edge2FacesType &edge2faces,
                                const Node2EdgesType &node2edges){
  Uint inode = edges[iedge].first[0];
  Uint jnode = edges[iedge].first[1];
  set<Uint> kfaces_set;
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
  vector<Uint> kfaces(kfaces_set.begin(), kfaces_set.end());
  return kfaces;
}

vector<double> areas(const vector<Uint> &kfaces,
                          Vector3d* x_rw,
                          const FacesType &faces, const EdgesType &edges){
  vector<double> a;
  for (vector<Uint>::const_iterator faceit=kfaces.begin();
       faceit != kfaces.end(); ++faceit){
    a.push_back(area(*faceit, x_rw, faces, edges));
  }
  return a;
}

bool collapse_edge(const Uint iedge,
                   FacesType &faces,
                   EdgesType &edges,
                   Edge2FacesType &edge2faces,
                   Node2EdgesType &node2edges,
                   vector<bool> &face_isactive,
                   vector<bool> &edge_isactive,
                   vector<bool> &node_isactive,
                   Vector3d* x_rw,
                   Vector3d* u_rw,
                   double* rho_rw, double* p_rw, double* c_rw,
                   double* H_rw, Vector3d* n_rw,
                   Vector3d* a_rw,
                   const bool do_output_all,
                   Interpol *intp,
                   const double U0, const int int_order){

  assert(edge_isactive[iedge]);

  Uint inode = edges[iedge].first[0];
  Uint jnode = edges[iedge].first[1];
  // double ds0 = edges[iedge].second;
  set<Uint> inodes;
  set<Uint> jnodes;
  map<Uint, Uint> icoedges;
  map<Uint, Uint> jcoedges;
  get_conodes(inodes, icoedges, inode, node2edges, edges);
  get_conodes(jnodes, jcoedges, jnode, node2edges, edges);
  vector<Uint> joint_nodes;
  set_intersection(inodes.begin(), inodes.end(),
                   jnodes.begin(), jnodes.end(),
                   back_inserter(joint_nodes));

  //cout << "iedge=" << iedge << endl;
  //print(inodes);
  //print(jnodes);
  //print(joint_nodes);

  // Check if it has too many/too few common joint nodes.
  if (joint_nodes.size() > 2){
    // Wrong number of joint nodes
    // cout << "Wrong number of joint nodes." << endl;
    return false;
  }
  Vector3d x;
  bool new_node_isgood = get_new_pos(x, x_rw, iedge, edges, edge2faces, node2edges);
  if (!new_node_isgood){
    // Corner edge
    // cout << "Corner node!" << endl;
    return false;
  }
  //cout << "New pos: " << x << " " << y << " " << z << endl;

  set<Uint> jfaces = get_cofaces(iedge, edge2faces, faces);

  //print(jfaces);

  double dA_res = 0.;
  double dA0_res = 0.;
  for (FacesListType::const_iterator faceit=edge2faces[iedge].begin();
       faceit != edge2faces[iedge].end(); ++faceit){
    dA_res += area(*faceit, x_rw, faces, edges);
    dA0_res += faces[*faceit].second;
    //faces[*faceit].second = 0.;
  }

  vector<Uint> kfaces = get_incident_faces(iedge, edges, edge2faces, node2edges);
  vector<double> dAs_old = areas(kfaces, x_rw, faces, edges);

  // assert(jfaces.size()==4 || jfaces.size()==2);

  if (!normals_are_ok(iedge, x, x_rw,
                      jfaces, faces, edges)){
    // Normals are not ok
    // cout << "Normals are not OK!" << endl;
    return false;
  }

  intp->probe(x);  //

  Uint irws[2] = {inode, jnode};
  for (Uint i=0; i<2; ++i){
    Uint irw = irws[i];
    x_rw[irw] = x;
    u_rw[irw] = U0*intp->get_u();

    if (do_output_all){
      rho_rw[irw] = intp->get_rho();
      p_rw[irw] = intp->get_p();
    }
    // Second-order terms
    if (int_order >= 2){
      double U02 = U0*U0;
      a_rw[irw] = U02*intp->get_Ju() + U0*intp->get_a();
    }
    c_rw[irw] = 0.5*(c_rw[inode]+c_rw[jnode]);
    H_rw[irw] = 0.5*(H_rw[inode]+H_rw[jnode]);
    n_rw[irw] = 0.5*(n_rw[inode]+n_rw[jnode]);
  }

  for (EdgesListType::const_iterator edgeit=node2edges[max(inode, jnode)].begin();
       edgeit != node2edges[max(inode, jnode)].end(); ++edgeit){
    for (Uint j=0; j<2; ++j){
      if (edges[*edgeit].first[j] == max(inode, jnode)){
        edges[*edgeit].first[j] = min(inode, jnode);
      }
    }
  }

  vector<Uint> iedges_vec;
  set_symmetric_difference(node2edges[inode].begin(),
                           node2edges[inode].end(),
                           node2edges[jnode].begin(),
                           node2edges[jnode].end(),
                           back_inserter(iedges_vec));
  set<Uint> iedges_set(iedges_vec.begin(), iedges_vec.end());
  //print(node2edges[inode]);
  //print(node2edges[jnode]);
  //print(iedges_set);

  map<Uint, Uint> replace_edges;
  for (vector<Uint>::const_iterator joint_node_it=joint_nodes.begin();
       joint_node_it != joint_nodes.end(); ++joint_node_it){
    Uint kedge_min = min(icoedges[*joint_node_it], jcoedges[*joint_node_it]);
    Uint kedge_max = max(icoedges[*joint_node_it], jcoedges[*joint_node_it]);
    replace_edges[kedge_max] = kedge_min;
    assert(edge_isactive[kedge_max]);
    edge_isactive[kedge_max] = false;
    if (contains(iedges_set, kedge_max))
      iedges_set.erase(kedge_max);
    //cout << "kedge_max=" << kedge_max << ", kedge_min=" << kedge_min << endl;
  }
  //print(iedges_set);

  assert(edge_isactive[iedge]);
  edge_isactive[iedge] = false;

  for (set<Uint>::const_iterator jfaceit=jfaces.begin();
       jfaceit != jfaces.end(); ++jfaceit){
    //cout << "faces[" << *jfaceit << "].first="
    //     << faces[*jfaceit].first[0] << " "
    //     << faces[*jfaceit].first[1] << " "
    //     << faces[*jfaceit].first[2] << endl;
    for (Uint j=0; j<3; ++j){
      if (contains(replace_edges, faces[*jfaceit].first[j])){
        faces[*jfaceit].first[j] = replace_edges[faces[*jfaceit].first[j]];
      }
    }
    //cout << "faces[" << *jfaceit << "].first="
    //     << faces[*jfaceit].first[0] << " "
    //     << faces[*jfaceit].first[1] << " "
    //     << faces[*jfaceit].first[2] << endl;
  }
  vector<double> dAs_new = areas(kfaces, x_rw, faces, edges);
  double wsum = 0.;
  vector<double> w_;
  for (Uint k=0; k < kfaces.size(); ++k){
    //Uint kface = kfaces[k];
    // cout << faces[kface].second << "+=" << dA0_res << "/" << dA_res << "*(" << dAs_new[k] << "-" << dAs_old[k] << ")" << endl;
    // cout << dA0_res << " " << dA_res << endl;
    // if (dA_res > 0.0){
    //     faces[kface].second += dA0_res/dA_res*(dAs_new[k]-dAs_old[k]);
    // faces[kface].second += dA0_res/kfaces.size();
    double w = dAs_new[k]-dAs_old[k]; // May give negative mass
    w = max(w, 0.0);
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

  //print(node2edges[min(inode, jnode)]);
  node2edges[min(inode, jnode)].assign(iedges_set.begin(), iedges_set.end());
  //print(node2edges[min(inode, jnode)]);
  node2edges[max(inode, jnode)].clear();
  assert(node_isactive[max(inode, jnode)]);
  node_isactive[max(inode, jnode)] = false;
  for (vector<Uint>::const_iterator joint_node_it=joint_nodes.begin();
       joint_node_it != joint_nodes.end(); ++joint_node_it){
    set<Uint> tmp_set(node2edges[*joint_node_it].begin(),
                      node2edges[*joint_node_it].end());
    //print(node2edges[*joint_node_it]);
    for (map<Uint, Uint>::const_iterator mit=replace_edges.begin();
         mit != replace_edges.end(); ++mit){
      tmp_set.erase(mit->first);
    }
    node2edges[*joint_node_it].assign(tmp_set.begin(), tmp_set.end());
    //print(node2edges[*joint_node_it]);
  }

  for (map<Uint, Uint>::const_iterator mit=replace_edges.begin();
       mit != replace_edges.end(); ++mit){
    Uint kedge_max = mit->first;
    Uint kedge_min = mit->second;
    vector<Uint> ifaces;
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

void remove_faces(FacesType &faces, const vector<bool> face_isactive){
  assert(faces.size() == face_isactive.size());
  for (int i=int(face_isactive.size())-1; i >= 0; --i){
    if (!face_isactive[i])
      faces.erase(faces.begin()+i);
  }
}

void remove_edges(FacesType &faces, EdgesType &edges,
                  const vector<bool> &edge_isactive){
  assert(edges.size() == edge_isactive.size());
  vector<Uint> used_edges;
  for (Uint i=0; i<edges.size(); ++i){
    used_edges.push_back(i);
  }
  for (int i=edge_isactive.size()-1; i >= 0; --i){
    if (!edge_isactive[i]){
      edges.erase(edges.begin()+i);
      used_edges.erase(used_edges.begin()+i);
    }
  }
  map<Uint, Uint> replace_edges;
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
}

void remove_nodes(EdgesType& edges,
                  Vector3d* x_rw,
                  Vector3d* u_rw,
                  double* rho_rw, double* p_rw, double* c_rw,
                  double* H_rw, Vector3d* n_rw,
                  Vector3d* a_rw,
                  Uint& Nrw,
                  const vector<bool> &node_isactive,
                  const bool do_output_all,
                  const int int_order
                  ){
  assert(Nrw == node_isactive.size());
  vector<Uint> used_nodes;
  for (Uint i=0; i<Nrw; ++i){
    if (node_isactive[i])
      used_nodes.push_back(i);
  }
  map<Uint, Uint> replace_nodes;
  for (Uint i=0; i<used_nodes.size(); ++i){
    replace_nodes[used_nodes[i]] = i;
    Uint j = used_nodes[i];
    x_rw[i] = x_rw[j];
    u_rw[i] = u_rw[j];

    if (do_output_all){
      rho_rw[i] = rho_rw[j];
      p_rw[i] = p_rw[j];
    }
    c_rw[i] = c_rw[j];
    H_rw[i] = H_rw[j];
    n_rw[i] = n_rw[j];

    if (int_order >= 2){
      a_rw[i] = a_rw[j];
    }
  }

  Nrw = used_nodes.size();

  for (Uint iedge=0; iedge < edges.size(); ++iedge){
    for (Uint j=0; j<2; ++j){
      Uint inode = edges[iedge].first[j];
      if (contains(replace_nodes, inode)){
        edges[iedge].first[j] = replace_nodes[inode];
      }
    }
  }
}

void remove_unused_edges(FacesType &faces, EdgesType &edges){
  vector<bool> edge_isactive(edges.size(), false);
  for (Uint i=0; i < faces.size(); ++i){
    for (Uint k=0; k < 3; ++k){
      edge_isactive[faces[i].first[k]] = true;
    }
  }
  remove_edges(faces, edges, edge_isactive);
}

void remove_unused_nodes(EdgesType &edges,
                         Vector3d* x_rw,
                         Vector3d* u_rw,
                         double* rho_rw, double* p_rw, double* c_rw,
                         double* H_rw, Vector3d* n_rw,
                         Vector3d* a_rw,
                         Uint &Nrw,
                         const bool do_output_all,
                         const int int_order){
  vector<bool> node_isactive(Nrw, false);
  for (Uint i=0; i < edges.size(); ++i){
    for (Uint k=0; k < 2; ++k){
      node_isactive[edges[i].first[k]] = true;
    }
  }
  remove_nodes(edges, x_rw, u_rw,
               rho_rw, p_rw, c_rw, H_rw, n_rw, a_rw, Nrw,
               node_isactive, do_output_all, int_order);
}

void assert_equal(const Edge2FacesType &a,
                  const Edge2FacesType &b){
  assert(a.size() == b.size());
  for (Uint i=0; i < a.size(); ++i){
    assert(a[i].size() == b[i].size());
  }
  exit(0);
}

Uint sheet_coarsening(FacesType &faces,
                      EdgesType &edges,
                      Edge2FacesType &edge2faces,
                      Node2EdgesType &node2edges,
                      Vector3d* x_rw,
                      Vector3d* u_rw,
                      double* rho_rw, double* p_rw, double* c_rw,
                      double* H_rw, Vector3d* n_rw,
                      Vector3d* a_rw,
                      Uint &Nrw, const double ds_min,
                      const bool do_output_all, Interpol *intp,
                      const double U0, const int int_order,
                      const double curv_refine_factor){
  bool changed;
  Uint n_coll = 0;
  Uint iedge;

  vector<bool> face_isactive(faces.size(), true);
  vector<bool> edge_isactive(edges.size(), true);
  vector<bool> node_isactive(Nrw, true);

  do {
    changed = false;
    iedge = 0;
    while (iedge < edges.size()){
      Uint inode = edges[iedge].first[0];
      Uint jnode = edges[iedge].first[1];
      // double ds0 = edges[iedge].second;
      double ds = dist(inode, jnode, x_rw);
      double kappa = sqrt(abs(H_rw[inode]*H_rw[jnode]));
      double ds_min_loc = ds_min/(1.0 + curv_refine_factor*kappa);
      if (ds < ds_min_loc && edge_isactive[iedge]){
        // cout << "TRYING to collapse edge " << iedge
        //      << ". ds = " << ds << " and ds_min = " << ds_min << endl;
        bool coll = collapse_edge(iedge, faces, edges,
                                  edge2faces,
                                  node2edges,
                                  face_isactive,
                                  edge_isactive,
                                  node_isactive,
                                  x_rw,
                                  u_rw,
                                  rho_rw, p_rw, c_rw, H_rw, n_rw,
                                  a_rw,
                                  do_output_all, intp,
                                  U0, int_order);
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
  remove_unused_edges(faces, edges);
  remove_unused_nodes(edges, x_rw,
                      u_rw,
                      rho_rw, p_rw, c_rw, H_rw, n_rw,
                      a_rw,
                      Nrw, do_output_all, int_order);

  compute_edge2faces(edge2faces, faces, edges);
  compute_node2edges(node2edges, edges, Nrw);
  return n_coll;
}

Uint strip_coarsening(EdgesType &edges,
                      Node2EdgesType &node2edges,
                      Vector3d* x_rw,
                      Vector3d* u_rw,
                      double* rho_rw, double* p_rw, double* c_rw,
                      double* H_rw, Vector3d* n_rw,
                      Vector3d* a_rw,
                      Uint &Nrw, const double ds_min,
                      const bool do_output_all, Interpol *intp,
                      const double U0, const int int_order,
                      const double curv_refine_factor){
  bool changed;
  Uint n_coll = 0;
  Uint iedge;

  vector<bool> edge_isactive(edges.size(), true);
  vector<bool> node_isactive(Nrw, true);

  do {
    changed = false;
    iedge = 0;
    while (iedge < edges.size()){
      if (edge_isactive[iedge]){
        Uint inode = edges[iedge].first[0];
        Uint jnode = edges[iedge].first[1];
        double ds0 = edges[iedge].second;
        double ds = dist(inode, jnode, x_rw);
        double kappa = sqrt(abs(H_rw[inode]*H_rw[jnode]));
        double ds_min_loc = ds_min/(1.0 + curv_refine_factor*kappa);
        if (ds < ds_min_loc && edge_isactive[iedge]){
          edge_isactive[iedge] = false;

          Uint new_inode = min(inode, jnode);
          Uint old_inode = max(inode, jnode);

          node_isactive[old_inode] = false;

          Vector3d pos_new;
          bool inode_is_border = node2edges[inode].size() > 1;
          bool jnode_is_border = node2edges[jnode].size() > 1;
          if ((inode_is_border && jnode_is_border) || (!inode_is_border && !jnode_is_border)){
            pos_new = 0.5*(x_rw[inode] + x_rw[jnode]);
          }
          else if (inode_is_border){
            pos_new = x_rw[inode];
          }
          else {
            pos_new = x_rw[jnode];
          }
          x_rw[inode] = pos_new;
          x_rw[jnode] = pos_new;

          vector<Uint> jedges(node2edges[old_inode].begin(),
                              node2edges[old_inode].end());
          Uint jedge = get_other(jedges[0], jedges[1], iedge);
          replace(edges[jedge].first.begin(), edges[jedge].first.end(),
                  old_inode, new_inode);
          vector<Uint> kedges(node2edges[new_inode].begin(),
                              node2edges[new_inode].end());
          Uint kedge = get_other(kedges[0], kedges[1], iedge);

          edges[jedge].second += ds0/2;
          edges[kedge].second += ds0/2;

          node2edges[old_inode].clear();
          replace(node2edges[new_inode].begin(),
                  node2edges[new_inode].end(), iedge, jedge);

          intp->probe(x_rw[new_inode]);
          u_rw[new_inode] = U0*intp->get_u();
          if (do_output_all){
            rho_rw[new_inode] = intp->get_rho();
            p_rw[new_inode] = intp->get_p();
          }
          c_rw[new_inode] = 0.5*(c_rw[inode]+c_rw[jnode]);
          H_rw[new_inode] = 0.5*(H_rw[inode]+H_rw[jnode]);
          n_rw[new_inode] = 0.5*(n_rw[inode]+n_rw[jnode]);
          n_rw[new_inode] /= n_rw[new_inode].norm();
          if (int_order >= 2){
            a_rw[new_inode] = U0*U0*intp->get_Ju() + U0*intp->get_a();
          }

          changed = true;
          ++n_coll;
          //--iedge;
        }
      }
      ++iedge;
    }
  } while(changed);

  FacesType faces;
  remove_edges(faces, edges, edge_isactive);
  remove_unused_nodes(edges, x_rw,
                      u_rw,
                      rho_rw, p_rw, c_rw, H_rw, n_rw,
                      a_rw,
                      Nrw, do_output_all, int_order);
  compute_node2edges(node2edges, edges, Nrw);
  return n_coll;
}

Uint coarsening(FacesType &faces,
                EdgesType &edges,
                Edge2FacesType &edge2faces,
                Node2EdgesType &node2edges,
                Vector3d* x_rw,
                Vector3d* u_rw,
                double* rho_rw, double* p_rw, double* c_rw,
                double* H_rw, Vector3d* n_rw,
                Vector3d* a_rw,
                Uint &Nrw, const double ds_min,
                const bool do_output_all,
                Interpol *intp,
                const double U0, const int int_order,
                const double curv_refine_factor){
  if (faces.size() > 0){
    return sheet_coarsening(faces, edges,
                            edge2faces, node2edges,
                            x_rw, u_rw,
                            rho_rw, p_rw, c_rw, H_rw, n_rw,
                            a_rw, Nrw, ds_min,
                            do_output_all, intp,
                            U0, int_order, curv_refine_factor);
  }
  else{
    return strip_coarsening(edges, node2edges,
                            x_rw, u_rw, rho_rw, p_rw, c_rw, H_rw, n_rw,
                            a_rw, Nrw, ds_min,
                            do_output_all, intp,
                            U0, int_order, curv_refine_factor);
  }
  return 0;
}

array<double, 3> mixed_area_contrib(const double ang0,
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
                           vector<double> &mixed_areas,
                           vector<Vector3d> &face_normals,
                           const FacesType &faces,
                           const EdgesType &edges,
                           const Edge2FacesType &edge2faces,
                           Vector3d* x_rw,
                           const Uint Nrw){
  interior_ang.clear();
  mixed_areas.clear();
  face_normals.clear();
  mixed_areas.assign(Nrw, 0.);
  // face2nodes.clear();
  set<Uint> not_visited;
  for (Uint iface=0; iface < faces.size(); ++iface){
    not_visited.insert(iface);
    set<Uint> inodes_set;
    for (Uint i=0; i < 3; ++i){
      for (Uint j=0; j < 2; ++j){
        inodes_set.insert(edges[faces[iface].first[i]].first[j]);
      }
    }
    assert(inodes_set.size()==3);
    vector<Uint> inodes(inodes_set.begin(), inodes_set.end());
    array<Vector3d, 3> v;
    for (Uint i=0; i<3; ++i){
      Uint inode = inodes[i];
      v[i] = x_rw[inode];
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

    array<double, 3> a_loc = mixed_area_contrib(ang0, ang1, ang2, s01, s02, s12);
    for (Uint i=0; i<3; ++i){
      mixed_areas[inodes[i]] += a_loc[i];
    }
    Vector3d n_loc = get_normal(iface, faces, edges, x_rw);
    face_normals.push_back(n_loc);
  }
  set<Uint> to_visit;
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

void compute_sheet_curv(double* H_rw,
                        Vector3d* n_rw,
                        const FacesType &faces,
                        const EdgesType &edges,
                        const Edge2FacesType &edge2faces,
                        const Node2EdgesType &node2edges,
                        Vector3d* x_rw,
                        const Uint Nrw,
                        const InteriorAnglesType &interior_ang,
                        const vector<double> &mixed_areas,
                        const vector<Vector3d> &face_normals
                        ){
  vector<double> edge_w(edges.size(), 0.);
  for (Uint iedge=0; iedge < edges.size(); ++iedge){
    for (FacesListType::const_iterator faceit = edge2faces[iedge].begin();
         faceit != edge2faces[iedge].end(); ++faceit){
      vector<Uint> other_edges;
      for (Uint i=0; i<3; ++i){
        Uint jedge = faces[*faceit].first[i];
        if (iedge != jedge)
          other_edges.push_back(jedge);
      }
      assert(other_edges.size()==2);
      Uint inode = get_intersection(edges[other_edges[0]].first,
                                    edges[other_edges[1]].first);
      assert(contains(interior_ang[*faceit], inode));
      map<Uint, double> angles = interior_ang[*faceit];
      edge_w[iedge] += 1.0/tan(angles[inode]);
    }
  }
  for (Uint irw=0; irw<Nrw; ++irw){
    n_rw[irw] = {0., 0., 0.};
  }
  for (Uint iface=0; iface < interior_ang.size(); ++iface){
    map<Uint, double> angles = interior_ang[iface];
    for (map<Uint, double>::const_iterator angit=angles.begin();
         angit != angles.end(); ++angit){
      n_rw[angit->first] += angit->second*face_normals[iface];
    }
  }
  for (Uint irw=0; irw<Nrw; ++irw){
    n_rw[irw] /= n_rw[irw].norm();
  }

  for (Uint inode=0; inode<Nrw; ++inode){
    Vector3d lapl_v(0., 0., 0.);
    for (EdgesListType::const_iterator edgeit=node2edges[inode].begin();
         edgeit != node2edges[inode].end(); ++edgeit){
      Uint jnode = get_other(edges[*edgeit].first[0],
                             edges[*edgeit].first[1], inode);
      Vector3d dv = x_rw[jnode]-x_rw[inode];
      lapl_v += edge_w[*edgeit]*dv;
    }
    lapl_v /= 2*mixed_areas[inode];
    H_rw[inode] = 0.5*lapl_v.dot(n_rw[inode]);
  }
}

void compute_strip_curv(double* H_rw,
                        Vector3d* n_rw,
                        const EdgesType &edges,
                        const Node2EdgesType &node2edges,
                        Vector3d* x_rw, const Uint Nrw){
  for (Uint inode=0; inode<Nrw; ++inode){
    H_rw[inode] = 0.;
    n_rw[inode] = {1., 0., 0.};
    if (node2edges[inode].size() == 2){
      EdgesListType::const_iterator edgeit = node2edges[inode].begin();
      Uint jnode = get_other(edges[*edgeit].first[0],
                             edges[*edgeit].first[1], inode);
      ++edgeit;
      Uint knode = get_other(edges[*edgeit].first[0],
                             edges[*edgeit].first[1], inode);

      double R = circumcenter(x_rw[inode], x_rw[jnode], x_rw[knode]);

      H_rw[inode] = 1./R;

      //Vector3d n_loc = dij/dij.norm() + dik/dik.norm();
      //if (n_loc.norm() > 0)
      //  n_rw[inode] = n_loc/n_loc.norm();
    }
  }
}

void compute_mean_curv(double* H_rw,
                       Vector3d* n_rw,
                       const FacesType &faces,
                       const EdgesType &edges,
                       const Edge2FacesType &edge2faces,
                       const Node2EdgesType &node2edges,
                       Vector3d* x_rw,
                       const Uint Nrw,
                       const InteriorAnglesType &interior_ang,
                       const vector<double> &mixed_areas,
                       const vector<Vector3d> &face_normals
                       ){
  if (faces.size() > 1){
    compute_sheet_curv(H_rw, n_rw,
                       faces, edges,
                       edge2faces, node2edges,
                       x_rw, Nrw, interior_ang,
                       mixed_areas, face_normals);
  }
  else {
    compute_strip_curv(H_rw, n_rw, edges, node2edges, x_rw, Nrw);
  }
}

#endif
