#include <iostream>
#include "io.hpp"
#include "utils.hpp"
#include "Parameters.hpp"
#include "Interpol.hpp"
#include "distribute.hpp"
#include <vector>
#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <sstream>
#include <random>
#include <math.h>
#include <set>
#include "H5Cpp.h"

using namespace H5;


void write_stats(ofstream &statfile,
                 const double t,
                 double* x_rw, double* y_rw, double* z_rw,
                 double* ux_rw, double* uy_rw, double* uz_rw,
                 const Uint Nrw,
                 FacesType &faces,
                 EdgesType &edges,
                 const double ds_max,
                 const bool do_dump_hist,
                 const string histfolder,
                 const unsigned long int n_accepted,
                 const unsigned long int n_declined
                 ){
  double x_mean = 0.;
  double y_mean = 0.;
  double z_mean = 0.;
  double dx2_mean = 0.;
  double dy2_mean = 0.;
  double dz2_mean = 0.;
  double ux_mean = 0.;
  double uy_mean = 0.;
  double uz_mean = 0.;
  for (Uint irw=0; irw < Nrw; ++irw){
    // Sample mean
    x_mean += x_rw[irw]/Nrw;
    y_mean += y_rw[irw]/Nrw;
    z_mean += z_rw[irw]/Nrw;
    ux_mean += ux_rw[irw]/Nrw;
    uy_mean += uy_rw[irw]/Nrw;
    uz_mean += uz_rw[irw]/Nrw;
  }
  for (Uint irw=0; irw < Nrw; ++irw){
    // Sample variance
    dx2_mean += pow(x_rw[irw]-x_mean, 2)/(Nrw-1);
    dy2_mean += pow(y_rw[irw]-y_mean, 2)/(Nrw-1);
    dz2_mean += pow(z_rw[irw]-z_mean, 2)/(Nrw-1);
  }
  statfile << t << "\t"                   //  1
           << x_mean << "\t"              //  2
           << dx2_mean << "\t"            //  3
           << y_mean << "\t"              //  4
           << dy2_mean << "\t"            //  5
           << z_mean << "\t"              //  6
           << dz2_mean << "\t"            //  7
           << ux_mean << "\t"             //  8
           << uy_mean << "\t"             //  9
           << uz_mean << "\t"             // 10
           << Nrw << "\t"                 // 11
           << n_accepted << "\t"          // 12
           << n_declined << "\t";          // 13

  if (edges.size() > 0){
    // Strip method
    double s = 0.;
    Uint n_too_long = 0;
    double logelong_wmean = 0.;
    double logelong_w0mean = 0.;
    double s0 = 0.;
    vector<array<double, 3>> logelong_vec;
    for (EdgesType::const_iterator edgeit = edges.begin();
         edgeit != edges.end(); ++edgeit){
      int inode = edgeit->first[0];
      int jnode = edgeit->first[1];
      double ds0 = edgeit->second;
      double ds = dist(inode, jnode, x_rw, y_rw, z_rw);
      if (ds > ds_max){
        ++n_too_long;
      }
      double logelong = log(ds/ds0);
      logelong_wmean += logelong*ds;
      logelong_w0mean += logelong*ds0;
      logelong_vec.push_back({logelong, ds, ds0});
      s += ds;
      s0 += ds0;
    }
    logelong_wmean /= s;
    logelong_w0mean /= s0;
    statfile << n_too_long << "\t";          // 14

    if (faces.size() == 0){
      double logelong_wvar = 0.;
      double logelong_w0var = 0.;
      for (vector<array<double, 3>>::const_iterator lit = logelong_vec.begin();
           lit != logelong_vec.end(); ++lit){
        logelong_wvar += pow((*lit)[0]-logelong_wmean, 2)*(*lit)[1];
        logelong_w0var += pow((*lit)[0]-logelong_w0mean, 2)*(*lit)[2];
      }
      logelong_wvar /= s;
      logelong_w0var /= s0;
      double s_ = 0.;
      double s0_ = 0.;
      if (do_dump_hist){
        sort(logelong_vec.begin(), logelong_vec.end());
        ofstream logelongfile(histfolder + "/logelong_t" + to_string(t) + ".hist");
        for (vector<array<double, 3>>::const_iterator lit = logelong_vec.begin();
             lit != logelong_vec.end(); ++lit){
          s_ += (*lit)[1];
          s0_ += (*lit)[2];
          logelongfile << (*lit)[0] << " "
                      << (*lit)[1] << " " << s_ << " " << s_/s << " "
                      << (*lit)[2] << " " << s0_ << " " << s0_/s0 << " "
                      << endl;
        }
        logelongfile.close();
      }
      statfile << s << "\t"                   // 15
               << s0 << "\t"                  // 16
               << logelong_wmean << "\t"       // 17
               << logelong_wvar << "\t"        // 18
               << logelong_w0mean << "\t"      // 19
               << logelong_w0var << "\t";      // 20
    }
  }
  if (faces.size() > 0){
    // Sheet method
    double A = 0.;
    double A0 = 0.;
    double logthickness_wmean = 0.;
    double logthickness_w0mean = 0.;
    vector<array<double, 3>> logthickness_vec;
    for (FacesType::const_iterator faceit = faces.begin();
         faceit != faces.end(); ++faceit){
      Uint iedge = faceit->first[0];
      Uint jedge = faceit->first[1];
      // Uint kedge = faceit->first[2];
      double dA0 = faceit->second;
      double dA = area(iedge, jedge, x_rw, y_rw, z_rw, edges);
      double logthickness = log(dA0/dA);
      logthickness_wmean += logthickness*dA;
      logthickness_w0mean += logthickness*dA0;
      logthickness_vec.push_back({logthickness, dA, dA0});
      A += dA;
      A0 += dA0;
    }
    logthickness_wmean /= A;
    logthickness_w0mean /= A0;
    double logthickness_wvar = 0.;
    double logthickness_w0var = 0.;
    for (vector<array<double, 3>>::const_iterator lit = logthickness_vec.begin();
         lit != logthickness_vec.end(); ++lit){
      logthickness_wvar += pow((*lit)[0]-logthickness_wmean, 2)*(*lit)[1];
      logthickness_w0var += pow((*lit)[0]-logthickness_w0mean, 2)*(*lit)[1];
    }
    logthickness_wvar /= A;
    logthickness_w0var /= A0;
    double A_ = 0.;
    double A0_ = 0.;
    if (do_dump_hist){
      sort(logthickness_vec.begin(), logthickness_vec.end());
      ofstream logthicknessfile(histfolder + "/logthickness_t" + to_string(t) + ".hist");
      for (vector<array<double, 3>>::const_iterator lit=logthickness_vec.begin();
           lit != logthickness_vec.end(); ++lit){
        A_ += (*lit)[1];
        A0_ += (*lit)[2];
        logthicknessfile << (*lit)[0] << " "
                         << (*lit)[1] << " " << A_ << " " << A_/A << " "
                         << (*lit)[2] << " " << A0_ << " " << A0_/A0 << " "
                         << endl;
      }
      logthicknessfile.close();
    }
    statfile << A << "\t"                   // 15
             << A0 << "\t"                  // 16
             << logthickness_wmean << "\t"  // 17
             << logthickness_wvar << "\t"   // 18
             << logthickness_w0mean << "\t" // 19
             << logthickness_w0var << "\t"; // 20
  }
  statfile << std::endl;
}

void write_stats_header(ofstream &statfile,
                        const FacesType &faces,
                        const EdgesType &edges){
  statfile << "# t" << "\t"                   //  1
           << "x_mean" << "\t"                //  2
           << "dx2_mean" << "\t"              //  3
           << "y_mean" << "\t"                //  4
           << "dy2_mean" << "\t"              //  5
           << "z_mean" << "\t"                //  6
           << "dz2_mean" << "\t"              //  7
           << "ux_mean" << "\t"               //  8
           << "uy_mean" << "\t"               //  9
           << "uz_mean" << "\t"               // 10
           << "Nrw" << "\t"                   // 11
           << "n_accepted" << "\t"            // 12
           << "n_declined" << "\t";           // 13
  if (edges.size() > 0){
    statfile << "n_too_long" << "\t";         // 14
  }
  if (faces.size() == 0){
    statfile << "s" << "\t"                   // 15
             << "s0" << "\t"                  // 16
             << "logelong_wmean" << "\t"      // 17
             << "logelong_wvar" << "\t"       // 18
             << "logelong_w0mean" << "\t"     // 19
             << "logelong_w0var" << "\t";     // 20
  }
  else {
    statfile << "A" << "\t"                   // 15
             << "A0" << "\t"                  // 16
             << "logthickness_wmean" << "\t"  // 17
             << "logthickness_wvar" << "\t"   // 18
             << "logthickness_w0mean" << "\t" // 19
             << "logthickness_w0var" << "\t"; // 20
  }
  statfile << endl;
}

void print_mesh(const FacesType& faces, const EdgesType& edges,
                const Edge2FacesType& edge2faces,
                double* x_rw, double* y_rw, double* z_rw, const Uint Nrw){
  cout << "==================" << endl;
  cout << "nodes" << endl;
  for (Uint irw=0; irw < Nrw; ++irw){
    cout << irw << " " << x_rw[irw] << " " << y_rw[irw] << " " << z_rw[irw] << endl;
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
               double* x_rw, double* y_rw, double* z_rw, const Uint Nrw){
  ofstream nodesf("mesh.node");
  for (Uint irw=0; irw < Nrw; ++irw){
    nodesf << x_rw[irw] << " " << y_rw[irw] << " " << z_rw[irw] << endl;
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
                 double* x_rw, double* y_rw, double* z_rw,
                 double* ux_rw, double* uy_rw, double* uz_rw,
                 const Uint Nrw){
  for (Uint irw=0; irw < Nrw; ++irw){
    pos_out << irw << " "
            << x_rw[irw] << " "
            << y_rw[irw] << " "
            << z_rw[irw] << " "
            << ux_rw[irw] << " "
            << uy_rw[irw] << " "
            << uz_rw[irw] << std::endl;
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
              double* x_rw, double* y_rw, double* z_rw,
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
      dA[iface] = area(iface, x_rw, y_rw, z_rw, faces, edges);
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
  if (edges.size() > 0){
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
      dl[iedge] = dist(edges[iedge].first[0], edges[iedge].first[1], x_rw, y_rw, z_rw);
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
                     double* x_rw, double* y_rw, double* z_rw,
                     double* ux_rw, double* uy_rw, double* uz_rw,
                     double* rho_rw, double* p_rw, double* c_rw,
                     double* Jux_rw, double* Juy_rw, double* Juz_rw,
                     Uint& Nrw, const bool do_output_all, Interpol &intp,
                     const double U0, const int int_order){
  x_rw[Nrw] = 0.5*(x_rw[inode]+x_rw[jnode]);
  y_rw[Nrw] = 0.5*(y_rw[inode]+y_rw[jnode]);
  z_rw[Nrw] = 0.5*(z_rw[inode]+z_rw[jnode]);

  c_rw[Nrw] = 0.5*(c_rw[inode]+c_rw[jnode]);

  intp.probe(x_rw[Nrw], y_rw[Nrw], z_rw[Nrw]);

  ux_rw[Nrw] = U0*intp.get_ux();
  uy_rw[Nrw] = U0*intp.get_uy();
  uz_rw[Nrw] = U0*intp.get_uz();
  if (do_output_all){
    rho_rw[Nrw] = intp.get_rho();
    p_rw[Nrw] = intp.get_p();
  }
  // Second-order terms
  if (int_order >= 2){
    double U02 = U0*U0;
    Jux_rw[Nrw] = U02*intp.get_Jux();
    Juy_rw[Nrw] = U02*intp.get_Juy();
    Juz_rw[Nrw] = U02*intp.get_Juz();
  }
  ++Nrw;
}

void sheet_refinement(FacesType &faces,
                      EdgesType &edges,
                      Edge2FacesType &edge2faces,
                      double* x_rw, double* y_rw, double* z_rw,
                      double* ux_rw, double* uy_rw, double* uz_rw,
                      double* rho_rw, double* p_rw, double* c_rw,
                      double* Jux_rw, double* Juy_rw, double* Juz_rw,
                      Uint &Nrw, const Uint Nrw_max, const double ds_max,
                      const bool do_output_all, Interpol &intp,
                      const double U0, const int int_order){
  bool changed;
  do {
    changed = false;

    vector<double> ds_;
    for (EdgesType::const_iterator edgeit = edges.begin();
         edgeit != edges.end(); ++edgeit){
      Uint inode = edgeit->first[0];
      Uint jnode = edgeit->first[1];
      double ds = dist(inode, jnode, x_rw, y_rw, z_rw);
      // cout << inode << " " << jnode << " " << ds << endl;
      ds_.push_back(ds);
    }
    vector<size_t> ids_ = argsort_descending(ds_);


    for (vector<size_t>::iterator itedge = ids_.begin();
         itedge != ids_.end(); ++itedge){

      Uint iedge = *itedge;
      if (Nrw >= Nrw_max){
        // No more points can fit
        break;
      }
      double ds = ds_[iedge];
      if (ds < ds_max)
        break;
      changed = true;

      Uint inode = edges[iedge].first[0];
      Uint jnode = edges[iedge].first[1];
      double ds0 = edges[iedge].second;

      Uint new_inode = Nrw;
      // Add point
      append_new_node(inode, jnode,
                      x_rw, y_rw, z_rw, ux_rw, uy_rw, uz_rw,
                      rho_rw, p_rw, c_rw, Jux_rw, Juy_rw, Juz_rw,
                      Nrw, do_output_all, intp, U0, int_order);
      //
      edges[iedge].first[1] = new_inode;
      Uint new_iedge = edges.size();
      edges.push_back({{new_inode, jnode}, ds0/2});
      edge2faces.push_back({});

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
      }
    }
  } while (changed);
}

void strip_refinement(EdgesType &edges,
                      double* x_rw, double* y_rw, double* z_rw,
                      double* ux_rw, double* uy_rw, double* uz_rw,
                      double* rho_rw, double* p_rw, double* c_rw,
                      double* Jux_rw, double* Juy_rw, double* Juz_rw,
                      Uint &Nrw, const Uint Nrw_max, const double ds_max,
                      const bool do_output_all, Interpol &intp,
                      const double U0, const int int_order){
  for (EdgesType::iterator edgeit = edges.begin();
       edgeit != edges.end(); ++edgeit){
    Uint inode = edgeit->first[0];
    Uint jnode = edgeit->first[1];
    double ds0 = edgeit->second;
    double ds = dist(inode, jnode, x_rw, y_rw, z_rw);
    if (ds > ds_max && Nrw < Nrw_max){
      edgeit->first[1] = Nrw;
      edgeit->second = ds0/2;
      edgeit = edges.insert(edgeit+1, {{Nrw, jnode}, ds0/2});
      edgeit -= 2;

      append_new_node(inode, jnode,
                      x_rw, y_rw, z_rw, ux_rw, uy_rw, uz_rw,
                      rho_rw, p_rw, c_rw, Jux_rw, Juy_rw, Juz_rw,
                      Nrw, do_output_all, intp, U0, int_order);
    }
  }
}

void refinement(FacesType &faces,
                EdgesType &edges,
                Edge2FacesType &edge2faces,
                double* x_rw, double* y_rw, double* z_rw,
                double* ux_rw, double* uy_rw, double* uz_rw,
                double* rho_rw, double* p_rw, double* c_rw,
                double* Jux_rw, double* Juy_rw, double* Juz_rw,
                Uint &Nrw, const Uint Nrw_max, const double ds_max,
                const bool do_output_all, Interpol &intp,
                const double U0, const int int_order){
  if (faces.size() > 0){
    sheet_refinement(faces, edges, edge2faces,
                     x_rw, y_rw, z_rw,
                     ux_rw, uy_rw, uz_rw,
                     rho_rw, p_rw, c_rw,
                     Jux_rw, Juy_rw, Juz_rw,
                     Nrw, Nrw_max, ds_max,
                     do_output_all, intp,
                     U0, int_order);
  }
  else {
    strip_refinement(edges,
                     x_rw, y_rw, z_rw,
                     ux_rw, uy_rw, uz_rw,
                     rho_rw, p_rw, c_rw,
                     Jux_rw, Juy_rw, Juz_rw,
                     Nrw, Nrw_max, ds_max,
                     do_output_all, intp,
                     U0, int_order);
  }
}

int main(int argc, char* argv[]){
  // Input parameters
  if (argc < 2) {
    std::cout << "Specify an input timestamps file." << std::endl;
    return 0;
  }
  Parameters prm(argc, argv);
  if (prm.restart_folder != ""){
    prm.parse_file(prm.restart_folder + "/Checkpoints/params.dat");
    prm.parse_cmd(argc, argv);
  }

  string infilename = string(argv[1]);

  Interpol intp(infilename);

  double Dm = prm.Dm;
  double dt = prm.dt;
  Uint Nrw = prm.Nrw;
  double U0 = prm.U0;

  bool refine = prm.refine;
  double ds_max = prm.ds_max;
  Uint Nrw_max = prm.Nrw_max;

  string folder = intp.get_folder();
  string rwfolder = create_folder(folder + "/RandomWalkers/");
  string newfolder;
  if (prm.restart_folder != ""){
    newfolder = prm.folder;
  }
  else {
    ostringstream ss_Dm, ss_dt, ss_Nrw;
    ss_Dm << std::scientific << std::setprecision(7) << Dm;
    ss_dt << std::scientific << std::setprecision(7) << dt;
    ss_Nrw << Nrw;
    newfolder = create_folder(rwfolder +
                              "/Dm" + ss_Dm.str() + // "_U" + to_string(prm.U0) +
                              "_dt" + ss_dt.str() +
                              "_Nrw" + ss_Nrw.str() + "/");
  }
  string posfolder = create_folder(newfolder + "Positions/");
  string checkpointsfolder = create_folder(newfolder + "Checkpoints/");
  string histfolder = create_folder(newfolder + "Histograms/");
  prm.folder = newfolder;

  prm.print();

  random_device rd;
  mt19937 gen(rd());

  double Lx = intp.get_Lx();
  double Ly = intp.get_Ly();
  double Lz = intp.get_Lz();
  prm.Lx = Lx;
  prm.Ly = Ly;
  prm.Lz = Lz;
  prm.nx = intp.get_nx();
  prm.ny = intp.get_ny();
  prm.nz = intp.get_nz();

  double U02 = U0*U0;

  double t0 = max(intp.get_t_min(), prm.t0);
  double T = min(intp.get_t_max(), prm.T);
  prm.t0 = t0;
  prm.T = T;

  if (prm.interpolation_test > 0){
    int n = 0;
    double x, y, z;
    double ux, uy, uz;
    double rho, p;
    double divu, vortz;

    intp.update(t0);

    ofstream nodalfile(newfolder + "/nodal_values.dat");
    for (Uint ix=0; ix<intp.get_nx(); ++ix){
      for (Uint iy=0; iy<intp.get_ny(); ++iy){
        for (Uint iz=0; iz<intp.get_nz(); ++iz){
          bool inside = intp.get_nodal_inside(ix, iy, iz);
          ux = intp.get_nodal_ux(ix, iy, iz);
          uy = intp.get_nodal_uy(ix, iy, iz);
          uz = intp.get_nodal_uz(ix, iy, iz);
          nodalfile << ix << " " << iy << " " << iz << " " << inside << " "
                    << ux << " " << uy << " " << uz << endl;
        }
      }
    }
    nodalfile.close();

    uniform_real_distribution<> uni_dist_x(0, Lx);
    uniform_real_distribution<> uni_dist_y(0, Ly);
    uniform_real_distribution<> uni_dist_z(0, Lz);

    ofstream ofile(newfolder + "/interpolation.dat");
    while (n < prm.interpolation_test){
      x = uni_dist_x(gen);
      y = uni_dist_y(gen);
      z = uni_dist_z(gen);

      intp.probe(x, y, z);
      if (intp.inside_domain()){
        ux = intp.get_ux();
        uy = intp.get_uy();
        uz = intp.get_uz();
        rho = intp.get_rho();
        p = intp.get_p();
        divu = intp.get_divu();
        vortz = intp.get_vortz();
        ofile << x  << " " << y  << " " << z  << " "
              << ux << " " << uy << " " << uz << " "
              << rho << " " << p << " " << divu << " "
              << vortz << endl;
      }
      ++n;
    }
    ofile.close();
  }

  normal_distribution<double> rnd_normal(0.0, 1.0);

  double* x_rw = new double[Nrw_max];
  double* y_rw = new double[Nrw_max];
  double* z_rw = new double[Nrw_max];

  double* c_rw = new double[Nrw_max];
  double* e_rw = new double[Nrw_max];

  double* ux_rw = new double[Nrw_max];
  double* uy_rw = new double[Nrw_max];
  double* uz_rw = new double[Nrw_max];
  double* rho_rw = new double[Nrw_max];
  double* p_rw = new double[Nrw_max];

  // Second-order terms
  double* Jux_rw = new double[Nrw_max];
  double* Juy_rw = new double[Nrw_max];
  double* Juz_rw = new double[Nrw_max];
  if (prm.int_order > 2){
    cout << "No support for such high temporal integration order." << endl;
    exit(0);
  }

  vector<array<double, 3>> pos_init;
  EdgesType edges;
  FacesType faces;

  if (prm.restart_folder != ""){
    string posfile = prm.restart_folder + "/Checkpoints/positions.pos";
    load_positions(posfile, pos_init, Nrw);
    string facefile = prm.restart_folder + "/Checkpoints/faces.face";
    load_faces(facefile, faces);
    string edgefile = prm.restart_folder + "/Checkpoints/edges.edge";
    load_edges(edgefile, edges);
    string colfile = prm.restart_folder + "/Checkpoints/colors.col";
    load_colors(colfile, c_rw, Nrw);
  }
  else {
    pos_init = initial_positions(prm.init_mode,
                                 prm.init_weight,
                                 Nrw,
                                 prm.x0, prm.y0, prm.z0,
                                 prm.La, prm.Lb,
                                 prm.ds_max,
                                 intp, gen, edges, faces);
  }
  // Compute edge2faces map
  Edge2FacesType edge2faces;
  for (EdgesType::const_iterator edgeit = edges.begin();
       edgeit != edges.end(); ++edgeit){
    edge2faces.push_back({});
  }
  for (Uint iface=0; iface < faces.size(); ++iface){
    for (Uint i=0; i < 3; ++i){
      edge2faces[faces[iface].first[i]].push_back(iface);
    }
  }

  for (Uint irw=0; irw < Nrw; ++irw){
    // Assign initial position
    x_rw[irw] = pos_init[irw][0];
    y_rw[irw] = pos_init[irw][1];
    z_rw[irw] = pos_init[irw][2];

    if (prm.restart_folder == ""){
      c_rw[irw] = double(irw)/(Nrw-1);
    }
    intp.update(t0);
    intp.probe(x_rw[irw], y_rw[irw], z_rw[irw]);
    ux_rw[irw] = U0*intp.get_ux();
    uy_rw[irw] = U0*intp.get_uy();
    uz_rw[irw] = U0*intp.get_uz();

    rho_rw[irw] = intp.get_rho();
    p_rw[irw] = intp.get_p();

    // Second-order terms
    if (prm.int_order >= 2){
      Jux_rw[irw] = U02*intp.get_Jux();
      Juy_rw[irw] = U02*intp.get_Juy();
      Juz_rw[irw] = U02*intp.get_Juz();
    }
  }
  // Initial refinement
  if (refine){
    cout << "initial refinement" << endl;

    bool do_output_all = true;
    refinement(faces, edges, edge2faces,
               x_rw, y_rw, z_rw,
               ux_rw, uy_rw, uz_rw,
               rho_rw, p_rw, c_rw,
               Jux_rw, Juy_rw, Juz_rw,
               Nrw, Nrw_max, ds_max,
               do_output_all, intp,
               U0, prm.int_order);

  }

  double eta_x, eta_y, eta_z;
  double dx_rw, dy_rw, dz_rw;
  int it = 0;
  double dt2 = dt*dt;

  double sqrt2Dmdt = sqrt(2*Dm*dt);
  if (prm.verbose){
    print_param("sqrt(2*Dm*dt)", sqrt2Dmdt);
    print_param("U*dt         ", U0*dt);
  }

  double t = t0;
  if (prm.restart_folder != ""){
    t = prm.t;
  }

  prm.dump(newfolder, t);

  unsigned long int n_accepted = prm.n_accepted;
  unsigned long int n_declined = prm.n_declined;


  H5File* h5f = new H5File(newfolder + "/data_from_t" + to_string(t) + ".h5", H5F_ACC_TRUNC);

  int int_stat_intv = int(prm.stat_intv/dt);
  int int_dump_intv = int(prm.dump_intv/dt);
  int int_checkpoint_intv = int(prm.checkpoint_intv/dt);
  int int_chunk_intv = int_dump_intv*prm.dump_chunk_size;
  int int_refine_intv = int(prm.refine_intv/dt);
  int int_hist_intv = int_stat_intv*prm.hist_chunk_size;

  string write_mode = prm.write_mode;

  ofstream statfile(newfolder + "/tdata_from_t" + to_string(t) + ".dat");
  write_stats_header(statfile, faces, edges);
  ofstream declinedfile(newfolder + "/declinedpos_from_t" + to_string(t) + ".dat");
  while (t <= T){
    intp.update(t);
    // Statistics
    if (it % int_stat_intv == 0){
      cout << "Time = " << t << endl;
      bool do_dump_hist = (int_hist_intv > 0 && it % int_hist_intv == 0);
      write_stats(statfile, t,
                  x_rw, y_rw, z_rw,
                  ux_rw, uy_rw, uz_rw,
                  Nrw, faces, edges,
                  ds_max,
                  do_dump_hist,
                  histfolder,
                  n_accepted, n_declined);
    }
    // Checkpoint
    if (it % int_checkpoint_intv == 0){
       prm.t = t;
       prm.n_accepted = n_accepted;
       prm.n_declined = n_declined;
       prm.Nrw = Nrw;
       prm.dump(checkpointsfolder);
       dump_positions(checkpointsfolder + "/positions.pos", x_rw, y_rw, z_rw, Nrw);
       dump_faces(checkpointsfolder + "/faces.face", faces);
       dump_edges(checkpointsfolder + "/edges.edge", edges);
       dump_colors(checkpointsfolder + "/colors.col", c_rw, Nrw);
    }
    // Refinement
    if (refine && it % int_refine_intv == 0){
      bool do_output_all = it % int_dump_intv == 0;
      refinement(faces, edges, edge2faces,
                 x_rw, y_rw, z_rw,
                 ux_rw, uy_rw, uz_rw,
                 rho_rw, p_rw, c_rw,
                 Jux_rw, Juy_rw, Juz_rw,
                 Nrw, Nrw_max, ds_max,
                 do_output_all, intp,
                 U0, prm.int_order);
    }
    // Dump detailed data
    if (it % int_dump_intv == 0){
      if (write_mode == "text"){
        string posfile = posfolder + "xy_t" + to_string(t) + ".pos";
        ofstream pos_out(posfile);
        posdata2txt(pos_out,
                    x_rw, y_rw, z_rw,
                    ux_rw, uy_rw, uz_rw, Nrw);
        pos_out.close();
      }
      else if (write_mode == "hdf5"){
        // Clear file if it exists, otherwise create
        if (it % int_chunk_intv == 0 && it > 0){
          h5f->close();
          h5f = new H5File(newfolder + "/data_from_t" + to_string(t) + ".h5", H5F_ACC_TRUNC);
        }
        string groupname = to_string(t);
        h5f->createGroup(groupname + "/");

        for (EdgesType::const_iterator edgeit = edges.begin();
             edgeit != edges.end(); ++edgeit){
          int inode = edgeit->first[0];
          int jnode = edgeit->first[1];
          double ds0 = edgeit->second;
          double ds = dist(inode, jnode, x_rw, y_rw, z_rw);
          if (e_rw[inode] <= 0.) e_rw[inode] = ds/ds0;
          else e_rw[inode] = 0.5*(e_rw[inode] + ds/ds0);
          if (e_rw[jnode] <= 0.) e_rw[jnode] = ds/ds0;
          else e_rw[jnode] = 0.5*(e_rw[jnode] + ds/ds0);
        }

        vector2hdf5(h5f, groupname + "/points", x_rw, y_rw, z_rw, Nrw);
        vector2hdf5(h5f, groupname + "/u", ux_rw, uy_rw, uz_rw, Nrw);
        scalar2hdf5(h5f, groupname + "/rho", rho_rw, Nrw);
        scalar2hdf5(h5f, groupname + "/p", p_rw, Nrw);
        scalar2hdf5(h5f, groupname + "/c", c_rw, Nrw);
        scalar2hdf5(h5f, groupname + "/e", e_rw, Nrw);
        mesh2hdf(h5f, groupname, x_rw, y_rw, z_rw, faces, edges);
      }
    }
    for (Uint irw=0; irw < Nrw; ++irw){
      dx_rw = ux_rw[irw]*dt;
      dy_rw = uy_rw[irw]*dt;
      dz_rw = uz_rw[irw]*dt;

      // Set elongation
      e_rw[irw] = 0.;

      // Second-order terms
      if (prm.int_order >= 2){
        dx_rw += 0.5*Jux_rw[irw]*dt2;
        dy_rw += 0.5*Juy_rw[irw]*dt2;
        dz_rw += 0.5*Juz_rw[irw]*dt2;
      }
      if (Dm > 0.0){
        eta_x = rnd_normal(gen);
        eta_y = rnd_normal(gen);
        eta_z = rnd_normal(gen);
        dx_rw += sqrt2Dmdt*eta_x;
        dy_rw += sqrt2Dmdt*eta_y;
        dz_rw += sqrt2Dmdt*eta_z;
      }
      intp.probe(x_rw[irw]+dx_rw,
                 y_rw[irw]+dy_rw,
                 z_rw[irw]+dz_rw);
      if (intp.inside_domain()){
        x_rw[irw] += dx_rw;
        y_rw[irw] += dy_rw;
        z_rw[irw] += dz_rw;

        ux_rw[irw] = U0*intp.get_ux();
        uy_rw[irw] = U0*intp.get_uy();
        uz_rw[irw] = U0*intp.get_uz();

        if ((it+1) % int_dump_intv == 0){
          rho_rw[irw] = intp.get_rho();
          p_rw[irw] = intp.get_p();
        }

        // Second-order terms
        if (prm.int_order >= 2){
          Jux_rw[irw] = U02*intp.get_Jux();
          Juy_rw[irw] = U02*intp.get_Juy();
          Juz_rw[irw] = U02*intp.get_Juz();
        }
        n_accepted++;
      }
      else {
        n_declined++;
        declinedfile << t << " "
                     << x_rw[irw]+dx_rw << " "
                     << y_rw[irw]+dy_rw << " "
                     << z_rw[irw]+dz_rw << std::endl;
      }
    }
    t += dt;
    it += 1;
  }

  // Final checkpoint
  prm.t = t;
  prm.n_accepted = n_accepted;
  prm.n_declined = n_declined;
  prm.Nrw = Nrw;
  prm.dump(checkpointsfolder);
  dump_positions(checkpointsfolder + "/positions.pos",
                 x_rw, y_rw, z_rw, Nrw);
  dump_faces(checkpointsfolder + "/faces.face", faces);
  dump_edges(checkpointsfolder + "/edges.edge", edges);
  dump_colors(checkpointsfolder + "/colors.col", c_rw, Nrw);

  // Close files
  statfile.close();
  declinedfile.close();
  if (write_mode == "hdf5"){
    h5f->close();
  }

  return 0;
}
