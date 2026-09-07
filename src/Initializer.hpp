#ifndef __DISTRIBUTE_HPP
#define __DISTRIBUTE_HPP

#include <iomanip>
#include <vector>
#include <random>
#include <algorithm>
#include "typedefs.hpp"
#include "Interpol.hpp"
#include "utils.hpp"
#include "Params.hpp"
#include "ParticleSet.hpp"
#include "mesh.hpp"

// true if init_mode starts with one of the given kinds
inline bool init_mode_is(const partrac::Params& prm, const std::vector<std::string>& kinds){
  const std::string kind = split_string(prm.get<std::string>("init_mode"), "_")[0];
  return std::find(kinds.begin(), kinds.end(), kind) != kinds.end();
}

// The initializers index the init_mode tokens without bounds checks
inline bool init_mode_shape_ok(const std::string& init_mode){
  if (init_mode.rfind("from", 0) == 0)
    return init_mode.rfind("from_file:", 0) == 0 && init_mode.size() > 10;
  const std::vector<std::string> key = split_string(init_mode, "_");
  if (key[0] == "point") return true;
  if (key[0] == "randomgaussianstrip") return key.size() >= 3;
  return key.size() >= 2;
}

// Parameters read by set_initial_state and the initializers it dispatches to.
// La, Lb, ds_init and init_weight are only read by some init_modes.
inline void add_initializer_params(partrac::Schema& s){
  s.require<std::string>("init_mode", "initial distribution");
  s.require<Uint>("Nrw", "number of particles");
  // Nrw stays the request; these record what happened
  s.runtime<Uint>("Nrw_init", 0, "particles the initializer placed");
  s.runtime<Uint>("Nrw_current", 0, "particles in the set when this was written");
  s.require<Uint>("Nrw_max", "max number of particles");
  s.require<double>("ds_max", "max edge length");
  s.require<double>("ds_min", "min edge length");
  s.require_if<double>("La",
                       [](const partrac::Params& p){
                         return init_mode_is(p, {"strip", "sheet", "ellipsoid",
                                                 "randomgaussianstrip",
                                                 "randomgaussiancircle"});
                       },
                       "init_mode is a strip, sheet, ellipsoid or gaussian",
                       "principal extent");
  s.require_if<double>("Lb",
                       [](const partrac::Params& p){
                         return init_mode_is(p, {"sheet", "ellipsoid",
                                                 "randomgaussianstrip",
                                                 "randomgaussiancircle"});
                       },
                       "init_mode is a sheet, ellipsoid or gaussian",
                       "second extent");
  s.require_if<double>("ds_init",
                       [](const partrac::Params& p){
                         return init_mode_is(p, {"sheet", "pair", "pairs", "points"});
                       },
                       "init_mode is a sheet, pair or points distribution",
                       "initial edge length");
  s.require_if<std::string>("init_weight",
                            [](const partrac::Params& p){
                              return init_mode_is(p, {"points"});
                            },
                            "init_mode is a points distribution",
                            "sampling weight");
  s.opt<double>("t0", 0.0, "start time");
  s.opt<double>("x0", 0.0, "initial position");
  s.opt<double>("y0", 0.0, "initial position");
  s.opt<double>("z0", 0.0, "initial position");
  s.opt<bool>("inject", false, "inject new particles");
  s.opt<bool>("clear_initial_edges", false, "drop the initial edges");
  s.check([](const partrac::Params& p){
            return init_mode_shape_ok(p.get<std::string>("init_mode"));
          },
          "init_mode is missing a direction: most modes need one, as in"
          " uniform_x; randomgaussianstrip needs two, as in"
          " randomgaussianstrip_x_y; from_file needs a path, as in"
          " from_file:positions.h5");
}

// TODO: Massive cleanup!
struct less_than_op {
  inline bool operator() (const Vector3d &a, const Vector3d &b){
    return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]) || (a[0] == b[0] && a[1] == b[1] && a[2] < b[2]);
  }
};

class Initializer {
public:
  Initializer(std::shared_ptr<Interpol> intp, partrac::Params& prm) : intp(intp), prm(prm) {
    x0 = {prm.get<double>("x0"), prm.get<double>("y0"), prm.get<double>("z0")};
    x_min = intp->get_x_min();
    x_max = intp->get_x_max();
    L = x_max - x_min;
    inject = prm.get<bool>("inject");
    clear_initial_edges = prm.get<bool>("clear_initial_edges");
  };
  ~Initializer() { nodes.clear(); edges.clear(); faces.clear(); };
  /*std::vector<Vector3d>::const_iterator node_begin() const { return nodes.begin(); };
  std::vector<Vector3d>::const_iterator node_end() const { return nodes.end(); };
  EdgesType::const_iterator edge_begin() const { return edges.begin(); };
  EdgesType::const_iterator edge_end() const { return edges.end(); };
  FacesType::const_iterator face_begin() const { return faces.begin(); };
  FacesType::const_iterator face_end() const { return faces.end(); };*/
  std::vector<Vector3d> nodes;
  EdgesType edges;
  FacesType faces;
  bool inject;
  bool clear_initial_edges;
protected:
  std::shared_ptr<Interpol> intp;
  partrac::Params& prm;
  Vector3d x0;
  Vector3d x_min;
  Vector3d x_max;
  Vector3d L;
};

class PointInitializer : public Initializer {
public:
  PointInitializer(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, partrac::Params& prm) : Initializer(intp, prm) {
    const bool inside = intp->locate(x0);
    for (Uint irw=0; irw < prm.get<Uint>("Nrw"); ++irw){
      if (inside){
        nodes.push_back(x0);
      }
    }
    edges.clear();
    faces.clear();
  };
};

class UniformInitializer : public Initializer {
public:
  UniformInitializer(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, partrac::Params& prm) : Initializer(intp, prm) {
    Vector3d x_a = x0;
    Vector3d x_b = x0;
    if (key[1] == "x"){
      x_a[0] = x_min[0];
      x_b[0] = x_max[0];
    }
    else if (key[1] == "y"){
      x_a[1] = x_min[1];
      x_b[1] = x_max[1];
    }
    else if (key[1] == "z"){
      x_a[2] = x_min[2];
      x_b[2] = x_max[2];
    }
    else {
      std::cout << "Unrecognized initialization..." << std::endl;
      exit(1);
    }
    Vector3d Dx = (x_b - x_a) / (prm.get<Uint>("Nrw")-1);
    for (Uint irw=0; irw < prm.get<Uint>("Nrw"); ++irw){
      Vector3d x = x_a + Dx * irw;
      if (intp->locate(x)){
        nodes.push_back(x);
      }
    }
    for (Uint irw=0; irw < nodes.size()-1; ++irw){
      if ((nodes[irw] - nodes[irw+1]).norm() < 1.5*Dx.norm()){
        edges.push_back({{irw, irw+1}, dist(nodes[irw], nodes[irw+1])});
      }
    }
  };
};

class StripInitializer : public Initializer {
public:
  StripInitializer(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, partrac::Params& prm) : Initializer(intp, prm) {
    double La = prm.get<double>("La");
    // double Lb = prm.Lb;
    Vector3d n(0., 0., 0.);
    if (key[1] == "x"){
      n[0] = 1.0;
    }
    if (key[1] == "y"){
      n[1] = 1.0;
    }
    if (key[1] == "z"){
      n[2] = 1.0;
    }

    Vector3d x00 = x0;
    x00[0] += - La/2*n[0];
    x00[1] += - La/2*n[1];
    x00[2] += - La/2*n[2];

    Vector3d x01 = x0;
    x01[0] += La/2*n[0];
    x01[1] += La/2*n[1];
    x01[2] += La/2*n[2];

    Uint irw = 0;
    bool this_inside = false;
    bool prev_inside = false;

    for (Uint i=0; i < prm.get<Uint>("Nrw"); ++i){
      double alpha = float(i)/(prm.get<Uint>("Nrw")-1);
      Vector3d x0i = alpha * x00 + (1.-alpha) * x01;
      // check if inside domain
      this_inside = intp->locate(x0i);
      if (this_inside){
        nodes.push_back(x0i);
        if (prev_inside)
          edges.push_back({{irw-1, irw}, dist(nodes[irw], x0i)});
        ++irw;
      }
      prev_inside = this_inside;
    }
    if (irw == 0) {
      std::cout << "Strip not inside domain" << std::endl;
      exit(1);
    }
  };
};

class SheetInitializer : public Initializer {
public:
  SheetInitializer(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, partrac::Params& prm) : Initializer(intp, prm) {
    double La = prm.get<double>("La");
    double Lb = prm.get<double>("Lb");
    Vector3d n(0., 0., 0.);
    Vector3d ta(0., 0., 0.);
    Vector3d tb(0., 0., 0.);
    if (key[1] == "xy"){
      n[2] = 1.0;
      ta[0] = 1.0;
      tb[1] = 1.0;
    }
    if (key[1] == "xz"){
      n[1] = 1.0;
      ta[0] = 1.0;
      tb[2] = 1.0;
    }
    if (key[1] == "yz"){
      n[0] = 1.0;
      ta[1] = 1.0;
      tb[2] = 1.0;
    }

    Vector3d x00 = x0;
    x00[0] += - La/2*ta[0] - Lb/2*tb[0];
    x00[1] += - La/2*ta[1] - Lb/2*tb[1];
    x00[2] += - La/2*ta[2] - Lb/2*tb[2];

    Vector3d x01 = x0;
    x01[0] += La/2*ta[0] + Lb/2*tb[0];
    x01[1] += La/2*ta[1] - Lb/2*tb[1];
    x01[2] += - La/2*ta[2] - Lb/2*tb[2];

    Vector3d x10 = x0;
    x10[0] += La/2*ta[0] + Lb/2*tb[0];
    x10[1] += La/2*ta[1] + Lb/2*tb[1];
    x10[2] += La/2*ta[2] + Lb/2*tb[2];

    Vector3d x11 = x0;
    x11[0] += - La/2*ta[0] - Lb/2*tb[0];
    x11[1] += - La/2*ta[1] + Lb/2*tb[1];
    x11[2] += - La/2*ta[2] + Lb/2*tb[2];


    nodes.push_back(x00);
    nodes.push_back(x01);
    nodes.push_back(x10);
    nodes.push_back(x11);

    edges.push_back({{0, 1}, dist(nodes[0], nodes[1])});
    edges.push_back({{0, 2}, dist(nodes[0], nodes[2])});
    edges.push_back({{1, 2}, dist(nodes[1], nodes[2])});
    edges.push_back({{2, 3}, dist(nodes[2], nodes[3])});
    edges.push_back({{3, 0}, dist(nodes[3], nodes[4])});

    faces.push_back({{0, 2, 1}, La*Lb/2});
    faces.push_back({{1, 3, 4}, La*Lb/2});

    ParticleSet pset_loc(intp, prm.get<Uint>("Nrw_max"));
    pset_loc.add(nodes, 0);

    Edge2FacesType edge2faces_loc;
    Node2EdgesType node2edges_loc;

    compute_edge2faces(edge2faces_loc, faces, edges);
    compute_node2edges(node2edges_loc, edges, pset_loc.N());

    NodesListType nodes_inlet_dummy;
    EdgesListType edges_inlet_dummy;

    Uint n_add = 0;
    Uint n_rem = 0;
    do {
      n_add = sheet_refinement(faces, edges, edge2faces_loc, node2edges_loc, edges_inlet_dummy,
                                    pset_loc, prm.get<double>("ds_init"), 0.0, false, false);

      std::cout << "Added " << n_add << " edges." << std::endl;
    } while (n_add > 0);

    for (Uint iedge=0; iedge<edges.size(); ++iedge){
      Uint inode = edges[iedge].first[0];
      Uint jnode = edges[iedge].first[1];
      edges[iedge].second = pset_loc.dist(inode, jnode);
    }
    for (Uint iface=0; iface<faces.size(); ++iface){
      Uint iedge = faces[iface].first[0];
      Uint jedge = faces[iface].first[1];
      faces[iface].second = pset_loc.triangle_area(iedge, jedge, edges);
    }

    compute_edge2faces(edge2faces_loc, faces, edges);
    compute_node2edges(node2edges_loc, edges, pset_loc.N());

    std::set<Uint> edges_to_remove;

    for (Uint irw=0; irw < pset_loc.N(); ++irw){
      int cell_id = pset_loc.get_cell_id(irw);
      bool inside = intp->locate(pset_loc.x(irw), 0., cell_id);
      if (!inside){
        edges_to_remove.insert(node2edges_loc[irw].begin(), node2edges_loc[irw].end());
      }
    }

    std::cout << "Marking faces." << std::endl;

    std::vector<bool> edge_isactive(edges.size(), true);
    std::vector<bool> face_isactive(faces.size(), true);

    for (auto & jedge : edges_to_remove ){
      edge_isactive[jedge] = false;
      for ( auto & jface : edge2faces_loc[jedge] ){
          face_isactive[jface] = false;
      }
    }

    

    std::cout << "Removing faces." << std::endl;

    remove_faces(faces, face_isactive);
    remove_edges(faces, edges, edge_isactive, edges_inlet_dummy);
  
    // remove_unused_edges(faces, edges, edges_inlet_dummy);
    remove_unused_nodes(edges, nodes_inlet_dummy, pset_loc);
    compute_edge2faces(edge2faces_loc, faces, edges);
    compute_node2edges(node2edges_loc, edges, pset_loc.N());

    std::cout << "Removed faces in solid." << std::endl;

    std::vector<Uint> nums = {}; // 7 ? but doesn't work properly


    for ( auto & num : nums ){
      for (Uint inode=0; inode < pset_loc.N(); ++inode){
        if (node2edges_loc[inode].size() == num){
          std::vector<Uint> free_edges;
          for (auto & iedge : node2edges_loc[inode]){
            if (edge2faces_loc[iedge].size() == 1){
              free_edges.push_back(iedge);
            }
          }

          if ( free_edges.size() == 2 ){
            std::vector<Uint> unique_nodes;
            std::set_symmetric_difference(
              edges[free_edges[0]].first.begin(), edges[free_edges[0]].first.end(),
              edges[free_edges[1]].first.begin(), edges[free_edges[1]].first.end(),
              back_inserter(unique_nodes)
            );

            Uint iedge = edges.size();
            edges.push_back({{unique_nodes[0], unique_nodes[1]}, dist(nodes[unique_nodes[0]], nodes[unique_nodes[1]])});
            
                    faces.push_back({{iedge, free_edges[0], free_edges[1]}, pset_loc.triangle_area(free_edges[0], free_edges[1], edges)});
          }
        }
      }

      compute_edge2faces(edge2faces_loc, faces, edges);
      compute_node2edges(node2edges_loc, edges, pset_loc.N());
    }

    //n_rem = sheet_coarsening(faces, edges, edge2faces_loc, node2edges_loc, edges_inlet_dummy, nodes_inlet_dummy,
    //                         pset_loc, prm.ds_max, 0.0);

    int attempt = 0;
    int max_attempts = 100;
    do {
      n_rem = sheet_coarsening(faces, edges, edge2faces_loc, node2edges_loc, edges_inlet_dummy, nodes_inlet_dummy,
                               pset_loc, attempt == 0 ? prm.get<double>("ds_min") : prm.get<double>("ds_min"), 0.0);

      n_add = sheet_refinement(faces, edges, edge2faces_loc, node2edges_loc, edges_inlet_dummy,
                                    pset_loc, prm.get<double>("ds_max"), 0.0, false, true);

      std::cout << "Added " << n_add << " and removed " << n_rem << " edges." << std::endl;

      ++attempt;
    } while ( (n_add > 0 || n_rem > 0) && attempt < max_attempts );

    nodes.clear();
    for (Uint irw=0; irw<pset_loc.N(); ++irw){
      nodes.push_back(pset_loc.x(irw));
    }

  };
};

class EllipsoidInitializer : public Initializer {
public:
  EllipsoidInitializer(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, partrac::Params& prm) : Initializer(intp, prm) {
    double La = prm.get<double>("La");
    double Lb = prm.get<double>("Lb");
    double lx2 = Lb*Lb;
    double ly2 = Lb*Lb;
    double lz2 = Lb*Lb;
    Vector3d n(0., 0., 0.);
    if (key[1] == "xy"){
      n[2] = 1.0;
      lz2 = La * La;
    }
    if (key[1] == "xz"){
      n[1] = 1.0;
      ly2 = La * La;
    }
    if (key[1] == "yz"){
      n[0] = 1.0;
      lx2 = La * La;
    }

    //FacesType faces_loc;
    //EdgesType edges_loc;
    Edge2FacesType edge2faces_loc;
    Node2EdgesType node2edges_loc;

    double R = sqrt(La*Lb);

    Vector3d x_c = {prm.get<double>("x0"), prm.get<double>("y0"), prm.get<double>("z0")};

    Vector3d x_0 = {prm.get<double>("x0") - R/sqrt(2.), prm.get<double>("y0") - R/sqrt(6.0),   prm.get<double>("z0") - R/sqrt(3.0)/2};
    Vector3d x_1 = {prm.get<double>("x0") + R/sqrt(2.), prm.get<double>("y0") - R/sqrt(6.0),   prm.get<double>("z0") - R/sqrt(3.0)/2};
    Vector3d x_2 = {prm.get<double>("x0"),              prm.get<double>("y0") + R*sqrt(2./3.), prm.get<double>("z0") - R/sqrt(3.0)/2};
    Vector3d x_3 = {prm.get<double>("x0"),              prm.get<double>("y0"),                 prm.get<double>("z0") + R*sqrt(3.0)/2};

    std::cout << x_0.norm() << std::endl;
    std::cout << x_1.norm() << std::endl;
    std::cout << x_2.norm() << std::endl;
    std::cout << x_3.norm() << std::endl;
    std::cout << (x_1-x_0).norm() << std::endl;
    std::cout << (x_2-x_0).norm() << std::endl;
    std::cout << (x_3-x_0).norm() << std::endl;
    std::cout << (x_2-x_1).norm() << std::endl;
    std::cout << (x_3-x_1).norm() << std::endl;
    std::cout << (x_3-x_2).norm() << std::endl;


    double t0 = 0.;  // Not used in practice
    int cell_id = -1;

    bool inside_0 = intp->locate(x_0, t0, cell_id);
    bool inside_1 = intp->locate(x_1, t0, cell_id);
    bool inside_2 = intp->locate(x_2, t0, cell_id);
    bool inside_3 = intp->locate(x_3, t0, cell_id);

    if (inside_0 && inside_1 && inside_2 && inside_3){
      std::cout << "Ellipsoid inside domain." << std::endl;
    }
    else {
      std::cout << "Ellipsoid not inside domain" << std::endl;
      
      

      exit(1);
    }

    std::vector<Vector3d> nodes_loc;
    nodes_loc.push_back(x_0);
    nodes_loc.push_back(x_1);
    nodes_loc.push_back(x_2);
    nodes_loc.push_back(x_3);

    ParticleSet pset_loc(intp, prm.get<Uint>("Nrw_max"));
    pset_loc.add(nodes_loc, 0);

    edges.push_back({{0, 1}, dist(nodes_loc[0], nodes_loc[1])});
    edges.push_back({{1, 2}, dist(nodes_loc[1], nodes_loc[2])});
    edges.push_back({{2, 0}, dist(nodes_loc[2], nodes_loc[0])});
    edges.push_back({{1, 3}, dist(nodes_loc[1], nodes_loc[3])});
    edges.push_back({{2, 3}, dist(nodes_loc[2], nodes_loc[3])});
    edges.push_back({{0, 3}, dist(nodes_loc[0], nodes_loc[3])});

    faces.push_back({{0, 1, 2}, 1.});
    faces.push_back({{0, 3, 5}, 1.});
    faces.push_back({{1, 4, 3}, 1.});
    faces.push_back({{2, 5, 4}, 1.});

    NodesListType nodes_inlet_dummy;
    EdgesListType edges_inlet_dummy;
    compute_edge2faces(edge2faces_loc, faces, edges);
    compute_node2edges(node2edges_loc, edges, pset_loc.N());

    Uint n_add, n_rem;
    do {
      n_add = sheet_refinement(faces, edges, edge2faces_loc, node2edges_loc, edges_inlet_dummy,
                                    pset_loc, prm.get<double>("ds_max"), 0.0, false);
      for (Uint irw=0; irw<pset_loc.N(); ++irw){
        Vector3d x = pset_loc.x(irw);
        Vector3d nn = (x - x_c)/ (x - x_c).norm();
        double rad = 1./sqrt(nn[0]*nn[0]/lx2 + nn[1]*nn[1]/ly2 + nn[2]*nn[2]/lz2);
        pset_loc.set_x(irw, x_c + rad * nn);
      }
      n_rem = sheet_coarsening(faces, edges, edge2faces_loc, node2edges_loc, edges_inlet_dummy, nodes_inlet_dummy,
                               pset_loc, prm.get<double>("ds_min"), 0.0);

      std::cout << "Added " << n_add << " and removed " << n_rem << " edges." << std::endl;
    } while (n_add > 0 || n_rem > 0);

    for (Uint iedge=0; iedge<edges.size(); ++iedge){
      Uint inode = edges[iedge].first[0];
      Uint jnode = edges[iedge].first[1];
      edges[iedge].second = pset_loc.dist(inode, jnode);
    }
    for (Uint iface=0; iface<faces.size(); ++iface){
      Uint iedge = faces[iface].first[0];
      Uint jedge = faces[iface].first[1];
      faces[iface].second = pset_loc.triangle_area(iedge, jedge, edges);
    }
    for (Uint irw=0; irw<pset_loc.N(); ++irw){
      nodes.push_back(pset_loc.x(irw));
    }
  };
};

class RandomPairsInitializer : public Initializer {
protected:
  std::mt19937 &gen;
public:
  RandomPairsInitializer(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, partrac::Params& prm, std::mt19937 &gen) : Initializer(intp, prm), gen(gen) {
    std::uniform_real_distribution<> uni_dist_x(x_min[0], x_max[0]);
    std::uniform_real_distribution<> uni_dist_y(x_min[1], x_max[1]);
    std::uniform_real_distribution<> uni_dist_z(x_min[2], x_max[2]);
    std::normal_distribution<double> rnd_normal(0.0, 1.0);

    Uint Npairs = (key[0] == "pair") ? 1 : prm.get<Uint>("Nrw")/2;

    std::cout << "Npairs = " << Npairs << std::endl;

    Vector3d x0_ = x0;
    Uint ipair=0;
    while (ipair < Npairs){
      if (key[0] == "pairs" && key.size() == 3){
        if (contains(key[2], "x")){
          x0_[0] = uni_dist_x(gen);
        }
        if (contains(key[2], "y")){
          x0_[1] = uni_dist_y(gen);
        }
        if (contains(key[2], "z")){
          x0_[2] = uni_dist_z(gen);
        }
      }
      Vector3d dx(0., 0., 0.);
      if (!contains(key[1], "x")){
        dx[0] = 0.;
      }
      else {
        dx[0] = rnd_normal(gen);
      }
      if (!contains(key[1], "y")){
        dx[1] = 0.;
      }
      else {
        dx[1] = rnd_normal(gen);
      }
      if (!contains(key[1], "z")){
        dx[2] = 0.;
      }
      else {
        dx[2] = rnd_normal(gen);
      }
      dx *= 0.5*prm.get<double>("ds_init")/dx.norm();

      Vector3d x_a = x0_ + dx;
      bool inside_a = intp->locate(x_a);
      Vector3d x_b = x0_ - dx;
      bool inside_b = intp->locate(x_b);
      if (inside_a && inside_b){
        //std::cout << "INSIDE" << std::endl;
        // std::cout << "INSIDE: " << x_a << " " << x_b << std::endl; 
        nodes.push_back(x_a);
        nodes.push_back(x_b);
        double ds0 = (x_a-x_b).norm();
        edges.push_back({{2*ipair, 2*ipair+1}, ds0});
        ++ipair;
      }
      else if (key[0] == "pair"){
        std::cout << "Pair not inside domain" << std::endl;
        exit(1);
      }
      //std::cout << x_a << " " << x_b << std::endl;
    }
  };
};

class RandomPointsInitializer : public Initializer {
protected:
  std::mt19937 &gen;
public:
  RandomPointsInitializer( const std::vector<std::string>& key
                         , std::shared_ptr<Interpol> intp
                         , partrac::Params& prm
                         , std::mt19937 &gen
                         ) : Initializer(intp, prm), gen(gen) {
    bool init_rand_x = false;
    bool init_rand_y = false;
    bool init_rand_z = false;
    if (contains(key[1], "x"))
      init_rand_x = true;
    if (contains(key[1], "y"))
      init_rand_y = true;
    if (contains(key[1], "z"))
      init_rand_z = true;

    // TODO: Factor out position generation
    double tol = 1e-12;
    Uint N_est = 1000000;
    double dx_est;
  
    double Lx = L[0];
    double Ly = L[1];
    double Lz = L[2];

    std::cout << "FML: " << Lx << " " << Ly << " " << Lz << std::endl;

    bool hasLx = Lx > tol and init_rand_x;
    bool hasLy = Ly > tol and init_rand_y;
    bool hasLz = Lz > tol and init_rand_z;

    if (hasLx && hasLy && hasLz){
      dx_est = pow(Lx*Ly*Lz/N_est, 1./3);
    }
    else if (hasLx && hasLy){
      dx_est = pow(Lx*Ly/N_est, 1./2);
    }
    else if (hasLx && hasLz){
      dx_est = pow(Lx*Lz/N_est, 1./2);
    }
    else if (hasLy && hasLz){
      dx_est = pow(Ly*Lz/N_est, 1./2);
    }
    else if (hasLx){
      dx_est = Lx/N_est;
    }
    else if (hasLy){
      dx_est = Ly/N_est;
    }
    else if (hasLz){
      dx_est = Lz/N_est;
    }
    else {
      std::cout << "Something is wrong with the domain!" << std::endl;
      exit(1);
    }
    Uint Nx = 1;
    Uint Ny = 1;
    Uint Nz = 1;
    if (hasLx) Nx = Lx/dx_est+1;
    if (hasLy) Ny = Ly/dx_est+1;
    if (hasLz) Nz = Lz/dx_est+1;
    double dx = Lx/Nx;
    double dy = Ly/Ny;
    double dz = Lz/Nz;

    double ww;

    std::cout << "dx: " << dx << " " << dy << " " << dz << std::endl;
    std::cout << "Nx: " << Nx << " " << Ny << " " << Nz << std::endl;

    std::vector<double> wei;
    std::vector<Vector3d> pos;
    for (Uint ix=0; ix<Nx; ++ix){
      for (Uint iy=0; iy<Ny; ++iy){
        for (Uint iz=0; iz<Nz; ++iz){
          Vector3d x = x0;
          if (hasLx) x[0] = x_min[0]+(ix+0.5)*dx;
          if (hasLy) x[1] = x_min[1]+(iy+0.5)*dy;
          if (hasLz) x[2] = x_min[2]+(iz+0.5)*dz;
          PointValues ptvals(intp->get_U0());
          intp->evaluate(x, ptvals);
          if (prm.get<std::string>("init_weight") == "ux"){
            ww = abs(ptvals.U[0]);
          }
          else if (prm.get<std::string>("init_weight") == "uy"){
            ww = abs(ptvals.U[1]);
          }
          else if (prm.get<std::string>("init_weight") == "uz"){
            ww = abs(ptvals.U[2]);
          }
          else if (prm.get<std::string>("init_weight") == "u"){
            ww = ptvals.U.norm();
          }
          else {
            ww = 1.;
          }

          wei.push_back(ww);
          pos.push_back(x);
        }
      }
    }
    std::uniform_real_distribution<> uni_dist_dx(-0.5*dx, 0.5*dx);
    std::uniform_real_distribution<> uni_dist_dy(-0.5*dy, 0.5*dy);
    std::uniform_real_distribution<> uni_dist_dz(-0.5*dz, 0.5*dz);
    std::discrete_distribution<Uint> discrete_dist(wei.begin(), wei.end());

    for (Uint irw=0; irw<prm.get<Uint>("Nrw"); ++irw){
      Vector3d x;
      bool inside = false;
      do {
        Uint ind = discrete_dist(gen);
        x = pos[ind];
        if (hasLx) x[0] += uni_dist_dx(gen);
        if (hasLy) x[1] += uni_dist_dy(gen);
        if (hasLz) x[2] += uni_dist_dz(gen);
        inside = intp->locate(x);
      } while (!inside);
      nodes.push_back(x);
    }

    sort(nodes.begin(), nodes.end(), less_than_op());

    for (Uint irw=1; irw < prm.get<Uint>("Nrw"); ++irw){
      double ds0 = dist(nodes[irw-1], nodes[irw]);
      if (ds0 < 10*prm.get<double>("ds_init"))  // 2 lattice units (before) --> 10 x ds_max (now)
        edges.push_back({{irw-1, irw}, ds0});
      // Needs customization for 2D/3D applications
    }
  };
};

class RandomGaussianStripInitializer : public Initializer {
protected:
  std::mt19937 &gen;
public:
  RandomGaussianStripInitializer( const std::vector<std::string>& key
                                , std::shared_ptr<Interpol> intp
                                , partrac::Params& prm
                                , std::mt19937 &gen
                                ) : Initializer(intp, prm), gen(gen) {
    edges.clear();
    faces.clear();

    double La = prm.get<double>("La");
    double sigma0 = prm.get<double>("Lb");
    
    Vector3d n(0., 0., 0.);
    if (key[1] == "x"){
      n[0] = 1.0;
    }
    if (key[1] == "y"){
      n[1] = 1.0;
    }
    if (key[1] == "z"){
      n[2] = 1.0;
    }

    bool init_rand_x = contains(key[2], "x");
    bool init_rand_y = contains(key[2], "y");
    bool init_rand_z = contains(key[2], "z");
    
    std::uniform_real_distribution<double> rnd_unit(0.0, 1.0);
    std::normal_distribution<double> rnd_normal(0.0, 1.0);

    Vector3d x00 = x0;
    Vector3d x01 = x0;
    for (Uint dim=0; dim < 3; ++dim){
      x00[dim] += -La/2*n[dim];
      x01[dim] += La/2*n[dim];
    }

    Uint failed_attempts = 0;
    Uint max_failed_attempts = 1000000; // Maybe not hardcode?

    Uint irw = 0;
    while (irw < prm.get<Uint>("Nrw") && failed_attempts < max_failed_attempts){
      double alpha = rnd_unit(gen);
      Vector3d xi = alpha * x00 + (1.-alpha) * x01;
      if (init_rand_x)
        xi[0] += sigma0 * rnd_normal(gen);
      if (init_rand_y)
        xi[1] += sigma0 * rnd_normal(gen);
      if (init_rand_z)
        xi[2] += sigma0 * rnd_normal(gen);
      // check if inside domain
      if (intp->locate(xi)){
        nodes.push_back(xi);
        ++irw;
        failed_attempts = 0;
      }
      else {
        ++failed_attempts;
      }
    }
    if (irw == 0) {
      std::cout << "No points inside domain" << std::endl;
      exit(1);
    }
  };
};

class RandomGaussianCircleInitializer : public Initializer {
protected:
  std::mt19937 &gen;
public:
  RandomGaussianCircleInitializer( const std::vector<std::string>& key
                                , std::shared_ptr<Interpol> intp
                                , partrac::Params& prm
                                , std::mt19937 &gen
                                ) : Initializer(intp, prm), gen(gen) {
    edges.clear();
    faces.clear();

    double R = prm.get<double>("La")/2;
    double sigma0 = prm.get<double>("Lb");
    
    Vector3d t1(0., 0., 0.);
    Vector3d t2(0., 0., 0.);
    if (key[1] == "x"){
      t1[1] = 1.0;
      t2[2] = 1.0;
    }
    if (key[1] == "y"){
      t1[0] = 1.0;
      t2[2] = 1.0;
    }
    if (key[1] == "z"){
      t1[0] = 1.0;
      t2[1] = 1.0;
    }

    std::uniform_real_distribution<double> rnd_unit(0.0, 1.0);
    std::normal_distribution<double> rnd_normal(0.0, 1.0);

    Uint failed_attempts = 0;
    Uint max_failed_attempts = 1000000; // Maybe not hardcode?

    Uint irw = 0;
    while (irw < prm.get<Uint>("Nrw") && failed_attempts < max_failed_attempts){
      double alpha1 = 1.;
      double alpha2 = 1.;
      while (pow(alpha1, 2) + pow(alpha2, 2) > 1){
        alpha1 = 2*rnd_unit(gen)-1;
        alpha2 = 2*rnd_unit(gen)-1;
      }
      
      Vector3d xi = x0 + R * (alpha1 * t1 + alpha2 * t2);
      for (Uint dim=0; dim<3; ++dim)
        xi[dim] += sigma0 * rnd_normal(gen);

      // check if inside domain
      if (intp->locate(xi)){
        nodes.push_back(xi);
        ++irw;
        failed_attempts = 0;
      }
      else {
        ++failed_attempts;
      }
    }
    if (irw == 0) {
      std::cout << "No points inside domain" << std::endl;
      exit(1);
    }
  };
};

class FileInitializer : public Initializer {
public:
  FileInitializer(const std::vector<std::string>& key_col, std::shared_ptr<Interpol> intp, partrac::Params& prm) : Initializer(intp, prm) {

    std::string h5filename = key_col[1];

    verify_file_exists(h5filename);

    H5::H5File h5file(h5filename, H5F_ACC_RDONLY);
    H5::DataSet dset_nodes = h5file.openDataSet("nodes");
    H5::DataSpace dspace_nodes = dset_nodes.getSpace();
    hsize_t dims_nodes[2];
    dspace_nodes.getSimpleExtentDims(dims_nodes, NULL);
    
    std::vector<double> nodes_buf(dims_nodes[0]*dims_nodes[1]);
    dset_nodes.read(nodes_buf.data(), H5::PredType::NATIVE_DOUBLE, dspace_nodes, dspace_nodes);

   
    h5file.close();

    bool this_inside = false;
    bool prev_inside = false; 
    Uint irw = 0;

    int cell_id = -1;
    for (Uint i=0; i < dims_nodes[0]; ++i){
      Vector3d xi = {prm.get<double>("x0"), prm.get<double>("y0"), prm.get<double>("z0")};
      for (Uint j=0; j < dims_nodes[1]; ++j){
        xi[j] = nodes_buf[i * dims_nodes[1] + j];
      }
      // check if inside domain
      this_inside = intp->locate(xi, prm.get<double>("t0"), cell_id);
      if (this_inside){
        nodes.push_back(xi);
        if (prev_inside)
          edges.push_back({{irw-1, irw}, dist(nodes[irw-1], xi)});
        ++irw;
      }
      prev_inside = this_inside;
    }
    if (irw == 0) {
      std::cout << "No points inside domain" << std::endl;
      exit(1);
    }
  };
};


#endif
