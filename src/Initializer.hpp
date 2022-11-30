#ifndef __DISTRIBUTE_HPP
#define __DISTRIBUTE_HPP

#include <iomanip>
#include <vector>
#include <random>
#include <algorithm>
#include "typedefs.hpp"
#include "Interpol.hpp"
#include "utils.hpp"
#include "MPIwrap.hpp"

// TODO: Massive cleanup!
struct less_than_op {
  inline bool operator() (const Vector3d &a, const Vector3d &b){
    return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]) || (a[0] == b[0] && a[1] == b[1] && a[2] < b[2]);
  }
};

class Initializer {
public:
  Initializer(std::shared_ptr<Interpol> intp, Parameters& prm, MPIwrap& mpi) : intp(intp), prm(prm), m_mpi(mpi) {
    x0 = {prm.x0, prm.y0, prm.z0};
    x_min = intp->get_x_min();
    x_max = intp->get_x_max();
    L = x_max - x_min;
    inject = prm.inject;
    clear_initial_edges = prm.clear_initial_edges;
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
  Parameters& prm;
  Vector3d x0;
  Vector3d x_min;
  Vector3d x_max;
  Vector3d L;
  MPIwrap& m_mpi;
};

class PointInitializer : public Initializer {
public:
  PointInitializer(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, Parameters& prm, MPIwrap& mpi) : Initializer(intp, prm, mpi) {
    intp->probe(x0);
    for (Uint irw=0; irw < prm.Nrw; ++irw){
      if (intp->inside_domain()){
        nodes.push_back(x0);
      }
    }
    edges.clear();
    faces.clear();
  };
};

class UniformInitializer : public Initializer {
public:
  UniformInitializer(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, Parameters& prm, MPIwrap& mpi) : Initializer(intp, prm, mpi) {
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
      exit(0);
    }
    Vector3d Dx = (x_b - x_a) / (prm.Nrw-1);
    for (Uint irw=0; irw < prm.Nrw; ++irw){
      Vector3d x = x_a + Dx * irw;
      intp->probe(x);
      if (intp->inside_domain()){
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
  StripInitializer(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, Parameters& prm, MPIwrap& mpi) : Initializer(intp, prm, mpi) {
    double La = prm.La;
    double Lb = prm.Lb;
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

    intp->probe(x00);
    bool inside_00 = intp->inside_domain();
    intp->probe(x01);
    bool inside_01 = intp->inside_domain();

    if (inside_00 && inside_01){
      std::cout << "Strip inside domain." << std::endl;
    }
    else {
      std::cout << "Strip not inside domain" << std::endl;
      exit(0);
    }

    nodes.push_back(x0);
    for (Uint irw=0; irw < prm.Nrw; ++irw){
      double alpha = float(irw+1)/prm.Nrw;
      Vector3d x0i = alpha * x00 + (1.-alpha) * x01;
      // check if inside domain
      nodes.push_back(x0i);
      edges.push_back({{irw, irw+1}, dist(nodes[irw], nodes[irw+1])});
    }
  };
};

class SheetInitializer : public Initializer {
public:
  SheetInitializer(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, Parameters& prm, MPIwrap& mpi) : Initializer(intp, prm, mpi) {
    double La = prm.La;
    double Lb = prm.Lb;
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

    intp->probe(x00);
    bool inside_00 = intp->inside_domain();
    intp->probe(x01);
    bool inside_01 = intp->inside_domain();
    intp->probe(x10);
    bool inside_10 = intp->inside_domain();
    intp->probe(x11);
    bool inside_11 = intp->inside_domain();
    if (inside_00 && inside_01 && inside_10 && inside_11){
      std::cout << "Sheet inside domain." << std::endl;
    }
    else {
      std::cout << "Sheet not inside domain" << std::endl;
      exit(0);
    }

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
  };
};

class EllipsoidInitializer : public Initializer {
public:
  EllipsoidInitializer(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, Parameters& prm, MPIwrap& mpi) : Initializer(intp, prm, mpi) {
    double La = prm.La;
    double Lb = prm.Lb;
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

    Vector3d x_0 = {-R/sqrt(2.),  -R/sqrt(6.0), -R/sqrt(3.0)/2};
    Vector3d x_1 = { R/sqrt(2.),  -R/sqrt(6.0), -R/sqrt(3.0)/2};
    Vector3d x_2 = {         0., R*sqrt(2./3.), -R/sqrt(3.0)/2};
    Vector3d x_3 = {         0.,            0.,  R*sqrt(3.0)/2};

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


    intp->probe(x_0);
    bool inside_0 = intp->inside_domain();
    intp->probe(x_1);
    bool inside_1 = intp->inside_domain();
    intp->probe(x_2);
    bool inside_2 = intp->inside_domain();
    intp->probe(x_3);
    bool inside_3 = intp->inside_domain();

    if (inside_0 && inside_1 && inside_2 && inside_3){
      std::cout << "Ellipsoid inside domain." << std::endl;
    }
    else {
      std::cout << "Ellipsoid not inside domain" << std::endl;
      exit(0);
    }

    std::vector<Vector3d> nodes_loc;
    nodes_loc.push_back(x_0);
    nodes_loc.push_back(x_1);
    nodes_loc.push_back(x_2);
    nodes_loc.push_back(x_3);

    ParticleSet pset_loc(intp, prm.Nrw_max, m_mpi);
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
                                    pset_loc, prm.ds_max, 0.0, false);
      for (Uint irw=0; irw<pset_loc.N(); ++irw){
        Vector3d x = pset_loc.x(irw);
        Vector3d nn = x / x.norm();
        double rad = 1./sqrt(nn[0]*nn[0]/lx2 + nn[1]*nn[1]/ly2 + nn[2]*nn[2]/lz2);
        pset_loc.set_x(irw, rad * nn);
      }
      n_rem = sheet_coarsening(faces, edges, edge2faces_loc, node2edges_loc, edges_inlet_dummy, nodes_inlet_dummy,
                               pset_loc, prm.ds_min, 0.0);

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
  RandomPairsInitializer(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, Parameters& prm, MPIwrap& mpi, std::mt19937 &gen) : Initializer(intp, prm, mpi), gen(gen) {
    std::uniform_real_distribution<> uni_dist_x(x_min[0], x_max[0]);
    std::uniform_real_distribution<> uni_dist_y(x_min[1], x_max[1]);
    std::uniform_real_distribution<> uni_dist_z(x_min[2], x_max[2]);
    std::normal_distribution<double> rnd_normal(0.0, 1.0);

    Uint Npairs = (key[0] == "pair") ? 1 : prm.Nrw/2;

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
      dx *= 0.5*prm.ds_init/dx.norm();

      Vector3d x_a = x0_ + dx;
      intp->probe(x_a);
      bool inside_a = intp->inside_domain();
      Vector3d x_b = x0_ - dx;
      intp->probe(x_b);
      bool inside_b = intp->inside_domain();
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
        exit(0);
      }
      //std::cout << x_a << " " << x_b << std::endl;
    }
  };
};

class RandomPointsInitializer : public Initializer {
protected:
  std::mt19937 &gen;
public:
  RandomPointsInitializer(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, Parameters& prm, MPIwrap& mpi, std::mt19937 &gen) : Initializer(intp, prm, mpi), gen(gen) {
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

    if (Lx > tol && Ly > tol && Lz > tol){
      dx_est = pow(Lx*Ly*Lz/N_est, 1./3);
    }
    else if (Lx > tol && Ly > tol){
      dx_est = pow(Lx*Ly/N_est, 1./2);
    }
    else if (Lx > tol && Lz > tol){
      dx_est = pow(Lx*Lz/N_est, 1./2);
    }
    else if (Ly > tol && Lz > tol){
      dx_est = pow(Ly*Lz/N_est, 1./2);
    }
    else if (Lx > tol){
      dx_est = Lx/N_est;
    }
    else if (Ly > tol){
      dx_est = Ly/N_est;
    }
    else if (Lz > tol){
      dx_est = Lz/N_est;
    }
    else {
      std::cout << "Something is wrong with the domain!" << std::endl;
      exit(0);
    }
    Uint Nx = 0;
    Uint Ny = 0;
    Uint Nz = 0;
    if (init_rand_x) Nx = Lx/dx_est+1;
    if (init_rand_y) Ny = Ly/dx_est+1;
    if (init_rand_z) Nz = Lz/dx_est+1;
    double dx = Lx/Nx;
    double dy = Ly/Ny;
    double dz = Lz/Nz;

    double ww;

    std::vector<double> wei;
    std::vector<Vector3d> pos;
    for (Uint ix=0; ix<Nx; ++ix){
      for (Uint iy=0; iy<Ny; ++iy){
        for (Uint iz=0; iz<Nz; ++iz){
          Vector3d x = x0;
          if (init_rand_x) x[0] = x_min[0]+(ix+0.5)*dx;
          if (init_rand_y) x[1] = x_min[1]+(iy+0.5)*dy;
          if (init_rand_z) x[2] = x_min[2]+(iz+0.5)*dz;
          intp->probe(x);
          if (prm.init_weight == "ux"){
            ww = abs(intp->get_ux());
          }
          else if (prm.init_weight == "uy"){
            ww = abs(intp->get_uy());
          }
          else if (prm.init_weight == "uz"){
            ww = abs(intp->get_uz());
          }
          else if (prm.init_weight == "u"){
            ww = sqrt(pow(intp->get_ux(), 2)
                      + pow(intp->get_uy(), 2)
                      + pow(intp->get_uz(), 2));
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

    for (Uint irw=0; irw<prm.Nrw; ++irw){
      Vector3d x;
      do {
        Uint ind = discrete_dist(gen);
        x = pos[ind];
        if (init_rand_x) x[0] += uni_dist_dx(gen);
        if (init_rand_y) x[1] += uni_dist_dy(gen);
        if (init_rand_z) x[2] += uni_dist_dz(gen);
        intp->probe(x);
        //std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
      } while (!intp->inside_domain());
      nodes.push_back(x);
    }

    sort(nodes.begin(), nodes.end(), less_than_op());

    for (Uint irw=1; irw < prm.Nrw; ++irw){
      double ds0 = dist(nodes[irw-1], nodes[irw]);
      if (ds0 < 10*prm.ds_init)  // 2 lattice units (before) --> 10 x ds_max (now)
        edges.push_back({{irw-1, irw}, ds0});
      // Needs customization for 2D/3D applications
    }
  };
};

/*std::vector<Vector3d> initial_positions(const std::string init_mode,
                                        const std::string init_weight,
                                        Uint &Nrw,
                                        const Vector3d &x0,
                                        const double La,
                                        const double Lb,
                                        const double ds,
                                        const double t0,
                                        Interpol *intp,
                                        std::mt19937 &gen,
                                        EdgesType &edges,
                                        FacesType &faces
                                   ){
  intp->update(t0);

  Vector3d x;
  double Lx = intp->get_Lx();
  double Ly = intp->get_Ly();
  double Lz = intp->get_Lz();
  Vector3d x_min = intp->get_x_min();
  Vector3d x_max = intp->get_x_max();

  Uint Nx = 1;
  Uint Ny = 1;
  Uint Nz = 1;

  std::vector<std::string> key = split_string(init_mode, "_");

  if (key.size() == 0){
    std::cout << "init_mode not specified." << std::endl;
    exit(0);
  }

  if (key[0] == "point"){
    std::vector<Vector3d> pos_init;
    intp->probe(x0);
    for (Uint irw=0; irw < Nrw; ++irw){
      if (intp->inside_domain()){
        pos_init.push_back(x0);
      }
    }
    edges.clear();
    return pos_init;
  }

  if (key[0] == "uniform"){
    std::vector<Vector3d> pos_init;

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
      exit(0);
    }
    Vector3d Dx = (x_b - x_a) / (Nrw-1);
    for (Uint irw=0; irw < Nrw; ++irw){
      x = x_a + Dx * irw;
      intp->probe(x);
      if (intp->inside_domain()){
        pos_init.push_back(x);
      }
    }
    for (Uint irw=0; irw < pos_init.size()-1; ++irw){
      if ((pos_init[irw] - pos_init[irw+1]).norm() < 1.5*Dx.norm()){
        edges.push_back({{irw, irw+1}, dist(pos_init[irw], pos_init[irw+1])});
      }
    }
    return pos_init;
  }

  if (key[0] == "sheet"){
    std::vector<Vector3d> pos_init;

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

    intp->probe(x00);
    bool inside_00 = intp->inside_domain();
    intp->probe(x01);
    bool inside_01 = intp->inside_domain();
    intp->probe(x10);
    bool inside_10 = intp->inside_domain();
    intp->probe(x11);
    bool inside_11 = intp->inside_domain();
    if (inside_00 && inside_01 && inside_10 && inside_11){
      std::cout << "Sheet inside domain." << std::endl;
    }
    else {
      std::cout << "Sheet not inside domain" << std::endl;
      exit(0);
    }

    pos_init.push_back(x00);
    pos_init.push_back(x01);
    pos_init.push_back(x10);
    pos_init.push_back(x11);

    edges.push_back({{0, 1}, dist(pos_init[0], pos_init[1])});
    edges.push_back({{0, 2}, dist(pos_init[0], pos_init[2])});
    edges.push_back({{1, 2}, dist(pos_init[1], pos_init[2])});
    edges.push_back({{2, 3}, dist(pos_init[2], pos_init[3])});
    edges.push_back({{3, 0}, dist(pos_init[3], pos_init[4])});

    faces.push_back({{0, 2, 1}, La*Lb/2});
    faces.push_back({{1, 3, 4}, La*Lb/2});

    Nrw = pos_init.size();
    return pos_init;
  }
  else if (key[0] == "pair" || key[0] == "pairs"){
    std::vector<Vector3d> pos_init;

    std::uniform_real_distribution<> uni_dist_x(x_min[0], x_max[0]);
    std::uniform_real_distribution<> uni_dist_y(x_min[1], x_max[1]);
    std::uniform_real_distribution<> uni_dist_z(x_min[2], x_max[2]);
    std::normal_distribution<double> rnd_normal(0.0, 1.0);

    //std::cout << "x_min = " << x_min << std::endl;
    //std::cout << "x_max = " << x_max << std::endl;

    // std::uniform_real_distribution<> uni_dist_theta(0., 2*M_PI);
    Uint Npairs = (key[0] == "pair") ? 1 : Nrw/2;

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
      dx *= 0.5*ds/dx.norm();

      Vector3d x_a = x0_ + dx;
      intp->probe(x_a);
      bool inside_a = intp->inside_domain();
      Vector3d x_b = x0_ - dx;
      intp->probe(x_b);
      bool inside_b = intp->inside_domain();
      if (inside_a && inside_b){
        //std::cout << "INSIDE" << std::endl;
        // std::cout << "INSIDE: " << x_a << " " << x_b << std::endl; 
        pos_init.push_back(x_a);
        pos_init.push_back(x_b);
        double ds0 = (x_a-x_b).norm();
        edges.push_back({{2*ipair, 2*ipair+1}, ds0});
        ++ipair;
      }
      else if (key[0] == "pair"){
        std::cout << "Pair not inside domain" << std::endl;
        exit(0);
      }
      //std::cout << x_a << " " << x_b << std::endl;
    }
    Nrw = 2*Npairs;
    return pos_init;
  }

  bool init_rand_x = false;
  bool init_rand_y = false;
  bool init_rand_z = false;
  if (init_mode == "line_x" ||
      init_mode == "plane_xy" ||
      init_mode == "plane_xz" ||
      init_mode == "volume"){
    init_rand_x = true;
  }
  if (init_mode == "line_y" ||
      init_mode == "plane_xy" ||
      init_mode == "plane_yz" ||
      init_mode == "volume"){
    init_rand_y = true;

  }
  if (init_mode == "line_z" ||
      init_mode == "plane_xz" ||
      init_mode == "plane_yz" ||
      init_mode == "volume"){
    init_rand_z = true;
  }

  // TODO: Factor out position generation
  double tol = 1e-12;
  Uint N_est = 1000000;
  double dx_est;
  if (Lx > tol && Ly > tol && Lz > tol){
    dx_est = pow(Lx*Ly*Lz/N_est, 1./3);
  }
  else if (Lx > tol && Ly > tol){
    dx_est = pow(Lx*Ly/N_est, 1./2);
  }
  else if (Lx > tol && Lz > tol){
    dx_est = pow(Lx*Lz/N_est, 1./2);
  }
  else if (Ly > tol && Lz > tol){
    dx_est = pow(Ly*Lz/N_est, 1./2);
  }
  else if (Lx > tol){
    dx_est = Lx/N_est;
  }
  else if (Ly > tol){
    dx_est = Ly/N_est;
  }
  else if (Lz > tol){
    dx_est = Lz/N_est;
  }
  else {
    std::cout << "Something is wrong with the domain!" << std::endl;
    exit(0);
  }
  if (init_rand_x) Nx = Lx/dx_est+1;
  if (init_rand_y) Ny = Ly/dx_est+1;
  if (init_rand_z) Nz = Lz/dx_est+1;
  double dx = Lx/Nx;
  double dy = Ly/Ny;
  double dz = Lz/Nz;

  double ww;

  std::vector<double> wei;
  std::vector<Vector3d> pos;
  for (Uint ix=0; ix<Nx; ++ix){
    for (Uint iy=0; iy<Ny; ++iy){
      for (Uint iz=0; iz<Nz; ++iz){
        x = x0;
        if (init_rand_x) x[0] = x_min[0]+(ix+0.5)*dx;
        if (init_rand_y) x[1] = x_min[1]+(iy+0.5)*dy;
        if (init_rand_z) x[2] = x_min[2]+(iz+0.5)*dz;
        intp->probe(x);
        if (init_weight == "ux"){
          ww = abs(intp->get_ux());
        }
        else if (init_weight == "uy"){
          ww = abs(intp->get_uy());
        }
        else if (init_weight == "uz"){
          ww = abs(intp->get_uz());
        }
        else if (init_weight == "u"){
          ww = sqrt(pow(intp->get_ux(), 2)
                    + pow(intp->get_uy(), 2)
                    + pow(intp->get_uz(), 2));
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

  std::vector<Vector3d> pos_init;
  for (Uint irw=0; irw<Nrw; ++irw){
    do {
      Uint ind = discrete_dist(gen);
      x = pos[ind];
      if (init_rand_x) x[0] += uni_dist_dx(gen);
      if (init_rand_y) x[1] += uni_dist_dy(gen);
      if (init_rand_z) x[2] += uni_dist_dz(gen);
      intp->probe(x);
      //std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
    } while (!intp->inside_domain());
    pos_init.push_back(x);
  }

  sort(pos_init.begin(), pos_init.end(), less_than_op());

  for (Uint irw=1; irw < Nrw; ++irw){
    double ds0 = dist(pos_init[irw-1], pos_init[irw]);
    if (ds0 < 10*ds)  // 2 lattice units (before) --> 10 x ds_max (now)
      edges.push_back({{irw-1, irw}, ds0});
    // Needs customization for 2D/3D applications
  }

  return pos_init;
}*/

#endif
