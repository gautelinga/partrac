#ifndef __EXP_INITIALIZER_HPP
#define __EXP_INITIALIZER_HPP

#include <iomanip>
#include <vector>
#include <random>
#include <algorithm>

#include "typedefs.hpp"
#include "../Interpol.hpp"
#include "utils.hpp"
#include "../MPIwrap.hpp"
//#include "particles.hpp"

// TODO: Massive cleanup!
struct less_than_op {
  inline bool operator() (const Vector &a, const Vector &b){
    return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]) || (a[0] == b[0] && a[1] == b[1] && a[2] < b[2]);
  }
};

class Initializer {
public:
  //Initializer(IntpType& intp, Parameters& prm, MPIwrap& mpi) : m_intp(intp), prm(prm), m_mpi(mpi) {
  Initializer(Parameters& prm, MPIwrap& mpi) : prm(prm), m_mpi(mpi) {
    x0 = {prm.x0, prm.y0, prm.z0};
    //x_min = intp.get_x_min();
    //x_max = intp.get_x_max();
    //L = x_max - x_min;
    inject = prm.inject;
    clear_initial_edges = prm.clear_initial_edges;
  };
  virtual ~Initializer() { nodes.clear(); edges.clear(); faces.clear(); };
  //virtual void probe(IntpType& intp) = 0;
  /*std::vector<Vector>::const_iterator node_begin() const { return nodes.begin(); };
  std::vector<Vector>::const_iterator node_end() const { return nodes.end(); };
  EdgesType::const_iterator edge_begin() const { return edges.begin(); };
  EdgesType::const_iterator edge_end() const { return edges.end(); };
  FacesType::const_iterator face_begin() const { return faces.begin(); };
  FacesType::const_iterator face_end() const { return faces.end(); };*/
  std::vector<Vector> nodes;
  EdgesType edges;
  FacesType faces;
  bool inject;
  bool clear_initial_edges;
  // IntpType& interpolator() { return m_intp; };
  template<typename T> void initialize(T& particles);
protected:
  //IntpType& m_intp;
  Parameters& prm;
  Vector x0;
  Vector x_min;
  Vector x_max;
  Vector L;
  MPIwrap& m_mpi;
};

template<typename T>
void Initializer::initialize(T& ps)
{
  std::cout << "Nodes: " << nodes.size() << std::endl;
  for (auto & x : nodes)
  {
    //std::cout << "adding particle" << std::endl;
    ps.add_particle(x);
  }
  ps.color_particles(0., 1.);
  for (auto & inedge : edges)
  {
    Uint i = inedge.first[0];
    Uint j = inedge.first[1];
    Real w = inedge.second;
    //ps.add_edge(ps.particle_ptr(i), ps.particle_ptr(j), w);
    //std::cout << "adding edge" << std::endl;
    ps.add_edge(i, j, w);
  }
  for (auto & inface : faces)
  {
    Uint i = inface.first[0];
    Uint j = inface.first[1];
    Uint k = inface.first[2];
    Real w = inface.second;
    //ps.add_face(ps.edge_ptr(i), ps.edge_ptr(j), ps.edge_ptr(k), w);
    ps.add_face(i, j, k, w);
  }
}

/*
template<class InterpolType>
class PointInitializer : public Initializer<InterpolType> {
public:
  PointInitializer<InterpolType>(const std::vector<std::string>& key, InterpolType& intp, Parameters& prm, MPIwrap& mpi)
   : Initializer<InterpolType>(intp, prm, mpi) {
    this->interpolator().probe(this->x0);
    for (Uint irw=0; irw < prm.Nrw; ++irw){
      if (this->interpolator().inside_domain()){
        this->nodes.push_back(this->x0);
      }
    }
    this->edges.clear();
    this->faces.clear();
  };
};
*/

class UniformInitializer : public Initializer {
protected:  
  std::vector<std::string> key;
public:
  UniformInitializer( const std::vector<std::string>& key
                    //, IntpType& intp
                    , Parameters& prm
                    , MPIwrap& mpi
                    //) : Initializer(intp, prm, mpi) {
                    ) : Initializer(prm, mpi), key(key) {
    //probe(intp);
  };
  template<typename IntpType> 
  void probe(IntpType& intp);
};

template<typename IntpType>
void UniformInitializer::probe(IntpType& intp){
  x_min = intp.get_x_min();
  x_max = intp.get_x_max();
  L = x_max - x_min;

  Vector x_a = this->x0;
  Vector x_b = this->x0;
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
  Vector Dx = (x_b - x_a) / (prm.Nrw-1);
  for (Uint irw=0; irw < prm.Nrw; ++irw){
    Vector x = x_a + Dx * irw;
    intp.probe(x);
    if (intp.inside_domain()){
      this->nodes.push_back(x);
    }
  }
  for (Uint irw=0; irw < this->nodes.size()-1; ++irw){
    if ((this->nodes[irw] - this->nodes[irw+1]).norm() < 1.5*Dx.norm()){
      this->edges.push_back({{irw, irw+1}, dist(this->nodes[irw], this->nodes[irw+1])});
    }
  }
};

/*
template <class InterpolType>
class SheetInitializer : Initializer<InterpolType> {
public:
  SheetInitializer(const std::vector<std::string>& key, InterpolType& intp, Parameters& prm, MPIwrap& mpi) : Initializer<InterpolType>(intp, prm, mpi) {
    Real La = prm.La;
    Real Lb = prm.Lb;
    Vector n(0., 0., 0.);
    Vector ta(0., 0., 0.);
    Vector tb(0., 0., 0.);
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

    Vector x00 = this->x0;
    x00[0] += - La/2*ta[0] - Lb/2*tb[0];
    x00[1] += - La/2*ta[1] - Lb/2*tb[1];
    x00[2] += - La/2*ta[2] - Lb/2*tb[2];

    Vector x01 = this->x0;
    x01[0] += La/2*ta[0] + Lb/2*tb[0];
    x01[1] += La/2*ta[1] - Lb/2*tb[1];
    x01[2] += - La/2*ta[2] - Lb/2*tb[2];

    Vector x10 = this->x0;
    x10[0] += La/2*ta[0] + Lb/2*tb[0];
    x10[1] += La/2*ta[1] + Lb/2*tb[1];
    x10[2] += La/2*ta[2] + Lb/2*tb[2];

    Vector x11 = this->x0;
    x11[0] += - La/2*ta[0] - Lb/2*tb[0];
    x11[1] += - La/2*ta[1] + Lb/2*tb[1];
    x11[2] += - La/2*ta[2] + Lb/2*tb[2];

    this->interpolator().probe(x00);
    bool inside_00 = this->interpolator().inside_domain();
    this->interpolator().probe(x01);
    bool inside_01 = this->interpolator().inside_domain();
    this->interpolator().probe(x10);
    bool inside_10 = this->interpolator().inside_domain();
    this->interpolator().probe(x11);
    bool inside_11 = this->interpolator().inside_domain();
    if (inside_00 && inside_01 && inside_10 && inside_11){
      std::cout << "Sheet inside domain." << std::endl;
    }
    else {
      std::cout << "Sheet not inside domain" << std::endl;
      exit(0);
    }

    this->nodes.push_back(x00);#include <random>

    this->nodes.push_back(x11);

    this->edges.push_back({{0, 1}, dist(this->nodes[0], this->nodes[1])});
    this->edges.push_back({{0, 2}, dist(this->nodes[0], this->nodes[2])});
    this->edges.push_back({{1, 2}, dist(this->nodes[1], this->nodes[2])});
    this->edges.push_back({{2, 3}, dist(this->nodes[2], this->nodes[3])});
    this->edges.push_back({{3, 0}, dist(this->nodes[3], this->nodes[4])});

    this->faces.push_back({{0, 2, 1}, La*Lb/2});
    this->faces.push_back({{1, 3, 4}, La*Lb/2});
  };
};

template <class InterpolType>
class EllipsoidInitializer : public Initializer<InterpolType> {
public:
  EllipsoidInitializer<InterpolType>(const std::vector<std::string>& key, InterpolType& intp, Parameters& prm, MPIwrap& mpi) : Initializer<InterpolType>(intp, prm, mpi) {
    Real La = prm.La;
    Real Lb = prm.Lb;
    Real lx2 = Lb*Lb;
    Real ly2 = Lb*Lb;
    Real lz2 = Lb*Lb;
    Vector n(0., 0., 0.);
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

    Real R = sqrt(La*Lb);

    Vector x_0 = {-R/sqrt(2.),  -R/sqrt(6.0), -R/sqrt(3.0)/2};
    Vector x_1 = { R/sqrt(2.),  -R/sqrt(6.0), -R/sqrt(3.0)/2};
    Vector x_2 = {         0., R*sqrt(2./3.), -R/sqrt(3.0)/2};
    Vector x_3 = {         0.,            0.,  R*sqrt(3.0)/2};

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


    this->interpolator().probe(x_0);
    bool inside_0 = this->interpolator().inside_domain();
    this->interpolator().probe(x_1);
    bool inside_1 = this->interpolator().inside_domain();
    this->interpolator().probe(x_2);
    bool inside_2 = this->interpolator().inside_domain();
    this->interpolator().probe(x_3);
    bool inside_3 = this->interpolator().inside_domain();

    if (inside_0 && inside_1 && inside_2 && inside_3){
      std::cout << "Ellipsoid inside domain." << std::endl;
    }
    else {
      std::cout << "Ellipsoid not inside domain" << std::endl;
      exit(0);
    }

    std::vector<Vector> nodes_loc;
    nodes_loc.push_back(x_0);
    nodes_loc.push_back(x_1);
    nodes_loc.push_back(x_2);
    nodes_loc.push_back(x_3);

    // Must pass interpolator through a dummy integrator
    Integrator<InterpolType> integ(intp);
    ParticleSet<Integrator<InterpolType>> pset_loc(integ, prm.Nrw_max, this->m_mpi);
    pset_loc.add(nodes_loc, 0);

    this->edges.push_back({{0, 1}, dist(nodes_loc[0], nodes_loc[1])});
    this->edges.push_back({{1, 2}, dist(nodes_loc[1], nodes_loc[2])});
    this->edges.push_back({{2, 0}, dist(nodes_loc[2], nodes_loc[0])});
    this->edges.push_back({{1, 3}, dist(nodes_loc[1], nodes_loc[3])});
    this->edges.push_back({{2, 3}, dist(nodes_loc[2], nodes_loc[3])});
    this->edges.push_back({{0, 3}, dist(nodes_loc[0], nodes_loc[3])});

    this->faces.push_back({{0, 1, 2}, 1.});
    this->faces.push_back({{0, 3, 5}, 1.});
    this->faces.push_back({{1, 4, 3}, 1.});
    this->faces.push_back({{2, 5, 4}, 1.});

    NodesListType nodes_inlet_dummy;
    EdgesListType edges_inlet_dummy;
    compute_edge2faces(edge2faces_loc, this->faces, this->edges);
    compute_node2edges(node2edges_loc, this->edges, pset_loc.N());

    Uint n_add, n_rem;
    do {
      n_add = sheet_refinement(this->faces, this->edges, edge2faces_loc, node2edges_loc, edges_inlet_dummy,
                               pset_loc, prm.ds_max, 0.0, false);
      for (Uint irw=0; irw<pset_loc.N(); ++irw){
        Vector x = pset_loc.x(irw);
        Vector nn = x / x.norm();
        Real rad = 1./sqrt(nn[0]*nn[0]/lx2 + nn[1]*nn[1]/ly2 + nn[2]*nn[2]/lz2);
        pset_loc.set_x(irw, rad * nn);
      }
      n_rem = sheet_coarsening(this->faces, this->edges, edge2faces_loc, node2edges_loc, edges_inlet_dummy, nodes_inlet_dummy,
                               pset_loc, prm.ds_min, 0.0);

      std::cout << "Added " << n_add << " and removed " << n_rem << " edges." << std::endl;
    } while (n_add > 0 || n_rem > 0);

    for (Uint iedge=0; iedge<this->edges.size(); ++iedge){
      Uint inode = this->edges[iedge].first[0];
      Uint jnode = this->edges[iedge].first[1];
      this->edges[iedge].second = pset_loc.dist(inode, jnode);
    }
    for (Uint iface=0; iface<this->faces.size(); ++iface){
      Uint iedge = this->faces[iface].first[0];
      Uint jedge = this->faces[iface].first[1];
      this->faces[iface].second = pset_loc.triangle_area(iedge, jedge, this->edges);
    }
    for (Uint irw=0; irw<pset_loc.N(); ++irw){
      this->nodes.push_back(pset_loc.x(irw));
    }
  };
};
*/
class RandomPairsInitializer : public Initializer {
protected:
  std::mt19937 &gen;
  std::vector<std::string> key;
public:
  RandomPairsInitializer( const std::vector<std::string>& key
                        //, IntpType& intp
                        , Parameters& prm
                        , MPIwrap& mpi
                        , std::mt19937 &gen
                        )
   //: Initializer(intp, prm, mpi), gen(gen) {
    : Initializer(prm, mpi), gen(gen), key(key) {
      //probe(intp);
  };
  template<typename IntpType>
  void probe(IntpType& intp);
};

template<typename IntpType>
void RandomPairsInitializer::probe(IntpType& intp){
  x_min = intp.get_x_min();
  x_max = intp.get_x_max();
  L = x_max - x_min;

  std::uniform_real_distribution<> uni_dist_x(x_min[0], x_max[0]);
  std::uniform_real_distribution<> uni_dist_y(x_min[1], x_max[1]);
  std::uniform_real_distribution<> uni_dist_z(x_min[2], x_max[2]);
  std::normal_distribution<Real> rnd_normal(0.0, 1.0);

  Uint Npairs = (key[0] == "pair") ? 1 : prm.Nrw/2;

  std::cout << "Npairs = " << Npairs << std::endl;

  Vector x0_ = this->x0;
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
    Vector dx(0., 0., 0.);
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

    Vector x_a = x0_ + dx;
    intp.probe(x_a);
    bool inside_a = intp.inside_domain();
    Vector x_b = x0_ - dx;
    intp.probe(x_b);
    bool inside_b = intp.inside_domain();
    if (inside_a && inside_b){
      //std::cout << "INSIDE" << std::endl;
      // std::cout << "INSIDE: " << x_a << " " << x_b << std::endl; 
      this->nodes.push_back(x_a);
      this->nodes.push_back(x_b);
      Real ds0 = (x_a-x_b).norm();
      this->edges.push_back({{2*ipair, 2*ipair+1}, ds0});
      ++ipair;
    }
    else if (key[0] == "pair"){
      std::cout << "Pair not inside domain" << std::endl;
      exit(0);
    }
    //std::cout << x_a << " " << x_b << std::endl;
  }
};

class RandomPointsInitializer : public Initializer {
protected:
  std::mt19937 &gen;
  std::vector<std::string> key;
public:
  RandomPointsInitializer( const std::vector<std::string>& key
                         //, IntpType& intp
                         , Parameters& prm
                         , MPIwrap& mpi
                         , std::mt19937 &gen
                         )
   //: Initializer(intp, prm, mpi), gen(gen){
    : Initializer(prm, mpi), gen(gen), key(key) {
      //probe(intp);
  };
  ~RandomPointsInitializer() { std::cout << "Destruct Initializer." << std::endl; };
  template<typename IntpType>
  void probe(IntpType& intp);
};

template<typename IntpType>
void RandomPointsInitializer::probe(IntpType& intp){
  x_min = intp.get_x_min();
  x_max = intp.get_x_max();
  L = x_max - x_min;

  std::uniform_real_distribution<> uni_dist_x(x_min[0], x_max[0]);
  std::uniform_real_distribution<> uni_dist_y(x_min[1], x_max[1]);
  std::uniform_real_distribution<> uni_dist_z(x_min[2], x_max[2]);

  Uint Nrw = prm.Nrw;

  Vector x0_ = this->x0;
  Uint irw=0;
  while (irw < Nrw){
    if (key[0] == "points" && key.size() == 2){
      if (contains(key[1], "x")){
        x0_[0] = uni_dist_x(gen);
      }
      if (contains(key[1], "y")){
        x0_[1] = uni_dist_y(gen);
      }
      if (contains(key[1], "z")){
        x0_[2] = uni_dist_z(gen);
      }
    }
    
    intp.probe(x0_);
    bool inside = intp.inside_domain();
    if (inside){
      this->nodes.push_back(x0_);
      ++irw;
    }
    else if (key[0] == "point"){
      std::cout << "Point not inside domain" << std::endl;
      exit(0);
    }
  }
};

class RandomGaussianStripInitializer : public Initializer {
protected:
  std::mt19937 &gen;
  std::vector<std::string> key;
public:
  RandomGaussianStripInitializer( const std::vector<std::string>& key
                                //, std::shared_ptr<Interpol> intp
                                , Parameters& prm, MPIwrap& mpi
                                , std::mt19937 &gen
                                ) : Initializer(prm, mpi), gen(gen), key(key) {};
  ~RandomGaussianStripInitializer(){ std::cout << "Destructing initializer!" << std::endl; };
  template<typename IntpType>
  void probe(IntpType& intp);
};

template<typename IntpType>
void RandomGaussianStripInitializer::probe(IntpType& intp){
  this->edges.clear();
  this->faces.clear();

  double La = prm.La;
  double sigma0 = prm.Lb;
  
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
  while (irw < prm.Nrw && failed_attempts < max_failed_attempts){
    double alpha = rnd_unit(gen);
    Vector3d xi = alpha * x00 + (1.-alpha) * x01;
    if (init_rand_x)
      xi[0] += sigma0 * rnd_normal(gen);
    if (init_rand_y)
      xi[1] += sigma0 * rnd_normal(gen);
    if (init_rand_z)
      xi[2] += sigma0 * rnd_normal(gen);
    // check if inside domain
    intp.probe(xi);
    if (intp.inside_domain()){
      this->nodes.push_back(xi);
      ++irw;
      failed_attempts = 0;
    }
    else {
      ++failed_attempts;
    }
  }
  if (irw == 0) {
    std::cout << "No points inside domain" << std::endl;
    exit(0);
  }
  prm.Nrw = irw;
};

class RandomGaussianCircleInitializer : public Initializer {
protected:
  std::mt19937 &gen;
  std::vector<std::string> key;
public:
  RandomGaussianCircleInitializer( const std::vector<std::string>& key
                                //, std::shared_ptr<Interpol> intp
                                , Parameters& prm, MPIwrap& mpi
                                , std::mt19937 &gen
                                ) : Initializer(prm, mpi), gen(gen), key(key) {};
  ~RandomGaussianCircleInitializer(){ std::cout << "Destructing initializer!" << std::endl; };
  template<typename IntpType>
  void probe(IntpType& intp);
};


template<typename IntpType>
void RandomGaussianCircleInitializer::probe(IntpType& intp){
  edges.clear();
  faces.clear();

  double R = prm.La/2;
  double sigma0 = prm.Lb;
  
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

  std::vector<bool> init_rand_x = {contains(key[2], "x"),
                                   contains(key[2], "y"), 
                                   contains(key[2], "z")};

  std::uniform_real_distribution<double> rnd_unit(0.0, 1.0);
  std::normal_distribution<double> rnd_normal(0.0, 1.0);

  Uint failed_attempts = 0;
  Uint max_failed_attempts = 1000000; // Maybe not hardcode?

  Uint irw = 0;
  while (irw < prm.Nrw && failed_attempts < max_failed_attempts){
    double alpha1 = 1.;
    double alpha2 = 1.;
    while (pow(alpha1, 2) + pow(alpha2, 2) > 1){
      alpha1 = 2*rnd_unit(gen)-1;
      alpha2 = 2*rnd_unit(gen)-1;
    }
    
    Vector3d xi = x0 + R * (alpha1 * t1 + alpha2 * t2);
    for (Uint dim=0; dim<3; ++dim)
      if (init_rand_x[dim])
        xi[dim] += sigma0 * rnd_normal(gen);

    // check if inside domain
    intp.probe(xi);
    if (intp.inside_domain()){
      this->nodes.push_back(xi);
      ++irw;
      failed_attempts = 0;
    }
    else {
      ++failed_attempts;
    }
  }
  if (irw == 0) {
    std::cout << "No points inside domain" << std::endl;
    exit(0);
  }
  prm.Nrw = irw;
};


/*{
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
    Real tol = 1e-12;
    Uint N_est = 10000000;
    Real dx_est;
  
    Real Lx = this->L[0];
    Real Ly = this->L[1];
    Real Lz = this->L[2];

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
    Real dx = Lx/Nx;
    Real dy = Ly/Ny;
    Real dz = Lz/Nz;

    Real ww;

    std::vector<Real> wei;
    std::vector<Vector> pos;
    for (Uint ix=0; ix<Nx; ++ix){
      for (Uint iy=0; iy<Ny; ++iy){
        for (Uint iz=0; iz<Nz; ++iz){
          Vector x = this->x0;
          if (init_rand_x) x[0] = this->x_min[0]+(ix+0.5)*dx;
          if (init_rand_y) x[1] = this->x_min[1]+(iy+0.5)*dy;
          if (init_rand_z) x[2] = this->x_min[2]+(iz+0.5)*dz;
          intp.probe(x);
          if (prm.init_weight == "ux"){
            ww = abs(intp.get_ux());
          }
          else if (prm.init_weight == "uy"){
            ww = abs(intp.get_uy());
          }
          else if (prm.init_weight == "uz"){
            ww = abs(intp.get_uz());
          }
          else if (prm.init_weight == "u"){
            ww = sqrt(pow(intp.get_ux(), 2)
                      + pow(intp.get_uy(), 2)
                      + pow(intp.get_uz(), 2));
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
      Vector x;
      do {
        Uint ind = discrete_dist(gen);
        x = pos[ind];
        if (init_rand_x) x[0] += uni_dist_dx(gen);
        if (init_rand_y) x[1] += uni_dist_dy(gen);
        if (init_rand_z) x[2] += uni_dist_dz(gen);
        intp.probe(x);
        //std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
      } while (!intp.inside_domain());
      this->nodes.push_back(x);
    }

    sort(this->nodes.begin(), this->nodes.end(), less_than_op());
    
    //for (Uint irw=1; irw < prm.Nrw; ++irw){
    //  Real ds0 = dist(this->nodes[irw-1], this->nodes[irw]);
    //  if (ds0 < 10*prm.ds_init)  // 2 lattice units (before) --> 10 x ds_max (now)
    //    this->edges.push_back({{irw-1, irw}, ds0});
    //  // Needs customization for 2D/3D applications
    //}
  };
};/*


/*std::vector<Vector> initial_positions(const std::string init_mode,
                                        const std::string init_weight,
                                        Uint &Nrw,
                                        const Vector &x0,
                                        const Real La,
                                        const Real Lb,
                                        const Real ds,
                                        const Real t0,
                                        Interpol *intp,
                                        std::mt19937 &gen,
                                        EdgesType &edges,
                                        FacesType &faces
                                   ){
  intp.update(t0);

  Vector x;
  Real Lx = intp.get_Lx();
  Real Ly = intp.get_Ly();
  Real Lz = intp.get_Lz();
  Vector x_min = intp.get_x_min();
  Vector x_max = intp.get_x_max();

  Uint Nx = 1;
  Uint Ny = 1;
  Uint Nz = 1;

  std::vector<std::string> key = split_string(init_mode, "_");

  if (key.size() == 0){
    std::cout << "init_mode not specified." << std::endl;
    exit(0);
  }

  if (key[0] == "point"){
    std::vector<Vector> pos_init;
    intp.probe(x0);
    for (Uint irw=0; irw < Nrw; ++irw){
      if (intp.inside_domain()){
        pos_init.push_back(x0);
      }
    }
    edges.clear();
    return pos_init;
  }

  if (key[0] == "uniform"){
    std::vector<Vector> pos_init;

    Vector x_a = x0;
    Vector x_b = x0;
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
    Vector Dx = (x_b - x_a) / (Nrw-1);
    for (Uint irw=0; irw < Nrw; ++irw){
      x = x_a + Dx * irw;
      intp.probe(x);
      if (intp.inside_domain()){
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
    std::vector<Vector> pos_init;

    Vector n(0., 0., 0.);
    Vector ta(0., 0., 0.);
    Vector tb(0., 0., 0.);
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

    Vector x00 = x0;
    x00[0] += - La/2*ta[0] - Lb/2*tb[0];
    x00[1] += - La/2*ta[1] - Lb/2*tb[1];
    x00[2] += - La/2*ta[2] - Lb/2*tb[2];

    Vector x01 = x0;
    x01[0] += La/2*ta[0] + Lb/2*tb[0];
    x01[1] += La/2*ta[1] - Lb/2*tb[1];
    x01[2] += - La/2*ta[2] - Lb/2*tb[2];

    Vector x10 = x0;
    x10[0] += La/2*ta[0] + Lb/2*tb[0];
    x10[1] += La/2*ta[1] + Lb/2*tb[1];
    x10[2] += La/2*ta[2] + Lb/2*tb[2];

    Vector x11 = x0;
    x11[0] += - La/2*ta[0] - Lb/2*tb[0];
    x11[1] += - La/2*ta[1] + Lb/2*tb[1];
    x11[2] += - La/2*ta[2] + Lb/2*tb[2];

    intp.probe(x00);
    bool inside_00 = intp.inside_domain();
    intp.probe(x01);
    bool inside_01 = intp.inside_domain();
    intp.probe(x10);
    bool inside_10 = intp.inside_domain();
    intp.probe(x11);
    bool inside_11 = intp.inside_domain();
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
    std::vector<Vector> pos_init;

    std::uniform_real_distribution<> uni_dist_x(x_min[0], x_max[0]);
    std::uniform_real_distribution<> uni_dist_y(x_min[1], x_max[1]);
    std::uniform_real_distribution<> uni_dist_z(x_min[2], x_max[2]);
    std::normal_distribution<Real> rnd_normal(0.0, 1.0);

    //std::cout << "x_min = " << x_min << std::endl;
    //std::cout << "x_max = " << x_max << std::endl;

    // std::uniform_real_distribution<> uni_dist_theta(0., 2*M_PI);
    Uint Npairs = (key[0] == "pair") ? 1 : Nrw/2;

    std::cout << "Npairs = " << Npairs << std::endl;

    Vector x0_ = x0;
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
      Vector dx(0., 0., 0.);
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

      Vector x_a = x0_ + dx;
      intp.probe(x_a);
      bool inside_a = intp.inside_domain();
      Vector x_b = x0_ - dx;
      intp.probe(x_b);
      bool inside_b = intp.inside_domain();
      if (inside_a && inside_b){
        //std::cout << "INSIDE" << std::endl;
        // std::cout << "INSIDE: " << x_a << " " << x_b << std::endl; 
        pos_init.push_back(x_a);
        pos_init.push_back(x_b);
        Real ds0 = (x_a-x_b).norm();
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
  Real tol = 1e-12;
  Uint N_est = 1000000;
  Real dx_est;
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
  Real dx = Lx/Nx;
  Real dy = Ly/Ny;
  Real dz = Lz/Nz;

  Real ww;

  std::vector<Real> wei;
  std::vector<Vector> pos;
  for (Uint ix=0; ix<Nx; ++ix){
    for (Uint iy=0; iy<Ny; ++iy){
      for (Uint iz=0; iz<Nz; ++iz){
        x = x0;
        if (init_rand_x) x[0] = x_min[0]+(ix+0.5)*dx;
        if (init_rand_y) x[1] = x_min[1]+(iy+0.5)*dy;
        if (init_rand_z) x[2] = x_min[2]+(iz+0.5)*dz;
        intp.probe(x);
        if (init_weight == "ux"){
          ww = abs(intp.get_ux());
        }
        else if (init_weight == "uy"){
          ww = abs(intp.get_uy());
        }
        else if (init_weight == "uz"){
          ww = abs(intp.get_uz());
        }
        else if (init_weight == "u"){
          ww = sqrt(pow(intp.get_ux(), 2)
                    + pow(intp.get_uy(), 2)
                    + pow(intp.get_uz(), 2));
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

  std::vector<Vector> pos_init;
  for (Uint irw=0; irw<Nrw; ++irw){
    do {
      Uint ind = discrete_dist(gen);
      x = pos[ind];
      if (init_rand_x) x[0] += uni_dist_dx(gen);
      if (init_rand_y) x[1] += uni_dist_dy(gen);
      if (init_rand_z) x[2] += uni_dist_dz(gen);
      intp.probe(x);
      //std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
    } while (!intp.inside_domain());
    pos_init.push_back(x);
  }

  sort(pos_init.begin(), pos_init.end(), less_than_op());

  for (Uint irw=1; irw < Nrw; ++irw){
    Real ds0 = dist(pos_init[irw-1], pos_init[irw]);
    if (ds0 < 10*ds)  // 2 lattice units (before) --> 10 x ds_max (now)
      edges.push_back({{irw-1, irw}, ds0});
    // Needs customization for 2D/3D applications
  }

  return pos_init;
}*/


#endif
