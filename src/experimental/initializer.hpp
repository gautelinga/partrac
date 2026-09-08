#ifndef __EXP_INITIALIZER_HPP
#define __EXP_INITIALIZER_HPP

#include <iomanip>
#include <vector>
#include <random>
#include <algorithm>

#include "typedefs.hpp"
#include "../Interpol.hpp"
#include "utils.hpp"
#include "../Params.hpp"
//#include "particles.hpp"

// Parameters read by the experimental initializers.
// As in src/Initializer.hpp, but the gaussian circle here also reads key[2].
// Exactly the tokens each shape reads; a longer init_mode would be truncated.
inline bool experimental_init_mode_shape_ok(const std::string& init_mode){
  const std::vector<std::string> key = split_string(init_mode, "_");
  if (!init_mode_dirs_ok(key)) return false;
  if (key[0] == "point") return key.size() == 1;
  if (key[0] == "pairs") return key.size() == 2 || key.size() == 3;
  // the apps dispatch on a substring, so plain strip is the gaussian one too
  if (contains(key[0], "strip") || contains(key[0], "circle"))
    return key.size() == 3;
  return key.size() == 2;
}

inline void add_experimental_initializer_params(partrac::Schema& s){
  s.require<std::string>("init_mode", "initial distribution");
  s.require<Uint>("Nrw", "number of particles");
  // Nrw stays the request; these record what happened
  s.runtime<Uint>("Nrw_init", 0, "particles the initializer placed");
  s.runtime<Uint>("Nrw_current", 0, "particles in the set when this was written");
  s.opt<double>("x0", 0.0, "initial position");
  s.opt<double>("y0", 0.0, "initial position");
  s.opt<double>("z0", 0.0, "initial position");
  s.opt<bool>("inject", false, "inject new particles");
  s.opt<bool>("clear_initial_edges", false, "drop the initial edges");
  s.check([](const partrac::Params& p){
            return experimental_init_mode_shape_ok(p.get<std::string>("init_mode"));
          },
          "init_mode has the wrong number of directions: most modes take one,"
            " as in points_xy; point takes none; pairs takes one or two;"
            " randomgaussianstrip and randomgaussiancircle take two, as in"
            " randomgaussianstrip_x_y");
}

// Parameters every app in this family reads
inline void add_common_app_params(partrac::Schema& s){
  s.require<double>("Dm", "molecular diffusivity");
  s.require<double>("dt", "timestep");
  s.require<double>("T", "final time");
  s.require<Uint>("Nrw_max", "max number of particles");
  s.opt<double>("t0", 0.0, "start time");
  s.opt<double>("U", 1.0, "velocity scale");
  s.opt<double>("dump_intv", 100.0, "dump interval");
  s.opt<double>("stat_intv", 100.0, "statistics interval");
  s.opt<double>("checkpoint_intv", 1000.0, "checkpoint interval");
  s.opt<int>("seed", 0, "random seed");
  s.opt<int>("dump_chunk_size", 0, "particles per dump chunk");
  s.opt<bool>("minimal_output", false, "dump less");
  s.opt<bool>("verbose", false, "print the parameters");
  s.opt<std::string>("tag", "", "appended to the folder name");
  // dump_intv and stat_intv become step counts, so they must not round to zero
  // an interval of 0 turns that output off; a negative one is a typo
  s.check([](const partrac::Params& p){
            for (const auto& key : {"checkpoint_intv", "dump_intv", "stat_intv"})
              if (p.get<double>(key) < 0.) return false;
            return true;
          },
          "intervals cannot be negative");
  s.finalize([](partrac::Params& p){
    const double dt = p.get<double>("dt");
    if (p.get<double>("dump_intv") > 0.)
      p.set<double>("dump_intv", std::max(p.get<double>("dump_intv"), dt));
    if (p.get<double>("stat_intv") > 0.)
      p.set<double>("stat_intv", std::max(p.get<double>("stat_intv"), dt));
    p.set<Uint>("Nrw_max", std::max(p.get<Uint>("Nrw_max"), p.get<Uint>("Nrw")));
  });
}

// Parameters used by the restart path
inline void add_restart_params(partrac::Schema& s){
  s.opt<bool>("random", true, "draw the seed randomly");
  s.opt<bool>("output_all_props", true, "dump all properties");
  s.opt<std::string>("restart_folder", "", "folder to restart from");
  s.runtime<std::string>("folder", "", "output folder");
  s.runtime<double>("t", 0.0, "current time");
}

namespace experimental {

// TODO: Massive cleanup!
struct less_than_op {
  inline bool operator() (const Vector &a, const Vector &b){
    return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]) || (a[0] == b[0] && a[1] == b[1] && a[2] < b[2]);
  }
};

class Initializer {
public:
  //Initializer(IntpType& intp, Parameters& prm) : m_intp(intp), prm(prm) {
  //Initializer(Parameters& prm) : prm(prm) { // {
  Initializer(partrac::Params& prm) : prm(prm) {
    x0 = {prm.get<double>("x0"), prm.get<double>("y0"), prm.get<double>("z0")};
    //x_min = intp.get_x_min();
    //x_max = intp.get_x_max();
    //L = x_max - x_min;
    inject = prm.get<bool>("inject");
    clear_initial_edges = prm.get<bool>("clear_initial_edges");
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
  partrac::Params& prm;
  Vector x0;
  Vector x_min;
  Vector x_max;
  Vector L;
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
  // Nrw is the request; these are what was placed
  prm.set<Uint>("Nrw_init", nodes.size());
  prm.set<Uint>("Nrw_current", nodes.size());
}


class UniformInitializer : public Initializer {
protected:  
  std::vector<std::string> key;
public:
  UniformInitializer( const std::vector<std::string>& key
                    //, IntpType& intp
                    , partrac::Params& prm
                    //
                    //) : Initializer(intp, prm) {
                    ) : Initializer(prm), key(key) {
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
    exit(1);
  }
  Vector Dx = (x_b - x_a) / (prm.get<Uint>("Nrw")-1);
  for (Uint irw=0; irw < prm.get<Uint>("Nrw"); ++irw){
    Vector x = x_a + Dx * irw;
    if (intp.locate(x)){
      this->nodes.push_back(x);
    }
  }
  for (Uint irw=0; irw < this->nodes.size()-1; ++irw){
    if ((this->nodes[irw] - this->nodes[irw+1]).norm() < 1.5*Dx.norm()){
      this->edges.push_back({{irw, irw+1}, dist(this->nodes[irw], this->nodes[irw+1])});
    }
  }
};

class RandomPairsInitializer : public Initializer {
protected:
  std::mt19937 &gen;
  std::vector<std::string> key;
public:
  RandomPairsInitializer( const std::vector<std::string>& key
                        //, IntpType& intp
                        , partrac::Params& prm
                        //
                        , std::mt19937 &gen
                        )
   //: Initializer(intp, prm), gen(gen) {
    : Initializer(prm), gen(gen), key(key) {
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

  Uint Npairs = (key[0] == "pair") ? 1 : prm.get<Uint>("Nrw")/2;

  std::cout << "Npairs = " << Npairs << std::endl;

  Vector x0_ = this->x0;
  // only this shape redraws the centre; any other keeps x0, so a centre
  // outside the domain can never yield a pair
  const bool centre_moves = (key[0] == "pairs" && key.size() == 3);
  if (!centre_moves && !intp.locate(x0_)){
    std::cout << "Pair centre is not inside the domain" << std::endl;
    exit(1);
  }
  Uint ipair=0;
  Uint failed_attempts = 0;
  Uint max_failed_attempts = 1000000; // as in the gaussian initializers
  while (ipair < Npairs && failed_attempts < max_failed_attempts){
    if (centre_moves){
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
    dx *= 0.5*prm.get<double>("ds_init")/dx.norm();

    Vector x_a = x0_ + dx;
    bool inside_a = intp.locate(x_a);
    Vector x_b = x0_ - dx;
    bool inside_b = intp.locate(x_b);
    if (inside_a && inside_b){
      //std::cout << "INSIDE" << std::endl;
      // std::cout << "INSIDE: " << x_a << " " << x_b << std::endl; 
      this->nodes.push_back(x_a);
      this->nodes.push_back(x_b);
      Real ds0 = (x_a-x_b).norm();
      this->edges.push_back({{2*ipair, 2*ipair+1}, ds0});
      ++ipair;
      failed_attempts = 0;
    }
    else if (key[0] == "pair"){
      std::cout << "Pair not inside domain" << std::endl;
      exit(1);
    }
    else {
      ++failed_attempts;
    }
    //std::cout << x_a << " " << x_b << std::endl;
  }
  if (ipair < Npairs){
    std::cout << "Could not place all pairs inside the domain" << std::endl;
    exit(1);
  }
};

class RandomPointsInitializer : public Initializer {
protected:
  std::mt19937 &gen;
  std::vector<std::string> key;
public:
  RandomPointsInitializer( const std::vector<std::string>& key
                         //, IntpType& intp
                         , partrac::Params& prm
                         //
                         , std::mt19937 &gen
                         )
   //: Initializer(intp, prm), gen(gen){
    : Initializer(prm), gen(gen), key(key) {
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

  Uint Nrw = prm.get<Uint>("Nrw");

  Vector x0_ = this->x0;
  Uint irw=0;
  while (irw < Nrw){
    if (key[0] == "points"){ // && key.size() == 2){
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
    
    bool inside = intp.locate(x0_);
    if (inside){
      this->nodes.push_back(x0_);
      ++irw;
    }
    else if (key[0] == "point"){
      std::cout << "Point not inside domain" << std::endl;
      exit(1);
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
                                , partrac::Params& prm
                                , std::mt19937 &gen
                                ) : Initializer(prm), gen(gen), key(key) {};
  ~RandomGaussianStripInitializer(){ std::cout << "Destructing initializer!" << std::endl; };
  template<typename IntpType>
  void probe(IntpType& intp);
};

template<typename IntpType>
void RandomGaussianStripInitializer::probe(IntpType& intp){
  this->edges.clear();
  this->faces.clear();

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
    if (intp.locate(xi)){
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
    exit(1);
  }
};

class RandomGaussianCircleInitializer : public Initializer {
protected:
  std::mt19937 &gen;
  std::vector<std::string> key;
public:
  RandomGaussianCircleInitializer( const std::vector<std::string>& key
                                //, std::shared_ptr<Interpol> intp
                                , partrac::Params& prm
                                , std::mt19937 &gen
                                ) : Initializer(prm), gen(gen), key(key) {};
  ~RandomGaussianCircleInitializer(){ std::cout << "Destructing initializer!" << std::endl; };
  template<typename IntpType>
  void probe(IntpType& intp);
};


template<typename IntpType>
void RandomGaussianCircleInitializer::probe(IntpType& intp){
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

  std::vector<bool> init_rand_x = {contains(key[2], "x"),
                                   contains(key[2], "y"), 
                                   contains(key[2], "z")};

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
      if (init_rand_x[dim])
        xi[dim] += sigma0 * rnd_normal(gen);

    // check if inside domain
    if (intp.locate(xi)){
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
    exit(1);
  }
};




}

#endif
