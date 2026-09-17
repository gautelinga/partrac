#ifndef __TOPOLOGY_HPP
#define __TOPOLOGY_HPP

#include <fstream>
#include <map>
#include <memory>
#include <random>
#include <string>
#include "H5Cpp.h"
#include "typedefs.hpp"
#include "Params.hpp"
#include "ParticleSet.hpp"
#include "stats.hpp"

class Initializer;

class Topology {
public:
  Topology(ParticleSet& ps, const partrac::Params& prm);
  int dim();
  int dim_settled();
  void check_topology();
  void check_dim();
  void compute_maps();
  void clear();
  Uint refine();
  Uint coarsen(const bool full);
  bool filter();
  Uint remove_beyond(const int, const double);
  Uint inject();
  void compute_interior();
  // Whether compute_interior fills H and n
  bool computes_curvature() const { return curv_refine_factor > 0.; };
  void remove_nodes_safe(std::vector<bool>&);
  bool resize(const double);
  bool resize_doublings(const double);
  void integrate_tau(const double, const double);
  EdgesType edges;
  FacesType faces;
  Edge2FacesType edge2faces;
  Node2EdgesType node2edges;
  std::vector<Vector3d> pos_inj;
  int dim0 = -1;   // -1 until set
  double Dm = 0.;  // not every app has an integrator
  EdgesType edges_inj;
  EdgesListType edges_inlet;
  NodesListType nodes_inlet;
  // Initial curvature?
  InteriorAnglesType interior_ang;
  std::vector<double> mixed_areas;
  std::vector<Vector3d> face_normals;
  // Sort particles by cell and renumber
  bool sort_by_cell();
  // Shuffle particles and renumber
  void shuffle(std::mt19937& gen);
  // Checkpoint t_loc
  bool records_t_loc = false;
  // Edge doublings (filaments)
  bool records_doublings = false;
  std::vector<Uint> doublings;
  std::vector<Uint>& edge_doublings(){ doublings.resize(edges.size(), 0); return doublings; }
  void write_checkpoint(const std::string& checkpointsfolder, const double t, partrac::Params& prm) const;
  void load_checkpoint(const std::string& checkpointsfolder, const partrac::Params& prm);
  void dump_hdf5(H5::H5File& h5f, const std::string& groupname, std::map<std::string, bool>& output_fields);
  void load_initial_state(std::shared_ptr<Initializer> init_state, partrac::Params& prm);
  // Mesh statistics
  std::vector<StatsColumn> stats_header_columns(const double ds_max){
    return mesh_stats_columns(0., ps, faces, edges, ds_max, 0, 0, dim_settled());
  }
  template<typename T>
  void write_statistics(std::ofstream &statfile, const double t, const double ds_max, //const bool do_dump_hist, const std::string histfolder, 
                        T& integrator);
                        //std::shared_ptr<Integrator> integrator);
private:
  ParticleSet& ps;
  void renumber(const std::vector<Uint>& old2new);
  double ds_min;
  double ds_max;
  double curv_refine_factor;
  bool cut_if_stuck;
  bool injecting;
  bool inject_edges;
  bool verbose;
  Uint filter_target;
};

template<typename T>
void Topology::write_statistics( std::ofstream &statfile
                               , const double t
                               , const double ds_max
                               , T& integrator){
                           //const bool do_dump_hist,
                           //const std::string histfolder,
                           //std::shared_ptr<Integrator> integrator){
  write_stats_row(statfile, mesh_stats_columns(t, ps, faces, edges, ds_max,
                                               integrator.get_accepted(),
                                               integrator.get_declined(),
                                               dim_settled()));
              //integrator->get_accepted(), integrator->get_declined());
}

#endif
