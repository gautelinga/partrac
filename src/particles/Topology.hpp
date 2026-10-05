#ifndef __TOPOLOGY_HPP
#define __TOPOLOGY_HPP

#include <fstream>
#include <memory>
#include <random>
#include <string>
#include "H5Cpp.h"
#include "typedefs.hpp"
#include "Params.hpp"
#include "ParticleSet.hpp"
#include "stats.hpp"

struct InitialState;

// A Topology's settings; a cloud leaves the remeshing ones at their defaults
struct TopologyOptions {
  double ds_min = 0.;
  double ds_max = 0.;
  double curv_refine_factor = 0.;
  double Dm = 0.;
  bool cut_if_stuck = true;
  bool injecting = false;
  bool inject_edges = true;
  bool verbose = false;
  Uint filter_target = 0;
};
// A cloud's options: lengths, injection, output, diffusivity
TopologyOptions cloud_options(const partrac::Params& prm);
// A line's or sheet's: the cloud's and the remeshing keys
TopologyOptions mesh_options(const partrac::Params& prm);

// The particles' connectivity, and what remeshes it
class Topology : public Connectivity {
public:
  Topology(ParticleSet& ps, const TopologyOptions& opts);
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
  bool computes_curvature() const { return opts.curv_refine_factor > 0.; };
  void remove_nodes_safe(std::vector<bool>&);
  bool resize(const double);
  bool resize_doublings(const double);
  void integrate_tau(const double, const double);
  int dim0 = -1;   // -1 until set
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
  void load_text_checkpoint(const std::string& checkpointsfolder, const partrac::Params& prm);
  bool has_t_loc(const partrac::Params& prm) const;   // t_loc in a checkpoint
  void dump_hdf5(H5::H5File& h5f, const std::string& groupname, const OutputFields& output_fields);
  void load_initial_state(const InitialState& init_state, partrac::Params& prm);
  // Mesh statistics
  std::vector<StatsColumn> stats_header_columns(const double ds_max){
    return mesh_stats_columns(0., ps, faces, edges, ds_max, 0, 0, dim_settled());
  }
  template<typename T>
  void write_statistics(std::ofstream &statfile, const double t, const double ds_max,
                        T& integrator);
private:
  ParticleSet& ps;
  void renumber(const std::vector<Uint>& old2new);
  TopologyOptions opts;
};

template<typename T>
void Topology::write_statistics( std::ofstream &statfile
                               , const double t
                               , const double ds_max
                               , T& integrator){
  write_stats_row(statfile, mesh_stats_columns(t, ps, faces, edges, ds_max,
                                               integrator.get_accepted(),
                                               integrator.get_declined(),
                                               dim_settled()));
}

#endif
