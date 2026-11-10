#ifndef __MESH_HPP
#define __MESH_HPP

#include <algorithm>
#include <string>
#include <vector>
#include "H5Cpp.h"
#include "typedefs.hpp"
#include "ParticleSet.hpp"

void mesh2hdf( H5::H5File& h5f, const std::string& groupname
             , const ParticleSet& ps
             , const FacesType& faces
             , const EdgesType& edges
             , const bool output_tau
             , const std::vector<Uint>* doublings = nullptr
             );

Uint get_common_entry(Uint kedge, Uint ledge,
                      EdgesType &edges);

// Edge lengths shrink where a strip bends: by 1/min(1 + factor kappa, cap); no limit with factor 0. Not on sheets
// (Topology::check_dim refuses): splits alone cannot grade a sheet that steeply
struct CurvLimit {
  double factor = 0.;
  double cap = 16.;
  double scale(const double kappa) const { return 1./std::min(1. + factor*kappa, cap); }
};

// An edge whose midpoint cannot be placed: kept as it is, cut, or the run stops
enum class StuckEdge { Keep, Cut, Stop };

// Inlet edges are split along with their template
Uint sheet_refinement(Connectivity& m,
                      ParticleSet& ps,
                      const double ds_max,
                      const StuckEdge stuck,
                      const bool check_if_inside=true);

Uint refinement(Connectivity& m,
                ParticleSet& ps,
                const double ds_max,
                const CurvLimit& curv,
                const StuckEdge stuck);

void compute_edge2faces(Edge2FacesType &edge2faces,
                        const FacesType &faces,
                        const EdgesType &edges);

void compute_node2edges(Node2EdgesType &node2edges,
                        const EdgesType &edges,
                        const Uint Nrw);

// Remove inactive entities and whatever they leave dangling
void remove_inactive(Connectivity& m,
                     std::vector<bool> &face_isactive,
                     std::vector<bool> &edge_isactive,
                     std::vector<bool> &node_isactive,
                     ParticleSet& ps);

// ds_max > 0: no collapse lengthens an edge past it
Uint sheet_coarsening(Connectivity& m,
                      ParticleSet& ps,
                      const double ds_min,
                      const double ds_max);

Uint coarsening(Connectivity& m,
                ParticleSet& ps,
                const double ds_min,
                const double ds_max,
                const CurvLimit& curv);

// Thin a strip to filter_target edges
bool strip_filtering(Connectivity& m,
                     ParticleSet& ps,
                     const Uint filter_target);

bool resizing(EdgesType &edges,
              Node2EdgesType &node2edges,
              ParticleSet& ps,
              const double ds);

// Halve long edges to ds; count halvings per edge
bool resizing_doublings(const EdgesType &edges, std::vector<Uint>& doublings, ParticleSet& ps, const double ds);

void compute_interior_prop(InteriorAnglesType &interior_ang,
                           std::vector<double> &mixed_areas,
                           std::vector<Vector3d> &face_normals,
                           const Connectivity& m,
                           ParticleSet& ps);

void compute_mean_curv(const Connectivity& m,
                       ParticleSet& ps,
                       const InteriorAnglesType &interior_ang,
                       const std::vector<double> &mixed_areas,
                       const std::vector<Vector3d> &face_normals);

// Inject a generation; with inject_edges, stitch it to the previous one
bool injection(Connectivity& m,
               ParticleSet &ps,
               const bool inject_edges,
               const bool verbose);

#endif
