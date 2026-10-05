#ifndef __MESH_HPP
#define __MESH_HPP

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

// An edge whose midpoint cannot be placed: kept as it is, cut, or the run stops
enum class StuckEdge { Keep, Cut, Stop };

// Inlet edges are split along with their template
Uint sheet_refinement(Connectivity& m,
                      ParticleSet& ps,
                      const double ds_max,
                      const double curv_refine_factor,
                      const StuckEdge stuck,
                      const bool check_if_inside=true);

Uint refinement(Connectivity& m,
                ParticleSet& ps,
                const double ds_max,
                const double curv_refine_factor,
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

Uint sheet_coarsening(Connectivity& m,
                      ParticleSet& ps,
                      const double ds_min,
                      const double curv_refine_factor);

Uint coarsening(Connectivity& m,
                ParticleSet& ps,
                const double ds_min,
                const double curv_refine_factor);

bool filtering(Connectivity& m,
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
