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

// Inlet edges are split along with their template
Uint sheet_refinement(FacesType &faces,
                      EdgesType &edges,
                      Edge2FacesType &edge2faces,
                      Node2EdgesType &node2edges,                    
                      EdgesListType &edges_inlet,                
                      NodesListType &nodes_inlet,
                      std::vector<Vector3d> &pos_inj,
                      EdgesType &edges_inj,
                      ParticleSet& ps,
                      const double ds_max,
                      const double curv_refine_factor,
                      const bool cut_if_stuck,
                      const bool check_if_inside=true);

Uint refinement(FacesType &faces,
                EdgesType &edges,
                Edge2FacesType &edge2faces,
                Node2EdgesType &node2edges,
                EdgesListType &edges_inlet,
                NodesListType &nodes_inlet,
                std::vector<Vector3d> &pos_inj,
                EdgesType &edges_inj,
                ParticleSet& ps, const double ds_max,
                const double curv_refine_factor,
                const bool cut_if_stuck);

void compute_edge2faces(Edge2FacesType &edge2faces,
                        const FacesType &faces,
                        const EdgesType &edges);

void compute_node2edges(Node2EdgesType &node2edges,
                        const EdgesType &edges,
                        const Uint Nrw);

// Remove inactive entities and whatever they leave dangling
void remove_inactive(FacesType &faces, EdgesType &edges,
                     Edge2FacesType &edge2faces, Node2EdgesType &node2edges,
                     EdgesListType &edges_inlet, NodesListType &nodes_inlet,
                     std::vector<bool> &face_isactive,
                     std::vector<bool> &edge_isactive,
                     std::vector<bool> &node_isactive,
                     ParticleSet& ps);

Uint sheet_coarsening(FacesType &faces,
                      EdgesType &edges,
                      Edge2FacesType &edge2faces,
                      Node2EdgesType &node2edges,
                      EdgesListType& edges_inlet,
                      NodesListType& nodes_inlet,
                      ParticleSet& ps,
                      const double ds_min,
                      const double curv_refine_factor);

Uint coarsening(FacesType &faces,
                EdgesType &edges,
                Edge2FacesType &edge2faces,
                Node2EdgesType &node2edges,
                EdgesListType& edges_inlet,
                NodesListType& nodes_inlet,
                ParticleSet& ps,
                const double ds_min,
                const double curv_refine_factor);

bool filtering(FacesType &faces,
               EdgesType &edges,
               Edge2FacesType &edge2faces,
               Node2EdgesType &node2edges,
               EdgesListType &edges_inlet,
               NodesListType &nodes_inlet,
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
                           const FacesType &faces,
                           const EdgesType &edges,
                           const Edge2FacesType &edge2faces,
                           ParticleSet& ps);

void compute_mean_curv(const FacesType &faces,
                       const EdgesType &edges,
                       const Edge2FacesType &edge2faces,
                       const Node2EdgesType &node2edges,
                       ParticleSet& ps,
                       const InteriorAnglesType &interior_ang,
                       const std::vector<double> &mixed_areas,
                       const std::vector<Vector3d> &face_normals
                       );

// Inject a generation; with inject_edges, stitch it to the previous one
bool injection(const std::vector<Vector3d> &pos_inj,
               const EdgesType &edges_inj,
               EdgesListType &edges_inlet,
               NodesListType &nodes_inlet,
               EdgesType &edges,
               FacesType &faces,
               Edge2FacesType& edge2faces,
               Node2EdgesType& node2edges,
               ParticleSet &ps,
               const bool inject_edges,
               const bool verbose
               );

template<typename Row>
Uint get_common_entry(const Row& iedges, const Row& jedges){
  std::vector<Uint> out;
  std::set_intersection(iedges.begin(), iedges.end(), jedges.begin(), jedges.end(),
                        std::back_inserter(out));
  assert(out.size() == 1);
  return *out.begin();
}

// GL: Template this
template<typename T>
void print(const T vec){
  for ( auto & elem : vec ){
    std::cout << elem << " ";
  }
  std::cout << std::endl;
}

#endif
