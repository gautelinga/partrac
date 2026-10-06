#ifndef __SIMPLEX_LOAD_HPP
#define __SIMPLEX_LOAD_HPP

// Reading a dolfin HDF5 mesh and checkpointed function without dolfin: the
// mesh arrays with the cells in Morton order, the stored per-cell dof table
// composed through the global cell ids, and the values by node. Kept out of
// the translation unit that holds the evaluation.

#include <cstdint>
#include <string>
#include <vector>

#include "Params.hpp"
#include "cell_walk.hpp"
#include "mesh_tables.hpp"
#include "typedefs.hpp"

namespace simplex_load {

// No cell of the mesh carries this global id
constexpr std::uint32_t no_cell = 0xffffffffu;

// A stored Lagrange element: the degree, and the components of a vector one
struct Element {
  int degree = 0;
  std::size_t ncomp = 1;
};

// The mesh file's arrays, with the cells in Morton order of their centroids
struct MeshData {
  std::vector<std::uint32_t> topo;     // ncells x nv vertex ids into coords
  std::vector<double> coords;          // npoints x gdim
  std::vector<std::uint32_t> row_of;   // global cell id -> cell here, or no_cell
  std::size_t ncells = 0, nverts = 0;
  Uint gdim = 3;
  Vector3d x_min = Vector3d::Zero(), x_max = Vector3d::Zero();
};

// Nodes a cell of nv vertices carries at this degree
inline std::size_t nodes_per_cell(const int nv, const int degree){
  return degree == 1 ? std::size_t(nv) : std::size_t(nv) + std::size_t(nv*(nv-1)/2);
}

// A topology and a geometry dataset, with the cells renumbered into Morton
// order of their centroids; perm[new] = old is written for the caller. row_of
// stays empty: only a file that labels its cells globally fills it.
void read_mesh_arrays(const std::string& path, const std::string& topology,
                      const std::string& geometry, int nv, MeshData& m,
                      std::vector<std::uint32_t>& perm);

// The bounds of the coordinates, then the cells renumbered into Morton order
// of their centroids; perm[new] = old is written for the caller
void morton_cells(MeshData& m, int nv, std::vector<std::uint32_t>& perm);

// mesh/topology, mesh/coordinates and mesh/cell_indices, then the cells
// renumbered into Morton order; perm[new] = old is written for the caller
void read_mesh(const std::string& path, int nv, MeshData& m,
               std::vector<std::uint32_t>& perm);

// The element of <field> from its signature attribute
Element read_element(const std::string& path, const std::string& field, int nv);

// <field>'s stored dof table, one row per cell of m in m's order, and its vector
void read_field(const std::string& path, const std::string& field, const MeshData& m,
                const Element& el, int nv,
                std::vector<std::uint32_t>& rows, std::vector<double>& vec);

// A later stamp's vector_0 through the mapping the first stamp built, so every
// node a stored dof feeds -- its periodic images included -- gets the new value
void read_vector(const std::string& path, const std::string& field,
                 const mesh_tables::DofNodes& map, std::vector<double>& values);

// The nodes as points, the vertices first and the edge midpoints after them,
// stride doubles a node and the first gdim of them written
std::vector<double> node_points(const MeshData& m, const std::vector<std::uint32_t>& edges,
                                std::size_t nedges, int nv, std::size_t stride);

// A node and its periodic images are one node, whether or not the file was
// written from a space that identified them: the images are paired with their
// master here, and the master's dof serves them all.
struct NodePairs {
  std::vector<std::uint32_t> master;        // empty where no axis is periodic
  std::vector<double> points;               // held here when the mesh's own array is not it
  const std::vector<double>* x = nullptr;   // the nodes as points, gdim doubles a node
  Uint gdim = 0;
  Vector3d period = Vector3d::Zero();
  // From the mesh and the edges of a quadratic field, if any axis is periodic
  void build(const MeshData& m, const std::vector<std::uint32_t>& edges, std::size_t nedges,
             int nv, const std::vector<bool>& periodic, const Vector3d& x_min,
             const Vector3d& x_max, double tol);
  // An image node's stored value against its master's, to a relative tolerance
  void check(const std::vector<double>& values, std::size_t ncomp, const std::string& path,
             const std::string& field) const;
  // A field's n nodes with every image sent to its master, then renumbered by
  // node_map; empty where neither applies
  std::vector<std::uint32_t> masters(const std::vector<std::uint32_t>& node_map,
                                     std::size_t n) const;
};

// The element a parameter file declares; the file's own signature decides the
// element, but a name no loader knows is still a mistake
int declared_degree(const std::string& space, const char* what);

// What a parameter file says about the mesh and the spaces
struct Request {
  std::string infilename;       // the parameter file, named in the errors
  std::string mesh_file;
  std::string field_file;       // the stamp or component the elements come from
  std::string u_field, p_field, phi_field;
  const char* what = "";        // the loader, named in the errors
  int nv = 3;
  int want_u = 0, want_p = 0;   // the declared degrees
  bool include_pressure = true;
  bool include_phi = false;     // a phase field, in the element its file gives
  std::size_t n_dofs_max = 0;   // what an evaluation's buffers hold
  std::vector<bool> periodic = {false, false, false};
  double periodic_tol = 1e-12;
};

// Everything a dolfin HDF5 loader builds before its own field values: the mesh
// arrays, the fields' elements and node tables, the facet table and the
// mesh scale. The values are the loader's own: a stamp, or a frequency component.
struct Tables {
  MeshData mesh;
  int nv = 3;
  Element el_u, el_p, el_phi;
  Uint ncoeffs_u = 0, ncoeffs_p = 0, ncoeffs_phi = 0;
  std::vector<std::uint32_t> edges;   // empty unless a field carries midside nodes
  std::size_t nedges = 0;
  NodePairs np;
  // The nodes along the cells' curve, by degree; empty where the file's own order serves
  std::vector<std::uint32_t> map_quad, map_lin;
  const std::vector<std::uint32_t>& node_order(const Element& el) const {
    return el.degree == 2 ? map_quad : map_lin;
  }
  CellDofs u_dofs, p_dofs, phi_dofs;
  std::vector<std::int32_t> facets;
  double hmin = 0.;
};

// What the parameter file of a dolfin HDF5 loader says: the mesh path, the
// field names and the declared degrees. The caller has already set what its own
// object knows -- the file names, the cell, the periodicity and which fields it
// wants.
void request_from_params(Request& r, const partrac::Params& prm, const std::string& folder);

// read_mesh, the elements, then tables_from_mesh
void build_tables(const Request& r, Tables& t);

// The edge table, the node order, the periodic nodes, the node tables, the
// facet table and the mesh scale, in that order, of the mesh in t.mesh and the
// elements in t. A caller that knows the periodic images sets t.np.master and
// leaves r.periodic false: the node tables then read the masters it gave, and
// the facets on a periodic side stay walls for it to pair.
void tables_from_mesh(const Request& r, Tables& t);

// A later stamp or component into a buffer the size of the first one's, through
// the mapping the first built
inline void read_into(const std::string& path, const std::string& field,
                      const mesh_tables::DofNodes& map,
                      const std::vector<double>& first, std::vector<double>& values){
  values.assign(first.size(), 0.);
  read_vector(path, field, map, values);
}

// One field's values by node, and the mapping a later stamp or component is
// read through; node_map, when given, is this field's node order and the
// mapping comes back renumbered by it
void read_field_by_node(const std::string& path, const std::string& field, const Tables& t,
                        const Element& el, const std::vector<std::uint32_t>& node_map,
                        std::vector<double>& values, mesh_tables::DofNodes& map);

// Mean over the cells of (largest - smallest vertex id)/nverts: small where the
// file numbers its vertices with locality, about 0.6 where it does not
double node_span(const MeshData& m, int nv);

// Above this span the nodes are renumbered along the cells' Morton curve
constexpr double node_span_max = 0.25;

// The nodes along the same Morton curve as the cells: node_map[old] = new,
// vertices numbered first and edges after them
std::vector<std::uint32_t> morton_node_order(const MeshData& m, const std::vector<std::uint32_t>& edges,
                                             std::size_t nedges, int nv);

// The vertices themselves renumbered, node_map[old] = new: the coordinates
// moved and the topology rewritten
void renumber_vertices(MeshData& m, int nv, const std::vector<std::uint32_t>& node_map);

// The shortest edge of the mesh, a parallel reduction over the cells
double shortest_edge(const MeshData& m, int nv);

// Everything a load produces that is deterministic from the input files, so a
// later run reads it back instead of rebuilding it. The cell tree and the
// per-cell geometry are not here: both are rebuilt from the topology.
struct CacheTables {
  std::vector<std::uint32_t> topo;                 // ncells x nv
  std::vector<double> coords;                      // npoints x gdim
  std::vector<std::int32_t> facets;                // ncells x nv
  std::vector<std::uint32_t> u_nodes, p_nodes, phi_nodes;   // the per-cell node tables
  std::vector<double> u_values, p_values, phi_values;       // the first stamp, by node
  mesh_tables::DofNodes u_map, p_map, phi_map;              // stored dof -> the node slots it feeds
  std::size_t ncells = 0, nverts = 0;
  Uint gdim = 3, ncoeffs_u = 0, ncoeffs_p = 0, ncoeffs_phi = 0;
  std::size_t ncomp_u = 1;
  double hmin = 0.;
  Vector3d x_min = Vector3d::Zero(), x_max = Vector3d::Zero();
  std::string stamp;                               // the field file the values came from
};

// The identity of the inputs a cache was built from: path, size and mtime each
std::string cache_key(const std::vector<std::string>& inputs);

// True when the cache at path matches the key and was read
bool cache_read(const std::string& path, const std::string& key, CacheTables& c);

// Beside the input, if the directory takes it; a failure to write is not an error
void cache_write(const std::string& path, const std::string& key, const CacheTables& c);

}  // namespace simplex_load

#endif
