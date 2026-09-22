#include "simplex_load.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <numeric>
#include <omp.h>
#include <sys/stat.h>
#include <unistd.h>

#include <array>

#include "Error.hpp"
#include "h5direct.hpp"
#include "mesh_tables.hpp"
#include "morton.hpp"
#include "phase_timing.hpp"

namespace simplex_load {

namespace {

// The cell name dolfin writes into a signature
const char* cell_name(const int nv){ return nv == 4 ? "tetrahedron" : "triangle"; }

// The first unsigned integer at or after pos, or none
bool read_uint(const std::string& s, std::size_t pos, std::size_t& value){
  while (pos < s.size() && !std::isdigit(static_cast<unsigned char>(s[pos]))) ++pos;
  if (pos >= s.size()) return false;
  value = 0;
  while (pos < s.size() && std::isdigit(static_cast<unsigned char>(s[pos]))){
    value = 10*value + std::size_t(s[pos] - '0');
    ++pos;
  }
  return true;
}

// Counts a CellPos, a topology index and a 32-bit table can hold
void check_counts(const std::size_t ncells, const std::size_t nverts, const std::string& path){
  if (ncells > std::size_t(std::numeric_limits<std::int32_t>::max()))
    partrac::fail(path, ": ", ncells, " cells, more than a cell id can hold");
  if (nverts > std::size_t(std::numeric_limits<std::uint32_t>::max()))
    partrac::fail(path, ": ", nverts, " vertices, more than a vertex id can hold");
}

}  // namespace

void read_mesh_arrays(const std::string& path, const std::string& topology,
                      const std::string& geometry, const int nv, MeshData& m,
                      std::vector<std::uint32_t>& perm){
  const partrac::H5Id file = partrac::h5_open_read(path);

  const partrac::H5DatasetInfo topo_info = partrac::h5_dataset_info(file, topology);
  if (topo_info.rank() != 2 || topo_info.cols() != std::size_t(nv))
    partrac::fail(path, ": ", topology, " has ", topo_info.cols(), " vertices a cell, not ", nv);
  const partrac::H5DatasetInfo coord_info = partrac::h5_dataset_info(file, geometry);
  if (coord_info.rank() != 2)
    partrac::fail(path, ": ", geometry, " is not a table of points");

  m.ncells = topo_info.rows();
  m.nverts = coord_info.rows();
  m.gdim = Uint(coord_info.cols());
  if (m.gdim < Uint(nv) - 1)
    partrac::fail(path, ": ", geometry, " has ", m.gdim, " columns, too few for a ", nv,
                  "-vertex cell");
  m.gdim = Uint(nv) - 1;   // a 2D mesh stored with three columns keeps the first two
  check_counts(m.ncells, m.nverts, path);

  partrac::h5_read(file, topology, m.topo);
  partrac::phase(topology.c_str());
  partrac::h5_read(file, geometry, m.coords, m.gdim);
  partrac::phase(geometry.c_str());

  bool out_of_range = false;
#pragma omp parallel for schedule(static) reduction(||: out_of_range)
  for (std::size_t i = 0; i < m.topo.size(); ++i)
    if (m.topo[i] >= m.nverts) out_of_range = true;
  if (out_of_range)
    partrac::fail(path, ": ", topology, " names a vertex outside the ", m.nverts, " coordinates");

  // Bounds, then the cells along the Morton curve of their centroids
  for (Uint d = 0; d < 3; ++d){ m.x_min[d] = 0.; m.x_max[d] = 0.; }
  for (Uint d = 0; d < m.gdim; ++d){
    double lo = m.coords.empty() ? 0. : m.coords[d], hi = lo;
#pragma omp parallel for schedule(static) reduction(min: lo) reduction(max: hi)
    for (std::size_t i = 0; i < m.nverts; ++i){
      const double v = m.coords[i*m.gdim + d];
      lo = std::min(lo, v);
      hi = std::max(hi, v);
    }
    m.x_min[d] = lo;
    m.x_max[d] = hi;
  }

  std::vector<double> centre(3*m.ncells, 0.);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < m.ncells; ++i)
    for (int k = 0; k < nv; ++k)
      for (Uint d = 0; d < m.gdim; ++d)
        centre[3*i + d] += m.coords[std::size_t(m.topo[i*std::size_t(nv) + std::size_t(k)])*m.gdim + d]/double(nv);
  const partrac::MortonBox box(m.x_min, m.x_max, int(m.gdim));
  perm = partrac::morton_order(centre.data(), m.ncells, 3, box);
  centre.clear();
  centre.shrink_to_fit();

  std::vector<std::uint32_t> topo(m.topo.size());
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < m.ncells; ++i)
    for (int k = 0; k < nv; ++k)
      topo[i*std::size_t(nv) + std::size_t(k)] =
        m.topo[std::size_t(perm[i])*std::size_t(nv) + std::size_t(k)];
  m.topo.swap(topo);
  m.row_of.clear();
  partrac::phase("cells into Morton order");
}

void read_mesh(const std::string& path, const int nv, MeshData& m,
               std::vector<std::uint32_t>& perm){
  read_mesh_arrays(path, "mesh/topology", "mesh/coordinates", nv, m, perm);

  const partrac::H5Id file = partrac::h5_open_read(path);
  std::vector<std::uint64_t> gid;
  partrac::h5_read(file, "mesh/cell_indices", gid);
  if (gid.size() != m.ncells)
    partrac::fail(path, ": mesh/cell_indices has ", gid.size(), " entries for ", m.ncells, " cells");

  // Global cell id -> the cell here, the composition the field file needs
  std::uint64_t gid_max = 0;
  for (std::size_t i = 0; i < m.ncells; ++i) gid_max = std::max(gid_max, gid[i]);
  if (m.ncells > 0 && gid_max >= std::uint64_t(no_cell))
    partrac::fail(path, ": mesh/cell_indices holds ", gid_max, ", which a 32-bit id cannot hold");
  m.row_of.assign(m.ncells ? std::size_t(gid_max) + 1 : 0, no_cell);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < m.ncells; ++i)
    m.row_of[std::size_t(gid[std::size_t(perm[i])])] = std::uint32_t(i);
  partrac::phase("mesh/cell_indices");
}

Element read_element(const std::string& path, const std::string& field, const int nv){
  const partrac::H5Id file = partrac::h5_open_read(path);
  const std::string sig = partrac::h5_read_string_attribute(file, field, "signature");
  partrac::phase("signature");
  Element el;
  const bool lagrange = sig.find("Lagrange") != std::string::npos || sig.find("'CG'") != std::string::npos;
  const std::size_t at_cell = sig.find(cell_name(nv));
  if (!lagrange || at_cell == std::string::npos)
    partrac::fail(path, ": '", field, "' is ", sig, ", not a Lagrange element on a ", cell_name(nv));
  std::size_t degree = 0;
  if (!read_uint(sig, at_cell + std::strlen(cell_name(nv)), degree))
    partrac::fail(path, ": '", field, "' is ", sig, ", which names no degree");
  if (degree != 1 && degree != 2)
    partrac::fail(path, ": '", field, "' is of degree ", degree, "; this loader reads P1 and P2");
  el.degree = int(degree);
  el.ncomp = 1;
  if (sig.compare(0, 13, "VectorElement") == 0){
    const std::size_t at_dim = sig.find("dim=");
    std::size_t ncomp = 0;
    if (at_dim != std::string::npos ? !read_uint(sig, at_dim, ncomp)
                                    : !read_uint(sig, sig.rfind(',') , ncomp))
      partrac::fail(path, ": '", field, "' is ", sig, ", which names no component count");
    if (ncomp != std::size_t(nv) - 1)
      partrac::fail(path, ": '", field, "' has ", ncomp, " components, not ", nv - 1);
    el.ncomp = ncomp;
  }
  else if (sig.compare(0, 13, "FiniteElement") != 0){
    partrac::fail(path, ": '", field, "' is ", sig, ", neither a scalar nor a vector element");
  }
  return el;
}

void read_field(const std::string& path, const std::string& field, const MeshData& m,
                const Element& el, const int nv,
                std::vector<std::uint32_t>& rows, std::vector<double>& vec){
  const partrac::H5Id file = partrac::h5_open_read(path);
  const std::size_t per_cell = nodes_per_cell(nv, el.degree)*el.ncomp;

  // The row offsets must be a constant stride, or the table cannot be reshaped
  std::vector<std::uint64_t> offsets;
  partrac::h5_read(file, field + "/x_cell_dofs", offsets);
  partrac::phase("x_cell_dofs");
  if (offsets.size() != m.ncells + 1)
    partrac::fail(path, ": '", field, "/x_cell_dofs' has ", offsets.size(), " entries for ",
                  m.ncells, " cells");
  std::size_t ragged = m.ncells;
#pragma omp parallel for schedule(static) reduction(min: ragged)
  for (std::size_t i = 0; i < m.ncells; ++i)
    if (offsets[i+1] - offsets[i] != std::uint64_t(per_cell)) ragged = std::min(ragged, i);
  if (ragged < m.ncells)
    partrac::fail(path, ": '", field, "/x_cell_dofs' is ragged; row ", ragged, " holds ",
                  offsets[ragged+1] - offsets[ragged], " dofs, not ", per_cell);
  offsets.clear();
  offsets.shrink_to_fit();

  const std::size_t n = m.ncells*per_cell;
  const partrac::H5DatasetInfo dof_info = partrac::h5_dataset_info(file, field + "/cell_dofs");
  if (dof_info.count() != n)
    partrac::fail(path, ": '", field, "/cell_dofs' holds ", dof_info.count(), " dofs, not ", n);
  std::vector<std::uint32_t> flat(n);
  if (dof_info.type_class == H5T_INTEGER && dof_info.type_size == 4 && dof_info.type_signed){
    // Same width on disk: read into the destination and convert it where it lies
    partrac::h5_fill<std::int32_t>(file, field + "/cell_dofs", dof_info,
                                   reinterpret_cast<std::int32_t*>(flat.data()), n);
    bool negative = false;
#pragma omp parallel for schedule(static) reduction(||: negative)
    for (std::size_t i = 0; i < n; ++i){
      std::int32_t v;
      std::memcpy(&v, &flat[i], sizeof v);
      if (v < 0) negative = true;
      else flat[i] = std::uint32_t(v);
    }
    if (negative)
      partrac::fail(path, ": '", field, "/cell_dofs' holds a negative dof index");
  }
  else {
    partrac::h5_read(file, field + "/cell_dofs", flat);
  }
  partrac::phase("cell_dofs");

  std::vector<std::uint64_t> cells;
  partrac::h5_read(file, field + "/cells", cells);
  partrac::phase("cells");
  if (cells.size() != m.ncells)
    partrac::fail(path, ": '", field, "/cells' names ", cells.size(), " cells, the mesh has ", m.ncells);

  partrac::h5_read(file, field + "/vector_0", vec);
  partrac::phase("vector_0");
  if (vec.size() > std::size_t(std::numeric_limits<std::uint32_t>::max()))
    partrac::fail(path, ": '", field, "/vector_0' holds ", vec.size(),
                  " values, more than a 32-bit dof index can hold");

  // Both files label their rows by global cell id, and neither is in that
  // order: the field's row i is the mesh's row_of[cells[i]]
  std::vector<std::uint32_t> dest(m.ncells);
  bool missing = false;
  bool identity = true;
#pragma omp parallel for schedule(static) reduction(||: missing) reduction(&&: identity)
  for (std::size_t i = 0; i < m.ncells; ++i){
    const std::uint64_t g = cells[i];
    const std::uint32_t r = g < m.row_of.size() ? m.row_of[std::size_t(g)] : no_cell;
    if (r == no_cell) missing = true;
    dest[i] = r;
    if (r != std::uint32_t(i)) identity = false;
  }
  if (missing){
    const std::size_t i = std::find(dest.begin(), dest.end(), no_cell) - dest.begin();
    partrac::fail(path, ": '", field, "/cells' names cell ", cells[i],
                  ", which the mesh file does not");
  }
  if (identity){
    // The writer's order is the mesh's: no row moves
    rows.swap(flat);
  }
  else {
    rows.assign(n, 0);
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < m.ncells; ++i)
      std::copy(flat.data() + i*per_cell, flat.data() + (i+1)*per_cell,
                rows.data() + std::size_t(dest[i])*per_cell);
  }
  partrac::phase("compose through cell_indices");
}

void read_vector(const std::string& path, const std::string& field,
                 const mesh_tables::DofNodes& map, std::vector<double>& values){
  const partrac::H5Id file = partrac::h5_open_read(path);
  std::vector<double> vec;
  partrac::h5_read(file, field + "/vector_0", vec);
  if (map.start.empty() || vec.size() != map.start.size() - 1)
    partrac::fail(path, ": '", field, "/vector_0' holds ", vec.size(),
                  " values, the first stamp ", map.start.empty() ? 0 : map.start.size() - 1);
  // A mapping built for another numbering would write outside the values
  std::uint32_t top = 0;
#pragma omp parallel for schedule(static) reduction(max: top)
  for (std::size_t q = 0; q < map.slot.size(); ++q)
    top = std::max(top, map.slot[q]);
  if (!map.slot.empty() && std::size_t(top) >= values.size())
    partrac::fail(path, ": '", field, "' is read into node slot ", top,
                  ", outside the ", values.size(), " the field holds");
  // Each dof owns its slots, so the threads never write the same one
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < vec.size(); ++i)
    for (std::size_t q = map.start[i]; q < map.start[i + 1]; ++q)
      values[map.slot[q]] = vec[i];
}

std::vector<double> node_points(const MeshData& m, const std::vector<std::uint32_t>& edges,
                                const std::size_t nedges, const int nv, const std::size_t stride){
  const int ne = nv*(nv-1)/2;
  const std::size_t n_total = m.nverts + nedges;
  std::vector<double> x(stride*n_total, 0.);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < m.nverts; ++i)
    for (Uint d = 0; d < m.gdim; ++d) x[stride*i + d] = m.coords[i*m.gdim + d];
  // Every cell sharing an edge writes the same midpoint, bit for bit, since the
  // endpoints are taken in vertex-id order; dolfin's local order, as mesh_tables has it
  const std::size_t n_loc = std::size_t(ne);
  std::vector<std::array<int, 2>> loc(n_loc);
  {
    int i = 0;
    for (int a = 0; a < nv; ++a)
      for (int b = a + 1; b < nv; ++b){ loc[std::size_t(ne - 1 - i)] = {a, b}; ++i; }
  }
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < (edges.empty() ? 0 : m.ncells); ++i){
    const std::uint32_t* row = m.topo.data() + i*std::size_t(nv);
    for (int e = 0; e < ne; ++e){
      const int a = loc[std::size_t(e)][0];
      const int b = loc[std::size_t(e)][1];
      const std::uint32_t lo = std::min(row[a], row[b]), hi = std::max(row[a], row[b]);
      double* p = x.data() + stride*(m.nverts + std::size_t(edges[i*std::size_t(ne) + std::size_t(e)]));
      for (Uint d = 0; d < m.gdim; ++d)
        p[d] = 0.5*(m.coords[std::size_t(lo)*m.gdim + d] + m.coords[std::size_t(hi)*m.gdim + d]);
    }
  }
  return x;
}

void NodePairs::build(const MeshData& m, const std::vector<std::uint32_t>& edges,
                      const std::size_t nedges, const int nv, const std::vector<bool>& periodic,
                      const Vector3d& x_min, const Vector3d& x_max, const double tol){
  bool any = false;
  for (Uint d = 0; d < m.gdim; ++d) any = any || periodic[d];
  if (!any) return;
  gdim = m.gdim;
  for (Uint d = 0; d < gdim; ++d)
    if (periodic[d]) period[d] = x_max[d] - x_min[d];
  if (!edges.empty()) points = node_points(m, edges, nedges, nv, gdim);
  x = edges.empty() ? &m.coords : &points;
  master = mesh_tables::match_periodic_vertices(*x, m.nverts + nedges, gdim, periodic,
                                                x_min, x_max, tol);
  partrac::phase("periodic nodes");
}

void NodePairs::check(const std::vector<double>& values, const std::size_t ncomp,
                      const std::string& path, const std::string& field) const {
  if (master.empty()) return;
  const std::size_t n = values.size()/ncomp;
  double scale = 0.;
  for (const double v : values) scale = std::max(scale, std::abs(v));
  const double tol = 1e-8*scale;
  // A disagreeing pair is a periodic constraint the solver missed; the master serves both
  double worst = 0.;
  std::size_t worst_node = 0, pairs = 0, off = 0;
  for (std::size_t i = 0; i < n; ++i){
    const std::size_t mst = master[i];
    if (mst == i) continue;
    ++pairs;
    double d_i = 0.;
    for (std::size_t c = 0; c < ncomp; ++c)
      d_i = std::max(d_i, std::abs(values[i*ncomp + c] - values[mst*ncomp + c]));
    if (d_i > tol) ++off;
    if (d_i > worst){ worst = d_i; worst_node = i; }
  }
  if (off == 0) return;
  // The axis the worst pair crosses, the first of them at a corner
  const std::size_t mst = master[worst_node];
  char axis = '?';
  for (Uint d = 0; d < gdim; ++d)
    if (period[d] > 0. && axis == '?'
        && std::abs((*x)[worst_node*gdim + d] - (*x)[mst*gdim + d]) > 0.5*period[d])
      axis = char('x' + d);
  std::cout << path << ": '" << field << "': " << off << " of " << pairs
            << " periodic node pairs disagree (worst " << worst << " along " << axis
            << ", tolerance " << tol << "); masters used" << std::endl;
}

std::vector<std::uint32_t> NodePairs::masters(const std::vector<std::uint32_t>& node_map,
                                              const std::size_t n) const {
  if (master.empty()) return node_map;
  std::vector<std::uint32_t> out(master.begin(), master.begin() + std::ptrdiff_t(n));
  if (!node_map.empty())
    for (std::size_t i = 0; i < n; ++i) out[i] = node_map[out[i]];
  return out;
}

double node_span(const MeshData& m, const int nv){
  if (m.ncells == 0 || m.nverts == 0) return 0.;
  double sum = 0.;
#pragma omp parallel for schedule(static) reduction(+: sum)
  for (std::size_t i = 0; i < m.ncells; ++i){
    const std::uint32_t* row = m.topo.data() + i*std::size_t(nv);
    const auto lohi = std::minmax_element(row, row + nv);
    sum += double(*lohi.second - *lohi.first);
  }
  return sum / (double(m.ncells) * double(m.nverts));
}

std::vector<std::uint32_t> morton_node_order(const MeshData& m, const std::vector<std::uint32_t>& edges,
                                             const std::size_t nedges, const int nv){
  const std::size_t n_total = m.nverts + nedges;
  const std::vector<double> x = node_points(m, edges, nedges, nv, 3);
  const partrac::MortonBox box(m.x_min, m.x_max, int(m.gdim));
  const std::vector<std::uint32_t> order = partrac::morton_order(x.data(), n_total, 3, box);
  std::vector<std::uint32_t> node_map(n_total);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < n_total; ++i) node_map[order[i]] = std::uint32_t(i);
  return node_map;
}

double shortest_edge(const MeshData& m, const int nv){
  double h = std::numeric_limits<double>::max();
#pragma omp parallel for schedule(static) reduction(min: h)
  for (std::size_t i = 0; i < m.ncells; ++i){
    const std::uint32_t* row = m.topo.data() + i*std::size_t(nv);
    for (int a = 0; a < nv; ++a)
      for (int b = a + 1; b < nv; ++b){
        double d2 = 0.;
        for (Uint k = 0; k < m.gdim; ++k){
          const double dd = m.coords[std::size_t(row[a])*m.gdim + k] - m.coords[std::size_t(row[b])*m.gdim + k];
          d2 += dd*dd;
        }
        h = std::min(h, d2);
      }
  }
  return m.ncells ? std::sqrt(h) : 0.;
}

int declared_degree(const std::string& space, const char* what){
  if (space == "P1") return 1;
  if (space == "P2") return 2;
  partrac::fail("unrecognized ", what, " element: ", space);
  return 0;
}

void request_from_params(Request& r, const partrac::Params& prm, const std::string& folder){
  r.mesh_file = folder + "/" + prm.get<std::string>("mesh");
  r.u_field = prm.get<std::string>("velocity_field");
  r.p_field = prm.get<std::string>("pressure_field");
  // Only the schemas that read a phase field declare its dataset
  if (r.include_phi) r.phi_field = prm.get<std::string>("phase_field");
  // The file's own signature decides the element; a name here it does not know
  // is still refused
  r.want_u = declared_degree(prm.get<std::string>("velocity_space"), "velocity");
  r.want_p = r.include_pressure ? declared_degree(prm.get<std::string>("pressure_space"), "pressure") : 0;
}

void build_tables(const Request& r, Tables& t){
  const int nv = r.nv;
  const int ne = nv*(nv-1)/2;
  t.nv = nv;
  MeshData& m = t.mesh;
  {
    std::vector<std::uint32_t> cell_perm;
    read_mesh(r.mesh_file, nv, m, cell_perm);
  }

  // The element comes from the file; the parameter file must not claim another
  t.el_u = read_element(r.field_file, r.u_field, nv);
  if (r.include_pressure)
    t.el_p = read_element(r.field_file, r.p_field, nv);
  if (r.want_u != t.el_u.degree)
    partrac::fail(r.infilename, ": velocity_space is P", r.want_u, ", but '", r.u_field,
                  "' in ", r.field_file, " is of degree ", t.el_u.degree);
  if (r.include_pressure && r.want_p != t.el_p.degree)
    partrac::fail(r.infilename, ": pressure_space is P", r.want_p, ", but '", r.p_field,
                  "' in ", r.field_file, " is of degree ", t.el_p.degree);
  if (t.el_u.ncomp != std::size_t(nv - 1))
    partrac::fail(r.field_file, ": '", r.u_field, "' has ", t.el_u.ncomp, " components, not ", nv - 1);
  if (r.include_pressure && t.el_p.ncomp != 1)
    partrac::fail(r.field_file, ": '", r.p_field, "' has ", t.el_p.ncomp, " components, not one");
  if (r.include_phi){
    t.el_phi = read_element(r.field_file, r.phi_field, nv);
    if (t.el_phi.ncomp != 1)
      partrac::fail(r.field_file, ": '", r.phi_field, "' has ", t.el_phi.ncomp, " components, not one");
  }
  t.ncoeffs_u = Uint(nodes_per_cell(nv, t.el_u.degree));
  t.ncoeffs_p = r.include_pressure ? Uint(nodes_per_cell(nv, t.el_p.degree)) : 0;
  t.ncoeffs_phi = r.include_phi ? Uint(nodes_per_cell(nv, t.el_phi.degree)) : 0;
  check_dofs_fit(t.ncoeffs_u, std::max(t.ncoeffs_p, t.ncoeffs_phi), r.n_dofs_max, r.what);

  // The edges, when any field carries midside nodes
  const bool quadratic = t.el_u.degree == 2 || (r.include_pressure && t.el_p.degree == 2)
                      || (r.include_phi && t.el_phi.degree == 2);
  if (quadratic)
    t.nedges = nv == 4 ? mesh_tables::build_edge_table<4>(m.topo, m.ncells, t.edges)
                       : mesh_tables::build_edge_table<3>(m.topo, m.ncells, t.edges);
  partrac::phase("edge table");

  // The nodes along the cells' Morton curve, where the file numbers them without locality
  const double span = node_span(m, nv);
  std::cout << "Mean vertex-id span of a cell: " << span << " of the vertices" << std::endl;
  if (span > node_span_max){
    std::cout << "Node order: renumbering nodes along the cells' curve" << std::endl;
    if (quadratic)
      t.map_quad = morton_node_order(m, t.edges, t.nedges, nv);
    if (t.el_u.degree == 1 || (r.include_pressure && t.el_p.degree == 1)
        || (r.include_phi && t.el_phi.degree == 1))
      t.map_lin = morton_node_order(m, std::vector<std::uint32_t>(), 0, nv);
    partrac::phase("node order");
  }

  // A node and its periodic images are one node, the master's dof serving all
  t.np.build(m, t.edges, t.nedges, nv, r.periodic, m.x_min, m.x_max, r.periodic_tol);

  const std::size_t n_u = m.nverts + (t.el_u.degree == 2 ? t.nedges : 0);
  const std::vector<std::uint32_t> nodes_u = t.np.masters(t.node_order(t.el_u), n_u);
  t.u_dofs.fill(m.topo, t.edges, m.ncells, nv, ne, t.el_u.degree == 2, m.nverts,
                nodes_u.empty() ? nullptr : nodes_u.data());
  t.u_dofs.check_stride(t.ncoeffs_u, r.what);
  if (r.include_pressure){
    const std::size_t n_p = m.nverts + (t.el_p.degree == 2 ? t.nedges : 0);
    const std::vector<std::uint32_t> nodes_p = t.np.masters(t.node_order(t.el_p), n_p);
    t.p_dofs.fill(m.topo, t.edges, m.ncells, nv, ne, t.el_p.degree == 2, m.nverts,
                  nodes_p.empty() ? nullptr : nodes_p.data());
    t.p_dofs.check_stride(t.ncoeffs_p, r.what);
  }
  if (r.include_phi){
    const std::size_t n_phi = m.nverts + (t.el_phi.degree == 2 ? t.nedges : 0);
    const std::vector<std::uint32_t> nodes_phi = t.np.masters(t.node_order(t.el_phi), n_phi);
    t.phi_dofs.fill(m.topo, t.edges, m.ncells, nv, ne, t.el_phi.degree == 2, m.nverts,
                    nodes_phi.empty() ? nullptr : nodes_phi.data());
    t.phi_dofs.check_stride(t.ncoeffs_phi, r.what);
  }
  partrac::phase("cell dofs");

  if (nv == 4)
    mesh_tables::build_facet_neighbours<4>(m.topo, m.ncells, m.coords, m.gdim, r.periodic,
                                           m.x_min, m.x_max, r.periodic_tol, t.facets);
  else
    mesh_tables::build_facet_neighbours<3>(m.topo, m.ncells, m.coords, m.gdim, r.periodic,
                                           m.x_min, m.x_max, r.periodic_tol, t.facets);
  partrac::phase("facet table");

  t.hmin = shortest_edge(m, nv);
}

void read_field_by_node(const std::string& path, const std::string& field, const Tables& t,
                        const Element& el, const std::vector<std::uint32_t>& node_map,
                        std::vector<double>& values, mesh_tables::DofNodes& map){
  const MeshData& m = t.mesh;
  const bool quadratic = el.degree == 2;
  const std::vector<std::uint32_t> no_edges;
  const std::vector<std::uint32_t>& edges = quadratic ? t.edges : no_edges;
  const std::size_t nedges = quadratic ? t.nedges : 0;
  std::vector<std::uint32_t> rows;
  std::vector<double> vec;
  read_field(path, field, m, el, t.nv, rows, vec);
  for (std::size_t i = 0; i < vec.size(); ++i)
    if (!std::isfinite(vec[i])) partrac::fail(path, ": ", field, " holds a non-finite value at dof ", i);
  if (t.nv == 4)
    mesh_tables::scatter_dofs_to_nodes<4>(m.topo, edges, m.ncells, m.nverts, nedges, el.ncomp,
                                          rows, vec, values, map);
  else
    mesh_tables::scatter_dofs_to_nodes<3>(m.topo, edges, m.ncells, m.nverts, nedges, el.ncomp,
                                          rows, vec, values, map);
  rows.clear();
  rows.shrink_to_fit();
  partrac::phase("scatter");
  t.np.check(values, el.ncomp, path, field);
  if (node_map.empty()) return;
  const std::size_t ncomp = el.ncomp;
  const std::size_t n_total = values.size()/ncomp;
  std::vector<double> moved(values.size());
#pragma omp parallel for schedule(static)
  for (std::size_t n = 0; n < n_total; ++n)
    for (std::size_t c = 0; c < ncomp; ++c)
      moved[std::size_t(node_map[n])*ncomp + c] = values[n*ncomp + c];
  values.swap(moved);
#pragma omp parallel for schedule(static)
  for (std::size_t q = 0; q < map.slot.size(); ++q)
    map.slot[q] = std::uint32_t(std::size_t(node_map[map.slot[q]/ncomp])*ncomp + map.slot[q] % ncomp);
  partrac::phase("node renumbering");
}

namespace {

// The cache's own format; a change here invalidates every cache written before
constexpr int cache_version = 3;

// Scalars, in the order cache_write puts them
enum Scalar { S_NCELLS, S_NVERTS, S_GDIM, S_NCOEFFS_U, S_NCOEFFS_P, S_NCOMP_U, S_HMIN,
              S_XMIN, S_XMAX = S_XMIN + 3, S_NCOEFFS_PHI = S_XMAX + 3, S_COUNT };

template<typename T>
void cache_put(const hid_t file, const char* name, const std::vector<T>& v){
  if (v.empty()) return;
  const hsize_t n = v.size();
  const partrac::H5Id space(H5Screate_simple(1, &n, nullptr), H5Sclose);
  const hid_t type = partrac::h5_native_type<T>();
  // the default creation list is contiguous and unfiltered, which is what h5direct reads
  const partrac::H5Id dset(H5Dcreate2(file, name, type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                           H5Dclose);
  if (!dset.valid() || H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) < 0)
    partrac::fail("mesh cache: cannot write '", name, "'");
}

// A cache whose arrays do not fit its scalars is not read; the load rebuilds
bool cache_bad(const std::string& path, const std::string& what){
  std::cout << "Mesh cache: " << path << " does not fit its counts (" << what
            << "), rebuilding" << std::endl;
  return false;
}

template<typename T>
bool cache_get(const hid_t file, const std::string& name, std::vector<T>& v, const std::size_t n){
  v.clear();
  if (n == 0) return true;
  partrac::h5_read(file, name, v);
  return v.size() == n;
}

// Every node a cell table names and every slot the mapping feeds is a node of
// the values: a truncated file would otherwise be written past its end
bool values_fit(const std::vector<double>& values, const std::vector<std::uint32_t>& nodes,
                const mesh_tables::DofNodes& map, const std::size_t ncomp){
  if (ncomp == 0 || values.empty() || values.size() % ncomp != 0)
    return false;
  const std::size_t n_total = values.size()/ncomp;
  for (const std::uint32_t n : nodes)
    if (std::size_t(n) >= n_total) return false;
  for (const std::uint32_t q : map.slot)
    if (std::size_t(q) >= values.size()) return false;
  return true;
}

// The dof -> nodes mapping as its two arrays; an empty one writes nothing
void cache_put_map(const hid_t file, const std::string& field, const mesh_tables::DofNodes& m){
  cache_put(file, (field + "_start").c_str(), m.start);
  cache_put(file, (field + "_slot").c_str(), m.slot);
}

bool cache_get_map(const hid_t file, const std::string& field, mesh_tables::DofNodes& m){
  partrac::h5_read(file, field + "_start", m.start);
  partrac::h5_read(file, field + "_slot", m.slot);
  return !m.start.empty() && m.start.back() == m.slot.size();
}

void cache_put_string(const hid_t file, const char* name, const std::string& value){
  const partrac::H5Id type(H5Tcopy(H5T_C_S1), H5Tclose);
  H5Tset_size(type, value.size() + 1);
  const partrac::H5Id space(H5Screate(H5S_SCALAR), H5Sclose);
  const partrac::H5Id attr(H5Acreate2(file, name, type, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
  if (!attr.valid() || H5Awrite(attr, type, value.c_str()) < 0)
    partrac::fail("mesh cache: cannot write the attribute '", name, "'");
}

}  // namespace

std::string cache_key(const std::vector<std::string>& inputs){
  std::string key = "partrac-simplex-cache v" + std::to_string(cache_version);
  for (const std::string& p : inputs){
    struct stat st;
    if (::stat(p.c_str(), &st) != 0)
      return std::string();
    key += "|" + p + ":" + std::to_string(std::uint64_t(st.st_size)) + ":"
         + std::to_string(std::int64_t(st.st_mtime)) + "." + std::to_string(long(st.st_mtim.tv_nsec));
  }
  return key;
}

bool cache_read(const std::string& path, const std::string& key, CacheTables& c){
  if (key.empty() || ::access(path.c_str(), R_OK) != 0)
    return false;
  const hid_t raw = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (raw < 0)
    return false;
  const partrac::H5Id file(raw, H5Fclose);
  if (H5Aexists(file, "key") <= 0)
    return false;
  if (partrac::h5_read_string_attribute(file, "", "key") != key)
    return false;
  c.stamp = partrac::h5_read_string_attribute(file, "", "stamp");

  std::vector<double> sc;
  partrac::h5_read(file, "scalars", sc);
  if (sc.size() != std::size_t(S_COUNT))
    partrac::fail("mesh cache: ", path, " holds ", sc.size(), " scalars, not ", int(S_COUNT));
  c.ncells = std::size_t(sc[S_NCELLS]);
  c.nverts = std::size_t(sc[S_NVERTS]);
  c.gdim = Uint(sc[S_GDIM]);
  c.ncoeffs_u = Uint(sc[S_NCOEFFS_U]);
  c.ncoeffs_p = Uint(sc[S_NCOEFFS_P]);
  c.ncoeffs_phi = Uint(sc[S_NCOEFFS_PHI]);
  c.ncomp_u = std::size_t(sc[S_NCOMP_U]);
  c.hmin = sc[S_HMIN];
  for (int d = 0; d < 3; ++d){ c.x_min[d] = sc[S_XMIN + d]; c.x_max[d] = sc[S_XMAX + d]; }

  const std::size_t nv = c.gdim + 1;
  if (!cache_get(file, "topology", c.topo, c.ncells*nv))
    return cache_bad(path, "the topology");
  if (!cache_get(file, "coordinates", c.coords, c.nverts*c.gdim))
    return cache_bad(path, "the coordinates");
  if (!cache_get(file, "facets", c.facets, c.ncells*nv))
    return cache_bad(path, "the facet table");
  if (!cache_get(file, "u_nodes", c.u_nodes, c.ncells*std::size_t(c.ncoeffs_u)))
    return cache_bad(path, "the velocity node table");
  partrac::h5_read(file, "u_values", c.u_values);
  if (!cache_get_map(file, "u", c.u_map) || !values_fit(c.u_values, c.u_nodes, c.u_map, c.ncomp_u))
    return cache_bad(path, "the velocity values");
  if (c.ncoeffs_p > 0){
    if (!cache_get(file, "p_nodes", c.p_nodes, c.ncells*std::size_t(c.ncoeffs_p)))
      return cache_bad(path, "the pressure node table");
    partrac::h5_read(file, "p_values", c.p_values);
    if (!cache_get_map(file, "p", c.p_map) || !values_fit(c.p_values, c.p_nodes, c.p_map, 1))
      return cache_bad(path, "the pressure values");
  }
  if (c.ncoeffs_phi > 0){
    if (!cache_get(file, "phi_nodes", c.phi_nodes, c.ncells*std::size_t(c.ncoeffs_phi)))
      return cache_bad(path, "the phase field node table");
    partrac::h5_read(file, "phi_values", c.phi_values);
    if (!cache_get_map(file, "phi", c.phi_map) || !values_fit(c.phi_values, c.phi_nodes, c.phi_map, 1))
      return cache_bad(path, "the phase field values");
  }
  return true;
}

void cache_write(const std::string& path, const std::string& key, const CacheTables& c){
  if (key.empty())
    return;
  const std::size_t slash = path.find_last_of('/');
  const std::string dir = slash == std::string::npos ? std::string(".") : path.substr(0, slash);
  if (::access(dir.c_str(), W_OK) != 0){
    std::cout << "Mesh cache: " << dir << " is not writable, not caching" << std::endl;
    return;
  }
  const hid_t raw = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (raw < 0){
    std::cout << "Mesh cache: cannot create " << path << ", not caching" << std::endl;
    return;
  }
  const partrac::H5Id file(raw, H5Fclose);
  std::vector<double> sc(std::size_t(S_COUNT), 0.);
  sc[S_NCELLS] = double(c.ncells);
  sc[S_NVERTS] = double(c.nverts);
  sc[S_GDIM] = double(c.gdim);
  sc[S_NCOEFFS_U] = double(c.ncoeffs_u);
  sc[S_NCOEFFS_P] = double(c.ncoeffs_p);
  sc[S_NCOEFFS_PHI] = double(c.ncoeffs_phi);
  sc[S_NCOMP_U] = double(c.ncomp_u);
  sc[S_HMIN] = c.hmin;
  for (int d = 0; d < 3; ++d){ sc[S_XMIN + d] = c.x_min[d]; sc[S_XMAX + d] = c.x_max[d]; }
  cache_put(file, "scalars", sc);
  cache_put(file, "topology", c.topo);
  cache_put(file, "coordinates", c.coords);
  cache_put(file, "facets", c.facets);
  cache_put(file, "u_nodes", c.u_nodes);
  cache_put(file, "u_values", c.u_values);
  cache_put_map(file, "u", c.u_map);
  cache_put(file, "p_nodes", c.p_nodes);
  cache_put(file, "p_values", c.p_values);
  cache_put_map(file, "p", c.p_map);
  cache_put(file, "phi_nodes", c.phi_nodes);
  cache_put(file, "phi_values", c.phi_values);
  cache_put_map(file, "phi", c.phi_map);
  cache_put_string(file, "key", key);
  cache_put_string(file, "stamp", c.stamp);
}

}  // namespace simplex_load
