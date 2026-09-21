#include "SimplexInterpol.hpp"
#include "loader_params.hpp"
#include "p12_eval.hpp"
#include "mesh_tables.hpp"
#include "simplex_load.hpp"
#include "phase_timing.hpp"
#include <array>
#include <cassert>
#include <cstdio>
#include <iostream>

namespace {

// One field's values by node, and the mapping a later stamp is read through.
// node_map, when given, is this field's node order; the stored dof table is
// dropped as soon as the values are out of it.
template<int NV>
void load_field(const std::string& path, const std::string& field,
                const simplex_load::MeshData& m, const simplex_load::Element& el,
                const std::vector<std::uint32_t>& edges, const std::size_t nedges,
                const std::vector<std::uint32_t>& node_map, const simplex_load::NodePairs& np,
                std::vector<double>& values, mesh_tables::DofNodes& map){
  const bool quadratic = el.degree == 2;
  const std::vector<std::uint32_t> no_edges;
  std::vector<std::uint32_t> rows;
  std::vector<double> vec;
  simplex_load::read_field(path, field, m, el, NV, rows, vec);
  mesh_tables::scatter_dofs_to_nodes<NV>(m.topo, quadratic ? edges : no_edges, m.ncells,
                                         m.nverts, quadratic ? nedges : 0, el.ncomp,
                                         rows, vec, values, map);
  rows.clear();
  rows.shrink_to_fit();
  partrac::phase("scatter");
  np.check(values, el.ncomp, path, field);
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

// The element a parameter file declares; the file's own signature decides the
// element, but a name this loader does not know is still a mistake
int declared_degree(const std::string& space, const char* what){
  if (space == "P1") return 1;
  if (space == "P2") return 2;
  partrac::fail("unrecognized ", what, " element: ", space);
  return 0;
}

}  // namespace

template<typename Cell>
SimplexInterpol<Cell>::SimplexInterpol(const std::string& infilename)
  : MeshCore<Cell>(infilename)
{
  constexpr int nv = Cell::n_verts;
  constexpr int ne = nv*(nv-1)/2;
  partrac::phase_begin("load");
  dolfin_params = partrac::parse_file_or_exit(dolfin_h5_schema(mode), infilename);

  const std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  ts.initialize(get_folder() + "/" + dolfin_params.template get<std::string>("timestamps"));

  read_mesh_params();
  u_field_ = dolfin_params.template get<std::string>("velocity_field");
  p_field_ = dolfin_params.template get<std::string>("pressure_field");
  partrac::phase("params");

  const std::string mesh_file = get_folder() + "/" + dolfin_params.template get<std::string>("mesh");
  // The first stamp, by the path update builds for every other one
  const std::string first = get_folder() + "/" + ts.get(ts.get_t_min()).prev.filename;
  // The file's signature decides the element; a name here it does not know is
  // still refused, whether the tables are rebuilt or read from the cache
  const int want_u = declared_degree(dolfin_params.template get<std::string>("velocity_space"), "velocity");
  const int want_p = include_pressure
    ? declared_degree(dolfin_params.template get<std::string>("pressure_space"), "pressure") : 0;

  // The native cache, opt-in: everything below is deterministic from the two
  // input files, so a later run reads the tables back instead of rebuilding
  const bool use_cache = dolfin_params.template get<bool>("mesh_cache");
  const std::string cache_file = mesh_file.substr(0, mesh_file.find_last_of('.'))
                               + "_partrac_" + mode + ".h5";
  std::string key;
  if (use_cache){
    key = simplex_load::cache_key({mesh_file, first});
    if (!key.empty()){
      // The facet table and the node tables are built from the periodicity too
      char tol[32];
      std::snprintf(tol, sizeof tol, "%.17g", periodic_tol);
      key += "|u=" + u_field_ + "|p=" + (include_pressure ? p_field_ : std::string())
           + "|periodic=" + (periodic[0] ? "1" : "0") + (periodic[1] ? "1" : "0")
                          + (periodic[2] ? "1" : "0")
           + "|periodic_tol=" + tol;
    }
  }
  simplex_load::CacheTables cache;
  const bool cached = use_cache && simplex_load::cache_read(cache_file, key, cache)
                   && (!include_pressure || cache.ncoeffs_p > 0);
  if (cached){
    std::cout << "Mesh cache: read from " << cache_file << std::endl;
    dim = cache.gdim;
    x_min = cache.x_min;
    x_max = cache.x_max;
    hmin_ = cache.hmin;
    ncoeffs_u = cache.ncoeffs_u;
    ncoeffs_p = include_pressure ? cache.ncoeffs_p : 0;
    check_dofs_fit(ncoeffs_u, ncoeffs_p, Cell::n_dofs_max, "SimplexInterpol");
    topo_ = std::move(cache.topo);
    coords_ = std::move(cache.coords);
    facet_neigh_ = std::move(cache.facets);
    u_dofs_.adopt(std::move(cache.u_nodes), ncoeffs_u);
    u_a_ = std::move(cache.u_values);
    u_map_ = std::move(cache.u_map);
    if (include_pressure){
      p_dofs_.adopt(std::move(cache.p_nodes), ncoeffs_p);
      p_a_ = std::move(cache.p_values);
      p_map_ = std::move(cache.p_map);
    }
    file_a_ = cache.stamp;
    ncells_ = topo_.size()/nv;
    nverts_ = coords_.size()/dim;
    set_period();
    partrac::phase("mesh cache");
  }
  else {
    simplex_load::MeshData m;
    std::vector<std::uint32_t> cell_perm;
    simplex_load::read_mesh(mesh_file, nv, m, cell_perm);
    cell_perm.clear();
    cell_perm.shrink_to_fit();
    dim = m.gdim;
    x_min = m.x_min;
    x_max = m.x_max;
    set_period();

    // The element comes from the file; the parameter file must not claim another
    const simplex_load::Element el_u = simplex_load::read_element(first, u_field_, nv);
    simplex_load::Element el_p;
    if (include_pressure)
      el_p = simplex_load::read_element(first, p_field_, nv);
    if (want_u != el_u.degree)
      partrac::fail(infilename, ": velocity_space is P", want_u, ", but '", u_field_,
                    "' in ", first, " is of degree ", el_u.degree);
    if (include_pressure && want_p != el_p.degree)
      partrac::fail(infilename, ": pressure_space is P", want_p, ", but '", p_field_,
                    "' in ", first, " is of degree ", el_p.degree);
    if (el_u.ncomp != std::size_t(D))
      partrac::fail(first, ": '", u_field_, "' has ", el_u.ncomp, " components, not ", D);
    if (include_pressure && el_p.ncomp != 1)
      partrac::fail(first, ": '", p_field_, "' has ", el_p.ncomp, " components, not one");
    ncoeffs_u = Uint(simplex_load::nodes_per_cell(nv, el_u.degree));
    ncoeffs_p = include_pressure ? Uint(simplex_load::nodes_per_cell(nv, el_p.degree)) : 0;
    check_dofs_fit(ncoeffs_u, ncoeffs_p, Cell::n_dofs_max, "SimplexInterpol");

    // The edges, when either field carries midside nodes
    std::vector<std::uint32_t> edges;
    std::size_t nedges = 0;
    const bool quadratic = el_u.degree == 2 || (include_pressure && el_p.degree == 2);
    if (quadratic)
      nedges = mesh_tables::build_edge_table<nv>(m.topo, m.ncells, edges);
    partrac::phase("edge table");

    // The nodes along the cells' Morton curve, where the file numbers them without locality
    std::vector<std::uint32_t> map_quad, map_lin;
    const double span = simplex_load::node_span(m, nv);
    std::cout << "Mean vertex-id span of a cell: " << span << " of the vertices" << std::endl;
    if (span > simplex_load::node_span_max){
      std::cout << "Node order: renumbering nodes along the cells' curve" << std::endl;
      if (quadratic)
        map_quad = simplex_load::morton_node_order(m, edges, nedges, nv);
      if (el_u.degree == 1 || (include_pressure && el_p.degree == 1))
        map_lin = simplex_load::morton_node_order(m, std::vector<std::uint32_t>(), 0, nv);
      partrac::phase("node order");
    }
    const std::vector<std::uint32_t>& map_u = el_u.degree == 2 ? map_quad : map_lin;
    const std::vector<std::uint32_t>& map_p = el_p.degree == 2 ? map_quad : map_lin;

    // A node and its periodic images are one node, the master's dof serving all
    simplex_load::NodePairs np;
    np.build(m, edges, nedges, nv, periodic, x_min, x_max, periodic_tol);

    load_field<nv>(first, u_field_, m, el_u, edges, nedges, map_u, np, u_a_, u_map_);
    if (include_pressure)
      load_field<nv>(first, p_field_, m, el_p, edges, nedges, map_p, np, p_a_, p_map_);
    file_a_ = first;

    const std::vector<std::uint32_t> nodes_u = np.masters(map_u, u_a_.size()/el_u.ncomp);
    u_dofs_.fill(m.topo, edges, m.ncells, nv, ne, el_u.degree == 2, m.nverts,
                 nodes_u.empty() ? nullptr : nodes_u.data());
    u_dofs_.check_stride(ncoeffs_u, "SimplexInterpol");
    if (include_pressure){
      const std::vector<std::uint32_t> nodes_p = np.masters(map_p, p_a_.size()/el_p.ncomp);
      p_dofs_.fill(m.topo, edges, m.ncells, nv, ne, el_p.degree == 2, m.nverts,
                   nodes_p.empty() ? nullptr : nodes_p.data());
      p_dofs_.check_stride(ncoeffs_p, "SimplexInterpol");
    }
    edges.clear();
    edges.shrink_to_fit();
    map_quad.clear();
    map_quad.shrink_to_fit();
    map_lin.clear();
    map_lin.shrink_to_fit();
    partrac::phase("cell dofs");

    mesh_tables::build_facet_neighbours<nv>(m.topo, m.ncells, m.coords, dim, periodic,
                                            x_min, x_max, periodic_tol, facet_neigh_);
    partrac::phase("facet table");

    hmin_ = simplex_load::shortest_edge(m, nv);
    ncells_ = m.ncells;
    nverts_ = m.nverts;
    topo_ = std::move(m.topo);
    coords_ = std::move(m.coords);

    if (use_cache){
      simplex_load::CacheTables out;
      out.topo = topo_;
      out.coords = coords_;
      out.facets = facet_neigh_;
      out.u_nodes = u_dofs_.table();
      out.u_values = u_a_;
      out.u_map = u_map_;
      out.p_nodes = p_dofs_.table();
      out.p_values = p_a_;
      out.p_map = p_map_;
      out.ncells = ncells_;
      out.nverts = nverts_;
      out.gdim = dim;
      out.ncoeffs_u = ncoeffs_u;
      out.ncoeffs_p = ncoeffs_p;
      out.ncomp_u = std::size_t(D);
      out.hmin = hmin_;
      out.x_min = x_min;
      out.x_max = x_max;
      out.stamp = file_a_;
      simplex_load::cache_write(cache_file, key, out);
      partrac::phase("mesh cache written");
    }
  }

  // The per-cell geometry and the tree, from the topology either way
  cells_.resize(ncells_);
  {
    const double* c = coords_.data();
    const std::size_t g = dim;
    const std::uint32_t* topo = topo_.data();
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t i = 0; i < std::ptrdiff_t(ncells_); ++i){
      const std::uint32_t* row = topo + std::size_t(i)*nv;
      if constexpr (nv == 4)
        cells_[i] = Cell(c + std::size_t(row[0])*g, c + std::size_t(row[1])*g,
                         c + std::size_t(row[2])*g, c + std::size_t(row[3])*g);
      else
        cells_[i] = Cell(c + std::size_t(row[0])*g, c + std::size_t(row[1])*g,
                         c + std::size_t(row[2])*g);
    }
  }
  partrac::phase("build cells");

  tree_ = std::make_unique<partrac::CellTree>(topo_.data(), ncells_, nv, coords_.data(),
                                              nverts_, dim, true);
  partrac::phase("cell tree");
  partrac::phase_total();

  std::cout << "Setting max threads: " << omp_get_max_threads() << std::endl;
}

template<typename Cell>
bool SimplexInterpol<Cell>::locate_tree(const Vector3d& xx, CellPos& pos)
{
  const int id = tree_->locate(xx);
  if (id < 0)
    return false;
  pos.id = id;
  // The exact test decided the cell; the barycentrics are the cell's own, as
  // every other path in this code computes them
  cells_[id].contains(xx, pos.bary);
  return true;
}

template<typename Cell>
void SimplexInterpol<Cell>::read_stamp(const std::string& filename, std::vector<double>& u_buf,
                                       std::vector<double>& p_buf)
{
  u_buf.assign(u_a_.size(), 0.);
  simplex_load::read_vector(filename, u_field_, u_map_, u_buf);
  if (include_pressure){
    p_buf.assign(p_a_.size(), 0.);
    simplex_load::read_vector(filename, p_field_, p_map_, p_buf);
  }
}

template<typename Cell>
void SimplexInterpol<Cell>::update(const double t)
{
  StampPair sp = ts.get(t);

  // Always load once; keep last bracket past t_max
  if ( !is_initialized || ((t_prev != sp.prev.t || t_next != sp.next.t) && t < ts.get_t_max()) ){
    partrac::phase_begin("update");
    const std::string prev_file = get_folder() + "/" + sp.prev.filename;
    const std::string next_file = get_folder() + "/" + sp.next.filename;

    // The stamp the next buffer holds may be the one wanted as previous
    if (file_b_ == prev_file && file_a_ != prev_file){
      std::cout << "Prev: Timestep = " << sp.prev.t << ", swapping... " << std::endl;
      u_a_.swap(u_b_);
      p_a_.swap(p_b_);
      file_a_.swap(file_b_);
    }
    if (file_a_ != prev_file){
      std::cout << "Prev: Timestep = " << sp.prev.t << ", filename = " << sp.prev.filename << std::endl;
      read_stamp(prev_file, u_a_, p_a_);
      file_a_ = prev_file;
    }
    u_prev_ = u_a_.data();
    p_prev_ = p_a_.data();

    std::cout << "Next: Timestep = " << sp.next.t << ", filename = " << sp.next.filename << std::endl;
    // Single stamp: the same values, not a copy of them
    if (next_file == prev_file){
      u_next_ = u_prev_;
      p_next_ = p_prev_;
    }
    else {
      if (file_b_ != next_file){
        read_stamp(next_file, u_b_, p_b_);
        file_b_ = next_file;
      }
      u_next_ = u_b_.data();
      p_next_ = p_b_.data();
    }

    partrac::phase("vector read");
    partrac::phase_total();

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  t_update = t;
}

template<typename Cell>
void SimplexInterpol<Cell>::evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  evaluate_impl<true>(x, t, pos, fields);
}

template<typename Cell>
void SimplexInterpol<Cell>::evaluate_motion(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  evaluate_impl<false>(x, t, pos, fields);
}

template<typename Cell>
template<bool Scalars>
void SimplexInterpol<Cell>::evaluate_impl(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  // Assuming inside fluid
  assert(t <= t_next && t >= t_prev);
  const double alpha_t = stamp_weight(t, t_prev, t_next);

  // Compute Pk-Pl basis at x
  const int id = pos.id;
  std::array<double, Cell::n_dofs_max> _Nu_, _Np_, _Nux_, _Nuy_, _Nuz_;   // _Nuz_ unused in 2D

  cell_basis(cells_[id], pos.bary, ncoeffs_u, _Nu_.data(), "u");
  if constexpr (Scalars)
    if (include_pressure)
      cell_basis(cells_[id], pos.bary, ncoeffs_p, _Np_.data(), "p");

  // Gathered by node: the D components of a node are consecutive
  std::array<double, Cell::n_dofs_max*3> u_prev_block, u_next_block;
  gather_stamps_nodes<Cell::n_verts, Cell::n_dofs_max, D>(u_dofs_[id], u_dofs_.stride(),
                      u_prev_, u_next_, u_prev_block.data(), u_next_block.data());

  // Evaluate
  const Vector3d U_prev = block_value<D>(_Nu_.data(), u_prev_block.data(), ncoeffs_u);
  const Vector3d U_next = block_value<D>(_Nu_.data(), u_next_block.data(), ncoeffs_u);

  // Update
  fields.U = alpha_t * U_next + (1-alpha_t) * U_prev;
  fields.A = stamp_rate(U_next, U_prev, t_prev, t_next);

  if constexpr (Scalars){
    if (include_pressure){
      std::array<double, Cell::n_dofs_max> p_prev_block, p_next_block;
      gather_stamps_nodes<Cell::n_verts, Cell::n_dofs_max, 1>(p_dofs_[id], p_dofs_.stride(),
                          p_prev_, p_next_, p_prev_block.data(), p_next_block.data());
      const double P_prev = block_scalar(_Np_.data(), p_prev_block.data(), ncoeffs_p);
      const double P_next = block_scalar(_Np_.data(), p_next_block.data(), ncoeffs_p);
      fields.P = alpha_t * P_next + (1-alpha_t) * P_prev;
    }
  }

  if (wants_gradient()){
    cell_deriv(cells_[id], pos.bary, ncoeffs_u, _Nux_.data(), _Nuy_.data(), _Nuz_.data(), "u");
    const Matrix3d gradU_prev = block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(), u_prev_block.data(), ncoeffs_u);
    const Matrix3d gradU_next = block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(), u_next_block.data(), ncoeffs_u);
    fields.gradU = alpha_t * gradU_next + (1-alpha_t) * gradU_prev;
    fields.gradA = stamp_rate(gradU_next, gradU_prev, t_prev, t_next);
  }
}

template class SimplexInterpol<Triangle>;
template class SimplexInterpol<Tet>;
