#include "Error.hpp"
#include "XDMFInterpol.hpp"
#include "loader_params.hpp"
#include "mesh_tables.hpp"
#include "p12_eval.hpp"
#include "phase_timing.hpp"
#include "simplex_load.hpp"
#include "xdmf_helpers.hpp"
#include <array>
#include <algorithm>
#include <cassert>

template<typename Cell>
XDMFInterpol<Cell>::XDMFInterpol(const std::string& infilename)
  : MeshCore<Cell>(infilename)
{
  constexpr int nv = Cell::n_verts;
  constexpr int ne = nv*(nv-1)/2;
  partrac::phase_begin("load");

  // Input file (e.g. dolfin_params.dat)
  dolfin_params = partrac::parse_file_or_exit(xdmf_schema(D == 2 ? "xdmftriangle" : "xdmftet"), infilename);

  const std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  // Set periodicity
  read_mesh_params();
  if constexpr (D == 3)
    periodic_tol = 1e-4;   // the XDMF tet meshes need a looser match
  include_phi = dolfin_params.template get<bool>("include_phi");
  // P2 near walls: edge (default) or none
  wall_p2_ = dolfin_params.template get<std::string>("wall_p2") == "none" ? WallP2::None : WallP2::Edge;

  // One xdmf per field; the velocity's carries the mesh
  const std::string xdmffilename_u = get_folder() + "/" + dolfin_params.template get<std::string>("u");
  const std::string xdmffilename_p = get_folder() + "/"
    + (dolfin_params.has("p") ? dolfin_params.template get<std::string>("p") : "");
  const std::string xdmffilename_phi = get_folder() + "/"
    + (dolfin_params.has("phi") ? dolfin_params.template get<std::string>("phi") : "");

  std::string h5filename_u, topology_path, geometry_path;
  ts.initialize(parse_xdmf(xdmffilename_u, h5filename_u, topology_path, geometry_path));

  std::cout << "mesh: " << h5filename_u << ": " << topology_path << " " << geometry_path << std::endl;

  if (include_pressure)
    ts.add("p", parse_xdmf(xdmffilename_p));
  if (include_phi)
    ts.add("phi", parse_xdmf(xdmffilename_phi));
  partrac::phase("params");

  // The mesh the xdmf names, its cells along the Morton curve
  simplex_load::MeshData m;
  std::vector<std::uint32_t> cell_perm;
  simplex_load::read_mesh_arrays(h5filename_u, topology_path, geometry_path, nv, m, cell_perm);
  cell_perm.clear();
  cell_perm.shrink_to_fit();
  dim = m.gdim;
  x_min = m.x_min;
  x_max = m.x_max;
  ncells_ = m.ncells;
  nverts_ = m.nverts;
  set_period();

  // A periodic space gave a vertex and its images one dof, so they read one value
  vclass_ = mesh_tables::match_periodic_vertices(m.coords, m.nverts, dim, periodic,
                                                 x_min, x_max, periodic_tol);
  partrac::phase("periodic vertices");

  // The dofs are the vertices: no dofmap, and one node table for every field
  ncoeffs_u = Uint(nv);
  ncoeffs_p = Uint(nv);
  check_dofs_fit(ncoeffs_u, ncoeffs_p, Cell::n_dofs_max, "XDMFInterpol");
  u_dofs_.fill(m.topo, std::vector<std::uint32_t>(), m.ncells, nv, ne, false, m.nverts,
               vclass_.data());
  u_dofs_.check_stride(ncoeffs_u, "XDMFInterpol");
  partrac::phase("cell dofs");

  mesh_tables::build_facet_neighbours<nv>(m.topo, m.ncells, m.coords, dim, periodic,
                                          x_min, x_max, periodic_tol, facet_neigh_);
  partrac::phase("facet table");

  // Identify edge cells
  near_wall::label_cell_type(cell_type_, facet_neigh_, nv);
  std::cout << "Built neighbour list" << std::endl;

  if (wall_p2_ == WallP2::Edge){
    near_wall::build_wall_edges<Cell>(m.topo, m.coords, ncells_, nverts_, dim, facet_neigh_,
                                      vclass_, wall_index_, wall_cells_, verbose);
    partrac::phase("wall edges");
  }

  hmin_ = simplex_load::shortest_edge(m, nv);
  topo_ = std::move(m.topo);
  coords_ = std::move(m.coords);

  // The per-cell geometry and the tree, from the topology
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

  u_prev_data_.resize(nverts_*dim);
  u_next_data_.resize(nverts_*dim);
  p_prev_data_.resize(nverts_);
  p_next_data_.resize(nverts_);
  phi_prev_data_.resize(nverts_);
  phi_next_data_.resize(nverts_);

  std::cout << "Setting max threads: " << omp_get_max_threads() << std::endl;
}

template<typename Cell>
bool XDMFInterpol<Cell>::locate_tree(const Vector3d& xx, CellPos& pos)
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
void XDMFInterpol<Cell>::read_stamp(const std::vector<std::string>& path, std::vector<double>& data,
                                    const int ncols)
{
  read_dataset_columns(path[0], path[1], data, ncols);
  if (data.size() != nverts_*std::size_t(ncols)){
    partrac::fail(path[0], ": '", path[1], "' holds ", data.size()/std::size_t(ncols),
                  " rows, the mesh has ", nverts_, " vertices");
  }
}

template<typename Cell>
void XDMFInterpol<Cell>::update(const double t)
{
  MultiStampPair sp = ts.get(t);

  // Always load once; keep last bracket past t_max
  if ( !is_initialized || ((t_prev != sp.prev.t || t_next != sp.next.t) && t < ts.get_t_max()) )
  {
    partrac::phase_begin("update");
    // Swap if possible
    if (is_initialized && t_next == sp.prev.t)
    {
      std::cout << "Prev: Timestep = " << sp.prev.t << ", swapping... "<< std::endl;
      u_prev_data_.swap(u_next_data_);
      rest_tol_prev_ = rest_tol_next_;
      if (include_pressure)
        p_prev_data_.swap(p_next_data_);
      if (include_phi)
        phi_prev_data_.swap(phi_next_data_);
    }
    else
    {
      auto u_path_prev = ts.get_path("u", sp.prev.it);
      std::cout << "Prev: Timestep = " << sp.prev.t << ", file = " << u_path_prev[0] << ":" << u_path_prev[1] << std::endl;
      read_stamp(u_path_prev, u_prev_data_, int(dim));
      if (wall_p2_ == WallP2::Edge)
        rest_tol_prev_ = rest_tol(u_prev_data_);

      if (include_pressure)
        read_stamp(ts.get_path("p", sp.prev.it), p_prev_data_, 1);

      if (include_phi)
        read_stamp(ts.get_path("phi", sp.prev.it), phi_prev_data_, 1);
    }

    auto u_path_next = ts.get_path("u", sp.next.it);
    std::cout << "Next: Timestep = " << sp.next.t << ", file = " << u_path_next[0] << ":" << u_path_next[1] << std::endl;
    read_stamp(u_path_next, u_next_data_, int(dim));

    if (include_pressure)
      read_stamp(ts.get_path("p", sp.next.it), p_next_data_, 1);

    if (include_phi)
      read_stamp(ts.get_path("phi", sp.next.it), phi_next_data_, 1);

    if (wall_p2_ == WallP2::Edge)
      rest_tol_next_ = rest_tol(u_next_data_);

    partrac::phase("vector read");
    partrac::phase_total();

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  t_update = t;

}

template<typename Cell>
double XDMFInterpol<Cell>::rest_tol(const std::vector<double>& u_data) const
{
  // Round-off of the largest velocity a cell can read: an image vertex's own
  // value is never gathered, its master's is
  double u_max = 0.;
  for (std::size_t v = 0; v < nverts_; ++v){
    if (vclass_[v] != v) continue;
    for (Uint c = 0; c < dim; ++c) u_max = std::max(u_max, std::abs(u_data[v*dim + c]));
  }
  return 1e-12*u_max;
}

template<typename Cell>
void XDMFInterpol<Cell>::evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  evaluate_impl<true>(x, t, pos, fields);
}

template<typename Cell>
void XDMFInterpol<Cell>::evaluate_motion(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  evaluate_impl<false>(x, t, pos, fields);
}

template<typename Cell>
template<bool Scalars>
void XDMFInterpol<Cell>::evaluate_impl(const Vector3d &x, const double tin, const CellPos& pos, PointValues& fields)
{
  const double _alpha_t = stamp_weight(tin, t_prev, t_next);

  // Compute P1 basis at x
  const int id = pos.id;
  std::array<double, Cell::n_dofs_max> _Nu_, _Np_, _Nux_, _Nuy_, _Nuz_;   // _Nuz_ unused in 2D

  cell_basis(cells_[id], pos.bary, ncoeffs_u, _Nu_.data(), "u");
  if constexpr (Scalars)
    if (include_pressure)
      cell_basis(cells_[id], pos.bary, ncoeffs_p, _Np_.data(), "p");

  // Restrict solution to cell: gathered by vertex, the D components consecutive
  std::array<double, Cell::n_dofs_max*3> u_prev_block, u_next_block;
  gather_stamps_nodes<Cell::n_verts, Cell::n_dofs_max, D>(u_dofs_[id], u_dofs_.stride(),
                      u_prev_data_.data(), u_next_data_.data(),
                      u_prev_block.data(), u_next_block.data());

  // Evaluate
  Vector3d U_prev = block_value<D>(_Nu_.data(), u_prev_block.data(), ncoeffs_u);
  Vector3d U_next = block_value<D>(_Nu_.data(), u_next_block.data(), ncoeffs_u);

  // Gradient
  Matrix3d gradU_prev = Matrix3d::Zero(), gradU_next = Matrix3d::Zero();
  if (wants_gradient()){
    cell_deriv(cells_[id], pos.bary, ncoeffs_u, _Nux_.data(), _Nuy_.data(), _Nuz_.data(), "u");
    gradU_prev = block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(), u_prev_block.data(), ncoeffs_u);
    gradU_next = block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(), u_next_block.data(), ncoeffs_u);
  }

  // Quadratic near walls
  std::array<double, D*Cell::n_dofs_max> u_prev_block_2, u_next_block_2;
  bool quad_prev = false;
  bool quad_next = false;
  if (wall_p2_ == WallP2::Edge && wall_index_[id] >= 0){
    const WallEdges& w = wall_cells_[wall_index_[id]];
    quad_prev = wall_block(u_prev_block.data(), u_prev_block_2.data(), w, rest_tol_prev_);
    quad_next = wall_block(u_next_block.data(), u_next_block_2.data(), w, rest_tol_next_);
  }

  if (quad_prev || quad_next){
    const Uint n2 = Cell::n_dofs_max;
    std::array<double, Cell::n_dofs_max> _Nu2_, _Nu2x_, _Nu2y_, _Nu2z_;   // _Nu2z_ unused in 2D
    cell_basis(cells_[id], pos.bary, n2, _Nu2_.data(), "u");
    if (wants_gradient())
      cell_deriv(cells_[id], pos.bary, n2, _Nu2x_.data(), _Nu2y_.data(), _Nu2z_.data(), "u");

    if (quad_prev){
      U_prev = block_value<D>(_Nu2_.data(), u_prev_block_2.data(), n2);
      if (wants_gradient())
        gradU_prev = block_gradient<D>(_Nu2x_.data(), _Nu2y_.data(), _Nu2z_.data(), u_prev_block_2.data(), n2);
    }
    if (quad_next){
      U_next = block_value<D>(_Nu2_.data(), u_next_block_2.data(), n2);
      if (wants_gradient())
        gradU_next = block_gradient<D>(_Nu2x_.data(), _Nu2y_.data(), _Nu2z_.data(), u_next_block_2.data(), n2);
    }
  }

  // Update
  fields.U = _alpha_t * U_next + (1-_alpha_t) * U_prev;
  fields.A = stamp_rate(U_next, U_prev, t_prev, t_next);

  if (wants_gradient()){
    fields.gradU = _alpha_t * gradU_next + (1-_alpha_t) * gradU_prev;
    fields.gradA = stamp_rate(gradU_next, gradU_prev, t_prev, t_next);
  }

  if constexpr (Scalars){
    if (include_pressure){
      std::array<double, Cell::n_dofs_max> p_prev_block, p_next_block;
      gather_stamps_nodes<Cell::n_verts, Cell::n_dofs_max, 1>(u_dofs_[id], u_dofs_.stride(),
                          p_prev_data_.data(), p_next_data_.data(),
                          p_prev_block.data(), p_next_block.data());
      const double P_prev = block_scalar(_Np_.data(), p_prev_block.data(), ncoeffs_p);
      const double P_next = block_scalar(_Np_.data(), p_next_block.data(), ncoeffs_p);
      fields.P = _alpha_t * P_next + (1-_alpha_t) * P_prev;
    }

    if (include_phi){
      std::array<double, Cell::n_dofs_max> phi_prev_block, phi_next_block;
      gather_stamps_nodes<Cell::n_verts, Cell::n_dofs_max, 1>(u_dofs_[id], u_dofs_.stride(),
                          phi_prev_data_.data(), phi_next_data_.data(),
                          phi_prev_block.data(), phi_next_block.data());
      const double Phi_prev = block_scalar(_Np_.data(), phi_prev_block.data(), ncoeffs_p);
      const double Phi_next = block_scalar(_Np_.data(), phi_next_block.data(), ncoeffs_p);
      fields.Phi = _alpha_t * Phi_next + (1-_alpha_t) * Phi_prev;
    }

    fields.cell_type = cell_type_[id];
  }
}

template<>
bool XDMFInterpol<Triangle>::wall_block(const double* u, double* u2, const WallEdges& w,
                                      const double tol) const
{
  // Walls: listed vertices at rest
  unsigned rest = 0;
  for (int k = 0; k < 3; ++k)
    if ((w.wall >> k & 1) && std::abs(u[k]) <= tol && std::abs(u[ncoeffs_u + k]) <= tol)
      rest |= 1u << k;
  if (rest == 0) return false;

  constexpr auto edge_ends = near_wall::WallRule<Triangle>::edge_ends;
  for (int k = 0; k < 3; ++k){
    const bool r = rest >> k & 1;
    u2[k] = r ? 0. : u[k];
    u2[6 + k] = r ? 0. : u[ncoeffs_u + k];
  }
  for (int e = 0; e < 3; ++e){
    const int a = edge_ends[e][0], c = edge_ends[e][1];
    const int m = Triangle::mid_[e];
    const bool wa = rest >> a & 1, wc = rest >> c & 1;
    const near_wall::WallRule<Triangle>::End& we = w.ends[2*e + (wa ? 0 : 1)];
    if (wa != wc){
      const int v = wa ? c : a;
      const double ux = u[v], uy = u[ncoeffs_u + v];
      u2[m] = we.mxx*ux + we.mxy*uy;
      u2[6 + m] = we.myx*ux + we.myy*uy;
    }
    else {
      u2[m] = 0.5*(u2[a] + u2[c]);
      u2[6 + m] = 0.5*(u2[6 + a] + u2[6 + c]);
    }
  }
  return true;
}

template<>
bool XDMFInterpol<Tet>::wall_block(const double* u, double* u2, const WallEdges& w,
                                 const double tol) const
{
  // Walls: listed vertices at rest
  unsigned rest = 0;
  for (int k = 0; k < 4; ++k)
    if ((w.wall >> k & 1) && std::abs(u[k]) <= tol && std::abs(u[ncoeffs_u + k]) <= tol
        && std::abs(u[2*ncoeffs_u + k]) <= tol)
      rest |= 1u << k;
  if (rest == 0) return false;

  constexpr auto edge_ends = near_wall::WallRule<Tet>::edge_ends;
  for (int k = 0; k < 4; ++k){
    const bool r = rest >> k & 1;
    for (int c = 0; c < 3; ++c)
      u2[10*c + k] = r ? 0. : u[c*ncoeffs_u + k];
  }
  for (int e = 0; e < 6; ++e){
    const int a = edge_ends[e][0], b = edge_ends[e][1];
    const int m = Tet::mid_[e];
    const bool wa = rest >> a & 1, wb = rest >> b & 1;
    const near_wall::WallRule<Tet>::End& we = w.ends[2*e + (wa ? 0 : 1)];
    if (wa != wb){
      const int v = wa ? b : a;
      const double ux = u[v], uy = u[ncoeffs_u + v], uz = u[2*ncoeffs_u + v];
      const double un = ux*we.nx + uy*we.ny + uz*we.nz;
      u2[m] = 0.5*ux + un*we.qx;
      u2[10 + m] = 0.5*uy + un*we.qy;
      u2[20 + m] = 0.5*uz + un*we.qz;
    }
    else {
      for (int c = 0; c < 3; ++c)
        u2[10*c + m] = 0.5*(u2[10*c + a] + u2[10*c + b]);
    }
  }
  return true;
}

template class XDMFInterpol<Triangle>;
template class XDMFInterpol<Tet>;
