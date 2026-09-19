#ifdef USE_DOLFIN
#include "Error.hpp"
#include "XDMFInterpol.hpp"
#include "loader_params.hpp"
#include "p12_eval.hpp"
#include "geometry.hpp"
#include <array>
#include <algorithm>
#include <limits>
#include <map>
#include <numeric>
#include <cassert>
#include "dolfin_elements/P1_2.h"
#include "dolfin_elements/vP1_2.h"
#include "dolfin_elements/P1_3.h"
#include "dolfin_elements/vP1_3.h"
#include "PeriodicBC.hpp"
#include "dolfin_helpers.hpp"
#include "xdmf_helpers.hpp"

// The element spaces, per dimension
template<>
void XDMFInterpol<Triangle>::make_spaces(std::shared_ptr<const dolfin::SubDomain> constrained_domain){
  u_space_ = std::make_shared<vP1_2::FunctionSpace>(mesh, constrained_domain);
  ncoeffs_u = 3;
  p_space_ = std::make_shared<P1_2::FunctionSpace>(mesh, constrained_domain);
  ncoeffs_p = 3;
}

template<>
void XDMFInterpol<Tet>::make_spaces(std::shared_ptr<const dolfin::SubDomain> constrained_domain){
  u_space_ = std::make_shared<vP1_3::FunctionSpace>(mesh, constrained_domain);
  ncoeffs_u = 4;
  p_space_ = std::make_shared<P1_3::FunctionSpace>(mesh, constrained_domain);
  ncoeffs_p = 4;
}

template<typename Cell>
XDMFInterpol<Cell>::XDMFInterpol(const std::string& infilename)
  : MeshInterpol<Cell>(infilename)
{

  // Input file (e.g. dolfin_params.dat)
  dolfin_params = partrac::parse_file_or_exit(xdmf_schema(D == 2 ? "xdmftriangle" : "xdmftet"), infilename);

  std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  // Set periodicity
  read_mesh_params();
  if constexpr (D == 3)
    periodic_tol = 1e-4;   // the XDMF tet meshes need a looser match
  include_phi = dolfin_params.template get<bool>("include_phi");

  // Make the mesh
  dolfin::Mesh mesh_in;

  // find velocity file
  std::string xdmffilename_u = get_folder() + "/" + dolfin_params.template get<std::string>("u");
  std::string xdmffilename_p = get_folder() + "/" + dolfin_params.template get<std::string>("p");
  std::string xdmffilename_phi = get_folder() + "/" + (dolfin_params.has("phi") ? dolfin_params.template get<std::string>("phi") : "");
  //if (include_phi)
  //std::string xdmffilename_phi = get_folder() + "/" + dolfin_params["phi"];

  std::vector<std::pair<double, std::vector<std::string>>> titems_u, titems_p, titems_phi;

  std::string topology_path, geometry_path;
  titems_u = parse_xdmf(xdmffilename_u, h5filename_u, topology_path, geometry_path);

  std::cout << "mesh: " << h5filename_u << ": " << topology_path << " " << geometry_path << std::endl; 
  
  ts.initialize(titems_u);

  if (include_pressure){
    titems_p = parse_xdmf(xdmffilename_p);
    // std::cout << "h5filename_p = " << h5filename_p << std::endl;
    ts.add("p", titems_p);
  }
  if (include_phi){
    titems_phi = parse_xdmf(xdmffilename_phi);
    // std::cout << "h5filename_phi = " << h5filename_phi << std::endl;
    ts.add("phi", titems_phi);
  }

  dolfin::HDF5File meshfile(MPI_COMM_WORLD, h5filename_u, "r");

  std::string cell_type_str = XDMFCell<Cell>::name;
  Uint gdim = D;

  std::unique_ptr<dolfin::CellType> cell_type(dolfin::CellType::create(cell_type_str));

  std::vector<std::int64_t> coords_shape = dolfin::HDF5Interface::get_dataset_shape(meshfile.h5_id(), geometry_path);

  meshfile.read(mesh_in, topology_path, geometry_path, gdim, *cell_type, -1, coords_shape[0], false);

  mesh = std::make_shared<dolfin::Mesh>(mesh_in);
  init_mesh_geometry();
  assert(gdim == dim);

  auto constrained_domain = std::make_shared<PeriodicBC>(periodic, x_min, x_max, dim);
  std::cout << "Made periodic domain." << std::endl;
  
  make_spaces(constrained_domain);

  // Precompute all cells
  // FIXME compute on the fly and save

  build_cells(*u_space_->dofmap());

  // Identify edge cells
  label_cell_type(cell_type_, facet_neigh_, Cell::n_verts);

  // P2 near walls: edge (default) or none
  wall_p2_ = dolfin_params.template get<std::string>("wall_p2") == "none" ? WallP2::None : WallP2::Edge;

  std::cout << "Built neighbour list" << std::endl;
  
  // const dolfin::GenericDofMap& u_dofmap = *u_space_->dofmap();
  auto xdof = p_space_->tabulate_dof_coordinates();

  std::map<std::array<double, D>, Uint> xmap;
  for ( Uint j=0; j < xdof.size()/dim; ++j ){
    std::array<double, D> key;
    for (int c = 0; c < D; ++c) key[c] = xdof[dim*j + c];
    xmap[key] = j;
  }

  std::vector<double> xdata;
  read_dataset_vector(h5filename_u, geometry_path, xdata, dim);

  // Periodic images of a dof's vertex are not in xmap
  const Uint unset = std::numeric_limits<Uint>::max();
  j2i.assign(xdof.size()/dim, unset);
  for ( Uint i=0; i < xdata.size()/dim; ++i ){
    std::array<double, D> key;
    for (int c = 0; c < D; ++c) key[c] = xdata[dim*i + c];
    const auto it = xmap.find(key);
    if (it != xmap.end())
      j2i[it->second] = i;
  }
  if (std::find(j2i.begin(), j2i.end(), unset) != j2i.end()){
    partrac::fail("XDMFInterpol: a dof has no vertex in ", h5filename_u);
  }

  u_prev_data_.resize(xdof.size());
  u_next_data_.resize(xdof.size());
  p_prev_data_.resize(xdof.size()/dim);
  p_next_data_.resize(xdof.size()/dim);
  phi_prev_data_.resize(xdof.size()/dim);
  phi_next_data_.resize(xdof.size()/dim);

  check_dofs_fit(ncoeffs_u, ncoeffs_p, Cell::n_dofs_max, "XDMFInterpol");

  std::cout << "Setting max threads: " << omp_get_max_threads() << std::endl;

  // Precomputing dofs

  const dolfin::GenericDofMap& u_dofmap = *u_space_->dofmap();
  const dolfin::GenericDofMap& p_dofmap = *p_space_->dofmap();
  
  u_dofs_.build(u_dofmap, dolfin_cells_, "XDMFInterpol");
  u_dofs_.check_stride(D*ncoeffs_u, "XDMFInterpol");

  p_dofs_.build(p_dofmap, dolfin_cells_, "XDMFInterpol");
  p_dofs_.check_stride(ncoeffs_p, "XDMFInterpol");

  if (wall_p2_ == WallP2::Edge)
    build_wall_edges(periodic_tol);
}

template<typename Cell>
void XDMFInterpol<Cell>::update(const double t)
{
  MultiStampPair sp = ts.get(t);

  // Always load once; keep last bracket past t_max
  if ( !is_initialized || ((t_prev != sp.prev.t || t_next != sp.next.t) && t < ts.get_t_max()) )
  {
    std::vector<double> data_;
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
      read_dataset_vector(u_path_prev[0], u_path_prev[1], data_, dim);
      reorder_indices(u_prev_data_, data_, j2i, dim);
      if (wall_p2_ == WallP2::Edge)
        rest_tol_prev_ = rest_tol(u_prev_data_);

      if (include_pressure){
        auto p_path_prev = ts.get_path("p", sp.prev.it);
        read_dataset_scalar(p_path_prev[0], p_path_prev[1], data_);
        reorder_indices(p_prev_data_, data_, j2i, 1);
      }

      if (include_phi){
        auto phi_path_prev = ts.get_path("phi", sp.prev.it);
        read_dataset_scalar(phi_path_prev[0], phi_path_prev[1], data_);
        reorder_indices(phi_prev_data_, data_, j2i, 1);
      }
    }

    auto u_path_next = ts.get_path("u", sp.next.it);
    std::cout << "Next: Timestep = " << sp.next.t << ", file = " << u_path_next[0] << ":" << u_path_next[1] << std::endl;
    read_dataset_vector(u_path_next[0], u_path_next[1], data_, dim);
    reorder_indices(u_next_data_, data_, j2i, dim);

    if (include_pressure){
      auto p_path_next = ts.get_path("p", sp.next.it);
      read_dataset_scalar(p_path_next[0], p_path_next[1], data_);
      reorder_indices(p_next_data_, data_, j2i, 1);
    }

    if (include_phi){
      auto phi_path_next = ts.get_path("phi", sp.next.it);
      read_dataset_scalar(phi_path_next[0], phi_path_next[1], data_);
      reorder_indices(phi_next_data_, data_, j2i, 1);
    }

    if (wall_p2_ == WallP2::Edge)
      rest_tol_next_ = rest_tol(u_next_data_);

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  t_update = t;

}

template<typename Cell>
double XDMFInterpol<Cell>::rest_tol(const std::vector<double>& u_data) const
{
  // Round-off of the largest velocity
  double u_max = 0.;
  for (const double u : u_data) u_max = std::max(u_max, std::abs(u));
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

  // Restrict solution to cell
  std::array<double, Cell::n_dofs_max*3> u_prev_block, u_next_block;
  gather_stamps<D*Cell::n_verts, D*Cell::n_dofs_max>(u_dofs_[id], u_dofs_.stride(), u_prev_data_, u_next_data_,
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
      gather_stamps<Cell::n_verts, Cell::n_dofs_max>(p_dofs_[id], p_dofs_.stride(), p_prev_data_, p_next_data_,
                    p_prev_block.data(), p_next_block.data());
      const double P_prev = block_scalar(_Np_.data(), p_prev_block.data(), ncoeffs_p);
      const double P_next = block_scalar(_Np_.data(), p_next_block.data(), ncoeffs_p);
      fields.P = _alpha_t * P_next + (1-_alpha_t) * P_prev;
    }

    if (include_phi){
      std::array<double, Cell::n_dofs_max> phi_prev_block, phi_next_block;
      gather_stamps<Cell::n_verts, Cell::n_dofs_max>(p_dofs_[id], p_dofs_.stride(), phi_prev_data_, phi_next_data_,
                    phi_prev_block.data(), phi_next_block.data());
      const double Phi_prev = block_scalar(_Np_.data(), phi_prev_block.data(), ncoeffs_p);
      const double Phi_next = block_scalar(_Np_.data(), phi_next_block.data(), ncoeffs_p);
      fields.Phi = _alpha_t * Phi_next + (1-_alpha_t) * Phi_prev;
    }

    fields.cell_type = cell_type_[id];
  }
}

template<>
void XDMFInterpol<Triangle>::build_wall_edges(const double tol)
{
  // Periodic images of a vertex share its P1 dof
  const dolfin::GenericDofMap& p_dofmap = *p_space_->dofmap();
  std::vector<std::size_t> vclass(mesh->num_vertices());
  for (Uint l = 0; l < mesh->num_cells(); ++l){
    const auto* vi = dolfin_cells_[l].entities(0);
    const auto dofs = p_dofmap.cell_dofs(dolfin_cells_[l].index());
    for (int k = 0; k < 3; ++k) vclass[vi[k]] = dofs[k];
  }
  const std::size_t nv = p_dofmap.global_dimension();

  // Wall normals at vertices, from non-periodic exterior facets
  std::vector<Vector3d> n_sum(nv, Vector3d::Zero());
  std::vector<Vector3d> n_first(nv, Vector3d::Zero());
  std::vector<std::uint8_t> n_count(nv, 0);
  std::vector<bool> corner(nv, false);
  std::vector<Vector3d> facet_normal(mesh->num_facets(), Vector3d::Zero());
  for (dolfin::FacetIterator f(*mesh); !f.end(); ++f){
    if (!f->exterior()) continue;
    const Vector3d pt(f->midpoint().coordinates());
    bool periodic_facet = false;
    for (Uint k = 0; k < dim; ++k)
      if (periodic[k] && (pt[k] < x_min[k] + tol || pt[k] > x_max[k] - tol))
        periodic_facet = true;
    if (periodic_facet) continue;
    const Vector3d n(f->normal().coordinates());
    facet_normal[f->index()] = n;
    for (dolfin::VertexIterator v(*f); !v.end(); ++v){
      const std::size_t iv = vclass[v->index()];
      if (n_count[iv] == 0) n_first[iv] = n;
      else if (n_count[iv] > 1 || n_first[iv].dot(n) < 0.5) corner[iv] = true;   // over 60 degrees
      n_sum[iv] += n;
      ++n_count[iv];
    }
  }

  // Side edges of wall cells: sum of facet normals, per fluid end and wall end
  const std::array<std::array<int, 2>, 3> edge_ends = {{{0, 1}, {0, 2}, {1, 2}}};
  const auto opposite_edges = [&](const dolfin::Cell& cell){
    // ei[k]: the edge opposite vertex k
    const auto* vi = cell.entities(0);
    const auto* ce = cell.entities(1);
    std::array<std::size_t, 3> ei{};
    for (int j = 0; j < 3; ++j){
      const auto* ev = dolfin::Edge(*mesh, ce[j]).entities(0);
      for (int k = 0; k < 3; ++k)
        if (ev[0] != vi[k] && ev[1] != vi[k]) ei[k] = ce[j];
    }
    return ei;
  };
  std::map<std::pair<std::size_t, std::size_t>, Vector3d> side_normal;
  for (Uint l = 0; l < mesh->num_cells(); ++l){
    const auto* vi = dolfin_cells_[l].entities(0);
    const auto ei = opposite_edges(dolfin_cells_[l]);
    for (int k = 0; k < 3; ++k){
      const Vector3d& n = facet_normal[ei[k]];
      if (n.isZero()) continue;
      for (int j = 0; j < 3; ++j){
        if (j == k) continue;
        side_normal.emplace(std::make_pair(vclass[vi[k]], vclass[vi[3 - j - k]]),
                            Vector3d::Zero().eval()).first->second += n;
      }
    }
  }

  wall_index_.assign(mesh->num_cells(), -1);
  wall_cells_.clear();
  for (Uint l = 0; l < mesh->num_cells(); ++l){
    const auto* vi = dolfin_cells_[l].entities(0);
    WallEdges w{};
    for (int k = 0; k < 3; ++k)
      if (n_count[vclass[vi[k]]] > 0) w.wall |= 1 << k;
    if (w.wall == 0) continue;

    for (int e = 0; e < 3; ++e){
      for (int o = 0; o < 2; ++o){
        WallEnd& we = w.ends[2*e + o];
        const std::size_t iw = vi[edge_ends[e][o]];
        const std::size_t iv = vi[edge_ends[e][1-o]];
        const std::size_t cw = vclass[iw];
        we = {0.5, 0., 0., 0.5};   // linear
        if (n_count[cw] == 0 || corner[cw]) continue;
        const Vector3d edge = Vector3d(dolfin::Vertex(*mesh, iv).point().coordinates())
                            - Vector3d(dolfin::Vertex(*mesh, iw).point().coordinates());
        // over two wall facets: their mean normal
        const auto side = side_normal.find({vclass[iv], cw});
        const Vector3d n = (side != side_normal.end() ? side->second : n_sum[cw]).normalized();
        const double delta = edge.dot(n);
        // v_n ~ delta^2, divergence-free wall cell
        Vector3d q = -0.25*n;
        if (std::abs(delta) > 0.1*edge.norm())
          q += (edge - delta*n)/(2*delta);
        // M = I/2 + q n^T
        we = {0.5 + q[0]*n[0], q[0]*n[1], q[1]*n[0], 0.5 + q[1]*n[1]};
      }
    }
    wall_index_[l] = wall_cells_.size();
    wall_cells_.push_back(w);
  }
  std::cout << "Wall cells: " << wall_cells_.size() << std::endl;
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

  const std::array<std::array<int, 2>, 3> edge_ends = {{{0, 1}, {0, 2}, {1, 2}}};
  for (int k = 0; k < 3; ++k){
    const bool r = rest >> k & 1;
    u2[k] = r ? 0. : u[k];
    u2[6 + k] = r ? 0. : u[ncoeffs_u + k];
  }
  for (int e = 0; e < 3; ++e){
    const int a = edge_ends[e][0], c = edge_ends[e][1];
    const int m = Triangle::mid_[e];
    const bool wa = rest >> a & 1, wc = rest >> c & 1;
    const WallEnd& we = w.ends[2*e + (wa ? 0 : 1)];
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
void XDMFInterpol<Tet>::build_wall_edges(const double tol)
{
  // Periodic images of a vertex share its P1 dof
  const dolfin::GenericDofMap& p_dofmap = *p_space_->dofmap();
  std::vector<std::size_t> vclass(mesh->num_vertices());
  for (Uint l = 0; l < mesh->num_cells(); ++l){
    const auto* vi = dolfin_cells_[l].entities(0);
    const auto dofs = p_dofmap.cell_dofs(dolfin_cells_[l].index());
    for (int k = 0; k < 4; ++k) vclass[vi[k]] = dofs[k];
  }
  const std::size_t nv = p_dofmap.global_dimension();

  // Wall normals at vertices, from non-periodic exterior facets
  std::vector<Vector3d> n_sum(nv, Vector3d::Zero());
  std::vector<Vector3d> n_first(nv, Vector3d::Zero());
  std::vector<std::uint8_t> n_count(nv, 0);
  std::vector<bool> corner(nv, false);
  std::vector<Vector3d> facet_normal(mesh->num_facets(), Vector3d::Zero());
  for (dolfin::FacetIterator f(*mesh); !f.end(); ++f){
    if (!f->exterior()) continue;
    const Vector3d pt(f->midpoint().coordinates());
    bool periodic_facet = false;
    for (Uint k = 0; k < dim; ++k)
      if (periodic[k] && (pt[k] < x_min[k] + tol || pt[k] > x_max[k] - tol))
        periodic_facet = true;
    if (periodic_facet) continue;
    const Vector3d n(f->normal().coordinates());
    facet_normal[f->index()] = n;
    for (dolfin::VertexIterator v(*f); !v.end(); ++v){
      const std::size_t iv = vclass[v->index()];
      if (n_count[iv] == 0) n_first[iv] = n;
      else if (n_first[iv].dot(n) < 0.5) corner[iv] = true;   // over 60 degrees
      n_sum[iv] += n;
      ++n_count[iv];
    }
  }

  // Side edges of wall cells: sum of facet normals, per fluid end and wall end
  std::map<std::pair<std::size_t, std::size_t>, Vector3d> side_normal;
  for (Uint l = 0; l < mesh->num_cells(); ++l){
    const auto* vi = dolfin_cells_[l].entities(0);
    const auto* fi = dolfin_cells_[l].entities(2);
    for (int j = 0; j < 4; ++j){
      const Vector3d& n = facet_normal[fi[j]];
      if (n.isZero()) continue;
      // the apex is the vertex off the facet
      const auto* fv = dolfin::Facet(*mesh, fi[j]).entities(0);
      for (int k = 0; k < 4; ++k){
        if (vi[k] == fv[0] || vi[k] == fv[1] || vi[k] == fv[2]) continue;
        for (int i = 0; i < 3; ++i)
          side_normal.emplace(std::make_pair(vclass[vi[k]], vclass[fv[i]]),
                              Vector3d::Zero().eval()).first->second += n;
      }
    }
  }

  const std::array<std::array<int, 2>, 6> edge_ends = {{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
  wall_index_.assign(mesh->num_cells(), -1);
  wall_cells_.clear();
  for (Uint l = 0; l < mesh->num_cells(); ++l){
    const auto* vi = dolfin_cells_[l].entities(0);
    WallEdges w{};
    for (int k = 0; k < 4; ++k)
      if (n_count[vclass[vi[k]]] > 0) w.wall |= 1 << k;
    if (w.wall == 0) continue;

    for (int e = 0; e < 6; ++e){
      for (int o = 0; o < 2; ++o){
        WallEnd& we = w.ends[2*e + o];
        const std::size_t iw = vi[edge_ends[e][o]];
        const std::size_t iv = vi[edge_ends[e][1-o]];
        const std::size_t cw = vclass[iw];
        we = {0., 0., 0., 0., 0., 0.};   // linear
        if (n_count[cw] == 0 || corner[cw]) continue;
        const Vector3d edge = Vector3d(dolfin::Vertex(*mesh, iv).point().coordinates())
                            - Vector3d(dolfin::Vertex(*mesh, iw).point().coordinates());
        // over several wall facets: their mean normal
        const auto side = side_normal.find({vclass[iv], cw});
        const Vector3d n = (side != side_normal.end() ? side->second : n_sum[cw]).normalized();
        const double delta = edge.dot(n);
        // v_n ~ delta^2, divergence-free wall cell
        Vector3d q = -0.25*n;
        if (std::abs(delta) > 0.1*edge.norm())
          q += (edge - delta*n)/(4*delta);
        we = {q[0], q[1], q[2], n[0], n[1], n[2]};
      }
    }
    wall_index_[l] = wall_cells_.size();
    wall_cells_.push_back(w);
  }
  std::cout << "Wall cells: " << wall_cells_.size() << std::endl;
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

  const std::array<std::array<int, 2>, 6> edge_ends = {{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
  for (int k = 0; k < 4; ++k){
    const bool r = rest >> k & 1;
    for (int c = 0; c < 3; ++c)
      u2[10*c + k] = r ? 0. : u[c*ncoeffs_u + k];
  }
  for (int e = 0; e < 6; ++e){
    const int a = edge_ends[e][0], b = edge_ends[e][1];
    const int m = Tet::mid_[e];
    const bool wa = rest >> a & 1, wb = rest >> b & 1;
    const WallEnd& we = w.ends[2*e + (wa ? 0 : 1)];
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

#endif
