#ifdef USE_DOLFIN
#include "geometry.hpp"
#include "loader_params.hpp"
#include "XDMFTriangleInterpol.hpp"
#include <array>
#include <algorithm>
#include <limits>
#include <map>
#include <numeric>
#include "Timestamps.hpp"
//#include "H5Cpp.h"
#include <cassert>
#include "dolfin_elements/P1_2.h"
#include "dolfin_elements/P2_2.h"
#include "dolfin_elements/vP1_2.h"
#include "dolfin_elements/vP2_2.h"
#include "PeriodicBC.hpp"
#include "dolfin_helpers.hpp"
#include "xdmf_helpers.hpp"

//using namespace H5;

XDMFTriangleInterpol::XDMFTriangleInterpol(const std::string& infilename)
  : Interpol(infilename)
{

  // Input file (e.g. dolfin_params.dat)
  dolfin_params = partrac::parse_file_or_exit(xdmf_schema("xdmftriangle"), infilename);

  std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  // Set periodicity
  if (dolfin_params.get<bool>("periodic_x")){
    periodic[0] = true;
  }
  if (dolfin_params.get<bool>("periodic_y")){
    periodic[1] = true;
  }
  include_pressure = !dolfin_params.get<bool>("ignore_pressure");
  include_phi = dolfin_params.get<bool>("include_phi");

  // Make the mesh
  dolfin::Mesh mesh_in;

  // find velocity file
  std::string xdmffilename_u = get_folder() + "/" + dolfin_params.get<std::string>("u");
  std::string xdmffilename_p = get_folder() + "/" + dolfin_params.get<std::string>("p");
  std::string xdmffilename_phi = get_folder() + "/" + (dolfin_params.has("phi") ? dolfin_params.get<std::string>("phi") : "");
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

  std::string cell_type_str = "triangle";
  std::unique_ptr<dolfin::CellType> cell_type(dolfin::CellType::create(cell_type_str));

  std::vector<std::int64_t> coords_shape = dolfin::HDF5Interface::get_dataset_shape(meshfile.h5_id(), geometry_path);

  Uint gdim = 2;

  meshfile.read(mesh_in, topology_path, geometry_path, gdim, *cell_type, -1, coords_shape[0], false);

  mesh = std::make_shared<dolfin::Mesh>(mesh_in);
  dim = mesh->geometry().dim();
  assert(dim == gdim);
  mesh->init();
  mesh->bounding_box_tree();

  std::vector<double> xx = mesh->coordinates();

  for (Uint i=0; i<dim; ++i){
    x_min[i] = xx[i];
    x_max[i] = xx[i];
  }

  for (Uint i=0; i<xx.size(); ++i){
    Uint i_loc = i % dim;
    x_min[i_loc] = std::min(x_min[i_loc], xx[i]);
    x_max[i_loc] = std::max(x_max[i_loc], xx[i]);
  }

  auto constrained_domain = std::make_shared<PeriodicBC>(periodic, x_min, x_max, dim);
  std::cout << "Made periodic domain." << std::endl;

  // Velocity
  u_space_ = std::make_shared<vP1_2::FunctionSpace>(mesh, constrained_domain);
  ncoeffs_u = 3;

  // Pressure
  p_space_ = std::make_shared<P1_2::FunctionSpace>(mesh, constrained_domain);
  ncoeffs_p = 3;

  // Precompute all triangles Taylor-Hood P2-P1
  // FIXME compute on the fly and save
  triangles_.resize(mesh->num_cells());
  dolfin_cells_.resize(mesh->num_cells());
  cell2cells_.resize(mesh->num_cells());

  const std::vector<std::uint32_t> order =
    cell_order(*u_space_->dofmap(), mesh->num_cells(), dolfin_params.get<std::string>("renumber_cells"), dolfin2local_);
  for (std::size_t l = 0; l < mesh->num_cells(); ++l)
  {
    dolfin::Cell dolfin_cell(*mesh, order[l]);
    triangles_[l] = Triangle(dolfin_cell);
    dolfin_cells_[l] = dolfin_cell;
  }
  // Build cell neighbour list for lookup speed
  build_neighbor_list(cell2cells_, mesh, dolfin_cells_,
                      dolfin2local_.empty() ? nullptr : &dolfin2local_);

  // Identify edge cells
  cell_type_.resize(mesh->num_cells());
  //label_cell_type(cell_type_, cell2cells_, dim);

  double tol = 1e-12; // heuristic

  apply_periodic_boundaries(cell2cells_, periodic, x_min, x_max, mesh, dolfin_cells_, dim, tol);

  label_cell_type(cell_type_, cell2cells_, dim);

  // P2 near walls: edge (default) or none
  wall_p2_ = dolfin_params.get<std::string>("wall_p2") == "none" ? WallP2::None : WallP2::Edge;
  if (wall_p2_ == WallP2::Edge)
    build_wall_edges(tol);

  std::cout << "Built neighbour list" << std::endl;
  
  // const dolfin::GenericDofMap& u_dofmap = *u_space_->dofmap();
  auto xdof = p_space_->tabulate_dof_coordinates();

  std::map<std::tuple<double, double>, Uint> xmap;
  for ( Uint j=0; j < xdof.size()/2; ++j ){
    //std::cout << " " << x_dof[i] << " " << x_dof[i+1] << std::endl;
    xmap[{xdof[2*j], xdof[2*j+1]}] = j;
  }

  std::vector<double> xdata;
  read_dataset_vector(h5filename_u, geometry_path, xdata, dim);

  // Periodic images of a dof's vertex are not in xmap
  const Uint unset = std::numeric_limits<Uint>::max();
  j2i.assign(xdof.size()/2, unset);
  for ( Uint i=0; i < xdata.size()/2; ++i ){
    const auto it = xmap.find({xdata[2*i], xdata[2*i+1]});
    if (it != xmap.end())
      j2i[it->second] = i;
  }
  if (std::find(j2i.begin(), j2i.end(), unset) != j2i.end()){
    std::cout << "XDMFTriangleInterpol: a dof has no vertex in " << h5filename_u << std::endl;
    exit(1);
  }

  u_prev_data_.resize(xdof.size());
  u_next_data_.resize(xdof.size());
  p_prev_data_.resize(xdof.size()/2);
  p_next_data_.resize(xdof.size()/2);
  phi_prev_data_.resize(xdof.size()/2);
  phi_next_data_.resize(xdof.size()/2);

  //u_prev_coefficients_.resize(dim*ncoeffs_u);
  //u_next_coefficients_.resize(dim*ncoeffs_u);

  //Nu_.resize(ncoeffs_u);
  //Nux_.resize(ncoeffs_u);
  //Nuy_.resize(ncoeffs_u);

  //if (include_pressure){
    //p_prev_coefficients_.resize(ncoeffs_p); // not needed?
    //p_next_coefficients_.resize(ncoeffs_p);
    
    //Np_.resize(ncoeffs_p);
  //}

  check_dofs_fit(ncoeffs_u, ncoeffs_p, Triangle::n_dofs_max, "XDMFTriangleInterpol");

  std::cout << "Setting max threads: " << omp_get_max_threads() << std::endl;

  found_.resize(omp_get_max_threads());

  // Precomputing dofs

  const dolfin::GenericDofMap& u_dofmap = *u_space_->dofmap();
  const dolfin::GenericDofMap& p_dofmap = *p_space_->dofmap();
  
  u_dofs_.build(u_dofmap, dolfin_cells_, "XDMFTriangleInterpol");
  u_dofs_.check_stride(2*ncoeffs_u, "XDMFTriangleInterpol");

  p_dofs_.build(p_dofmap, dolfin_cells_, "XDMFTriangleInterpol");
  p_dofs_.check_stride(ncoeffs_p, "XDMFTriangleInterpol");
}

void XDMFTriangleInterpol::update(const double t)
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



void XDMFTriangleInterpol::build_wall_edges(const double tol)
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

double XDMFTriangleInterpol::rest_tol(const std::vector<double>& u_data) const
{
  // Round-off of the largest velocity
  double u_max = 0.;
  for (const double u : u_data) u_max = std::max(u_max, std::abs(u));
  return 1e-12*u_max;
}

bool XDMFTriangleInterpol::wall_block(const double* u, double* u2, const WallEdges& w,
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

Vector3d XDMFTriangleInterpol::_modx(const Vector3d &x){
  Vector3d x_loc = x;
  for (std::size_t i=0; i<dim; ++i){
    if (periodic[i]){
      x_loc[i] = x_min[i] + modulox(x[i]-x_min[i], x_max[i]-x_min[i]);
    }
    else {
      x_loc[i] = x[i];
    }
  }
  return x_loc;
}


bool XDMFTriangleInterpol::locate(const Vector3d &x, const double t, CellPos& pos)
{
  assert(t <= t_next && t >= t_prev);
  const Vector3d xx = _modx(x);
  return locate_in_cells(triangles_, cell2cells_, *mesh, dim, xx, pos,
                         found_, dolfin2local_.empty() ? nullptr : &dolfin2local_);
}

void XDMFTriangleInterpol::evaluate(const Vector3d &x, const double tin, const CellPos& pos, PointValues& fields)
{
  const double _alpha_t = stamp_weight(tin, t_prev, t_next);

  // Compute Pk-Pl basis at x
  const int id = pos.id;
  const double r1 = pos.bary[0], r2 = pos.bary[1], r3 = pos.bary[2];

  std::array<double, Triangle::n_dofs_max> _Nu_;
  std::array<double, Triangle::n_dofs_max> _Np_;
  std::array<double, Triangle::n_dofs_max> _Nux_;
  std::array<double, Triangle::n_dofs_max> _Nuy_;

  triangles_[id].linearbasis(r1, r2, r3, _Nu_.data());

  if (include_pressure){
    triangles_[id].linearbasis(r1, r2, r3, _Np_.data());
  }

  std::array<double, Triangle::n_dofs_max*3> u_prev_block;
  std::array<double, Triangle::n_dofs_max*3> u_next_block;

  // Restrict solution to cell
  //const dolfin::GenericDofMap& u_dofmap = *u_space_->dofmap();
  //auto u_dofs = u_dofmap.cell_dofs(dolfin_cells_[id].index());

  const std::uint32_t* u_dofs = u_dofs_[id];
  for (std::size_t i = 0; i < u_dofs_.stride(); ++i){
      u_prev_block[i] = u_prev_data_[u_dofs[i]];
      u_next_block[i] = u_next_data_[u_dofs[i]];
  }

  // Evaluate
  Vector3d U_prev = {std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, u_prev_block.begin(), 0.0),
                     std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, &u_prev_block[ncoeffs_u], 0.0),
                     0.0};
  Vector3d U_next = {std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, u_next_block.begin(), 0.0),
                     std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, &u_next_block[ncoeffs_u], 0.0),
                     0.0 };

  // Gradient
  Matrix3d gradU_prev = Matrix3d::Zero(), gradU_next = Matrix3d::Zero();
  if (wants_gradient()){
    triangles_[id].linearderiv(r1, r2, r3, _Nux_.data(), _Nuy_.data());

    gradU_prev <<
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, u_prev_block.begin(), 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, u_prev_block.begin(), 0.0),
      0.0,
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, &u_prev_block[ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, &u_prev_block[ncoeffs_u], 0.0),
      0.0,
      0.0,
      0.0,
      0.0;
    gradU_next << 
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, u_next_block.begin(), 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, u_next_block.begin(), 0.0),
      0.0,
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, &u_next_block[ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, &u_next_block[ncoeffs_u], 0.0),
      0.0,
      0.0,
      0.0,
      0.0;
  }

  // Quadratic near walls
  std::array<double, 12> u_prev_block_2;
  std::array<double, 12> u_next_block_2;
  bool quad_prev = false;
  bool quad_next = false;
  if (wall_p2_ == WallP2::Edge && wall_index_[id] >= 0){
    const WallEdges& w = wall_cells_[wall_index_[id]];
    quad_prev = wall_block(u_prev_block.data(), u_prev_block_2.data(), w, rest_tol_prev_);
    quad_next = wall_block(u_next_block.data(), u_next_block_2.data(), w, rest_tol_next_);
  }

  if (quad_prev || quad_next){
    const Uint ncoeffs_u_2 = 6;
    std::array<double, 6> _Nu2_;
    std::array<double, 6> _Nu2x_;
    std::array<double, 6> _Nu2y_;
    triangles_[id].quadbasis(r1, r2, r3, _Nu2_.data());
    if (wants_gradient())
      triangles_[id].quadderiv(r1, r2, r3, _Nu2x_.data(), _Nu2y_.data());

    const auto quad = [&](const std::array<double, 12>& b, Vector3d& U, Matrix3d& gradU){
      U = {std::inner_product(_Nu2_.data(), _Nu2_.data()+ncoeffs_u_2, b.begin(), 0.0),
           std::inner_product(_Nu2_.data(), _Nu2_.data()+ncoeffs_u_2, &b[ncoeffs_u_2], 0.0),
           0.0};
      if (wants_gradient()){
        gradU <<
          std::inner_product(_Nu2x_.data(), _Nu2x_.data()+ncoeffs_u_2, b.begin(), 0.0),
          std::inner_product(_Nu2y_.data(), _Nu2y_.data()+ncoeffs_u_2, b.begin(), 0.0),
          0.0,
          std::inner_product(_Nu2x_.data(), _Nu2x_.data()+ncoeffs_u_2, &b[ncoeffs_u_2], 0.0),
          std::inner_product(_Nu2y_.data(), _Nu2y_.data()+ncoeffs_u_2, &b[ncoeffs_u_2], 0.0),
          0.0,
          0.0,
          0.0,
          0.0;
      }
    };
    if (quad_prev) quad(u_prev_block_2, U_prev, gradU_prev);
    if (quad_next) quad(u_next_block_2, U_next, gradU_next);
  }

  // Update
  fields.U = _alpha_t * U_next + (1-_alpha_t) * U_prev;
  fields.A = stamp_rate(U_next, U_prev, t_prev, t_next);

  if (wants_gradient()){
    fields.gradU = _alpha_t * gradU_next + (1-_alpha_t) * gradU_prev;
    fields.gradA = stamp_rate(gradU_next, gradU_prev, t_prev, t_next);
  }

  if (include_pressure){
    std::array<double, Triangle::n_dofs_max> p_prev_block;
    std::array<double, Triangle::n_dofs_max> p_next_block;
    
    const std::uint32_t* p_dofs = p_dofs_[id];
    for (std::size_t i = 0; i < p_dofs_.stride(); ++i){
        p_prev_block[i] = p_prev_data_[p_dofs[i]];
        p_next_block[i] = p_next_data_[p_dofs[i]];
    }
    
    // Evaluate
    double P_prev = std::inner_product(_Np_.data(), _Np_.data()+ncoeffs_p, p_prev_block.begin(), 0.0);
    double P_next = std::inner_product(_Np_.data(), _Np_.data()+ncoeffs_p, p_next_block.begin(), 0.0);
    fields.P = _alpha_t * P_next + (1-_alpha_t) * P_prev;
  }

  if (include_phi){
    std::array<double, Triangle::n_dofs_max> phi_prev_block;
    std::array<double, Triangle::n_dofs_max> phi_next_block;
    
    const std::uint32_t* p_dofs = p_dofs_[id];
    for (std::size_t i = 0; i < p_dofs_.stride(); ++i){
        phi_prev_block[i] = phi_prev_data_[p_dofs[i]];
        phi_next_block[i] = phi_next_data_[p_dofs[i]];
    }

    double Phi_prev = std::inner_product(_Np_.data(), _Np_.data()+ncoeffs_p, phi_prev_block.begin(), 0.0);
    double Phi_next = std::inner_product(_Np_.data(), _Np_.data()+ncoeffs_p, phi_next_block.begin(), 0.0);
    fields.Phi = _alpha_t * Phi_next + (1-_alpha_t) * Phi_prev;
  }

  // cell_type
  fields.cell_type = cell_type_[id];
}

void XDMFTriangleInterpol::enable_reflection()
{
  build_facet_neighbours(facet_neigh_, mesh, dolfin_cells_,
                         dolfin2local_.empty() ? nullptr : &dolfin2local_,
                         periodic, x_min, x_max, dim, 1e-12);
  period_ = periodic_lengths(periodic, x_min, x_max, dim);
  can_reflect = true;
}

bool XDMFTriangleInterpol::reflect(const Vector3d& x, Vector3d& dx, CellPos& pos)
{
  return reflect_in_cells(triangles_, facet_neigh_, 3, period_, x, dx, pos,
                          [this](const Vector3d& p){ return _modx(p); });
}

#endif
