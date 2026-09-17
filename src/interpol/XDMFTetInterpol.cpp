#ifdef USE_DOLFIN
#include "geometry.hpp"
#include "loader_params.hpp"
#include "XDMFTetInterpol.hpp"
#include <array>
#include <algorithm>
#include <limits>
#include <map>
#include <numeric>
#include "Timestamps.hpp"
//#include "H5Cpp.h"
#include <cassert>

#include "dolfin_elements/P1_3.h"
//#include "dolfin_elements/P2_3.h"
#include "dolfin_elements/vP1_3.h"
//#include "dolfin_elements/vP2_3.h"
#include "PeriodicBC.hpp"
#include "dolfin_helpers.hpp"
#include "xdmf_helpers.hpp"

//using namespace H5;

XDMFTetInterpol::XDMFTetInterpol(const std::string& infilename)
  : Interpol(infilename)
{

  // Input file (e.g. dolfin_params.dat)
  dolfin_params = partrac::parse_file_or_exit(xdmf_schema("xdmftet"), infilename);

  std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  // Set periodicity
  if (dolfin_params.get<bool>("periodic_x")){
    periodic[0] = true;
  }
  if (dolfin_params.get<bool>("periodic_y")){
    periodic[1] = true;
  }
  if (dolfin_params.get<bool>("periodic_z")){
    periodic[2] = true;
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

  std::string cell_type_str = "tetrahedron";
  Uint gdim = 3;

  std::unique_ptr<dolfin::CellType> cell_type(dolfin::CellType::create(cell_type_str));

  std::vector<std::int64_t> coords_shape = dolfin::HDF5Interface::get_dataset_shape(meshfile.h5_id(), geometry_path);

  meshfile.read(mesh_in, topology_path, geometry_path, gdim, *cell_type, -1, coords_shape[0], false);

  mesh = std::make_shared<dolfin::Mesh>(mesh_in);
  dim = mesh->geometry().dim();
  assert(gdim == dim);
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
  u_space_ = std::make_shared<vP1_3::FunctionSpace>(mesh, constrained_domain);
  //u_space_ = std::make_shared<vP1_3::FunctionSpace>(mesh);
  ncoeffs_u = 4;

  // Pressure
  p_space_ = std::make_shared<P1_3::FunctionSpace>(mesh, constrained_domain);
  //p_space_ = std::make_shared<P1_3::FunctionSpace>(mesh);
  ncoeffs_p = 4;

  // Precompute all tets Taylor-Hood P2-P1
  // FIXME compute on the fly and save
  tets_.resize(mesh->num_cells());
  dolfin_cells_.resize(mesh->num_cells());
  cell2cells_.resize(mesh->num_cells());

  const std::vector<std::uint32_t> order =
    cell_order(*u_space_->dofmap(), mesh->num_cells(), dolfin_params.get<std::string>("renumber_cells"), dolfin2local_);
  for (std::size_t l = 0; l < mesh->num_cells(); ++l)
  {
    dolfin::Cell dolfin_cell(*mesh, order[l]);
    tets_[l] = Tet(dolfin_cell);
    dolfin_cells_[l] = dolfin_cell;
  }
  // Build cell neighbour list for lookup speed
  build_neighbor_list(cell2cells_, mesh, dolfin_cells_,
                      dolfin2local_.empty() ? nullptr : &dolfin2local_);

  // Identify edge cells
  cell_type_.resize(mesh->num_cells());

  double tol = 1e-4; // 1e-12; // heuristic

  apply_periodic_boundaries(cell2cells_, periodic, x_min, x_max, mesh, dolfin_cells_, dim, tol);

  label_cell_type(cell_type_, cell2cells_, dim);

  cell_normal_.resize(mesh->num_cells());
  cell_facet_midpoint_.resize(mesh->num_cells());

  for ( Uint i=0; i < mesh->num_cells(); ++i)
  {
    if (cell_type_[i] == 1)
    {
      auto facets = dolfin_cells_[i].entities(dim-1);
      for ( std::size_t j = 0; j < dolfin_cells_[i].num_entities(dim-1); ++j ){
        dolfin::Facet dolfin_facet(*mesh, facets[j]);

        if (dolfin_facet.exterior()){
          Vector3d pt(dolfin_facet.midpoint().coordinates());

          bool periodic_facet = false;
          for ( Uint k=0; k < static_cast<Uint>(dim); ++k)
          {
            if (periodic[k] && (pt[k] < x_min[k] + tol || pt[k] > x_max[k] - tol))
            {
              periodic_facet = true;
              break;
            }
          }
          if (!periodic_facet){
            Vector3d n_loc (dolfin_facet.normal().coordinates());
            cell_normal_[i] = n_loc;
            cell_facet_midpoint_[i] = pt;
          }
        }
      }
    }
  }

  // P2 near walls: edge (default) or none
  wall_p2_ = dolfin_params.get<std::string>("wall_p2") == "none" ? WallP2::None : WallP2::Edge;

  std::cout << "Built neighbour list" << std::endl;
  
  // const dolfin::GenericDofMap& u_dofmap = *u_space_->dofmap();
  auto xdof = p_space_->tabulate_dof_coordinates();

  std::map<std::tuple<double, double, double>, Uint> xmap;
  for ( Uint j=0; j < xdof.size()/dim; ++j ){
    xmap[{xdof[dim*j], xdof[dim*j+1], xdof[dim*j+2]}] = j;
  }

  std::vector<double> xdata;
  read_dataset_vector(h5filename_u, geometry_path, xdata, dim);

  // Periodic images of a dof's vertex are not in xmap
  const Uint unset = std::numeric_limits<Uint>::max();
  j2i.assign(xdof.size()/dim, unset);
  for ( Uint i=0; i < xdata.size()/dim; ++i ){
    const auto it = xmap.find({xdata[dim*i], xdata[dim*i+1], xdata[dim*i+2]});
    if (it != xmap.end())
      j2i[it->second] = i;
  }
  if (std::find(j2i.begin(), j2i.end(), unset) != j2i.end()){
    std::cout << "XDMFTetInterpol: a dof has no vertex in " << h5filename_u << std::endl;
    exit(1);
  }

  u_prev_data_.resize(xdof.size());
  u_next_data_.resize(xdof.size());
  p_prev_data_.resize(xdof.size()/dim);
  p_next_data_.resize(xdof.size()/dim);
  phi_prev_data_.resize(xdof.size()/dim);
  phi_next_data_.resize(xdof.size()/dim);

  check_dofs_fit(ncoeffs_u, ncoeffs_p, Tet::n_dofs_max, "XDMFTetInterpol");

  std::cout << "Setting max threads: " << omp_get_max_threads() << std::endl;

  found_.resize(omp_get_max_threads());

  // Precomputing dofs

  const dolfin::GenericDofMap& u_dofmap = *u_space_->dofmap();
  const dolfin::GenericDofMap& p_dofmap = *p_space_->dofmap();
  
  u_dofs_.build(u_dofmap, dolfin_cells_, "XDMFTetInterpol");
  u_dofs_.check_stride(3*ncoeffs_u, "XDMFTetInterpol");

  p_dofs_.build(p_dofmap, dolfin_cells_, "XDMFTetInterpol");
  p_dofs_.check_stride(ncoeffs_p, "XDMFTetInterpol");

  if (wall_p2_ == WallP2::Edge)
    build_wall_edges(tol);
}

void XDMFTetInterpol::update(const double t)
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



void XDMFTetInterpol::build_wall_edges(const double tol)
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

double XDMFTetInterpol::rest_tol(const std::vector<double>& u_data) const
{
  // Round-off of the largest velocity
  double u_max = 0.;
  for (const double u : u_data) u_max = std::max(u_max, std::abs(u));
  return 1e-12*u_max;
}

bool XDMFTetInterpol::wall_block(const double* u, double* u2, const WallEdges& w,
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

Vector3d XDMFTetInterpol::_modx(const Vector3d &x){
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


bool XDMFTetInterpol::locate(const Vector3d &x, const double t, CellPos& pos)
{
  assert(t <= t_next && t >= t_prev);
  const Vector3d xx = _modx(x);
  return locate_in_cells(tets_, cell2cells_, *mesh, dim, xx, pos,
                         found_, dolfin2local_.empty() ? nullptr : &dolfin2local_);
}

void XDMFTetInterpol::evaluate(const Vector3d &x, const double tin, const CellPos& pos, PointValues& fields)
{
  const double _alpha_t = stamp_weight(tin, t_prev, t_next);

  // Compute Pk-Pl basis at x
  const int id = pos.id;
  const double r1 = pos.bary[0], r2 = pos.bary[1], r3 = pos.bary[2], r4 = pos.bary[3];

  std::array<double, Tet::n_dofs_max> _Nu_;
  std::array<double, Tet::n_dofs_max> _Np_;
  std::array<double, Tet::n_dofs_max> _Nux_;
  std::array<double, Tet::n_dofs_max> _Nuy_;
  std::array<double, Tet::n_dofs_max> _Nuz_;

  tets_[id].linearbasis(r1, r2, r3, r4, _Nu_.data());

  if (include_pressure){
    tets_[id].linearbasis(r1, r2, r3, r4, _Np_.data());
  }

  std::array<double, Tet::n_dofs_max*3> u_prev_block;
  std::array<double, Tet::n_dofs_max*3> u_next_block;

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
                     std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, &u_prev_block[2*ncoeffs_u], 0.0)};
  Vector3d U_next = {std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, u_next_block.begin(), 0.0),
                     std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, &u_next_block[ncoeffs_u], 0.0),
                     std::inner_product(_Nu_.data(), _Nu_.data()+ncoeffs_u, &u_next_block[2*ncoeffs_u], 0.0)};

  // Gradient
  Matrix3d gradU_prev = Matrix3d::Zero(), gradU_next = Matrix3d::Zero();
  if (wants_gradient()){
    tets_[id].linearderiv(r1, r2, r3, r4, _Nux_.data(), _Nuy_.data(), _Nuz_.data());

    gradU_prev <<
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, u_prev_block.begin(), 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, u_prev_block.begin(), 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, u_prev_block.begin(), 0.0),
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, &u_prev_block[ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, &u_prev_block[ncoeffs_u], 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, &u_prev_block[ncoeffs_u], 0.0),
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, &u_prev_block[2*ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, &u_prev_block[2*ncoeffs_u], 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, &u_prev_block[2*ncoeffs_u], 0.0);
    gradU_next << 
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, u_next_block.begin(), 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, u_next_block.begin(), 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, u_next_block.begin(), 0.0),
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, &u_next_block[ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, &u_next_block[ncoeffs_u], 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, &u_next_block[ncoeffs_u], 0.0),
      std::inner_product(_Nux_.data(), _Nux_.data()+ncoeffs_u, &u_next_block[2*ncoeffs_u], 0.0),
      std::inner_product(_Nuy_.data(), _Nuy_.data()+ncoeffs_u, &u_next_block[2*ncoeffs_u], 0.0),
      std::inner_product(_Nuz_.data(), _Nuz_.data()+ncoeffs_u, &u_next_block[2*ncoeffs_u], 0.0);
  }


  // Quadratic near walls
  std::array<double, 30> u_prev_block_2;
  std::array<double, 30> u_next_block_2;
  bool quad_prev = false;
  bool quad_next = false;
  if (wall_p2_ == WallP2::Edge && wall_index_[id] >= 0){
    const WallEdges& w = wall_cells_[wall_index_[id]];
    quad_prev = wall_block(u_prev_block.data(), u_prev_block_2.data(), w, rest_tol_prev_);
    quad_next = wall_block(u_next_block.data(), u_next_block_2.data(), w, rest_tol_next_);
  }

  if (quad_prev || quad_next){
    const Uint n2 = 10;
    std::array<double, 10> _Nu2_;
    std::array<double, 10> _Nu2x_{};
    std::array<double, 10> _Nu2y_{};
    std::array<double, 10> _Nu2z_{};
    tets_[id].quadbasis(r1, r2, r3, r4, _Nu2_.data());
    if (wants_gradient())
      tets_[id].quadderiv(r1, r2, r3, r4, _Nu2x_.data(), _Nu2y_.data(), _Nu2z_.data());

    const auto quad = [&](const std::array<double, 30>& b, Vector3d& U, Matrix3d& gradU){
      for (int c = 0; c < 3; ++c)
        U[c] = std::inner_product(_Nu2_.begin(), _Nu2_.end(), &b[n2*c], 0.0);
      if (wants_gradient()){
        for (int c = 0; c < 3; ++c){
          gradU(c, 0) = std::inner_product(_Nu2x_.begin(), _Nu2x_.end(), &b[n2*c], 0.0);
          gradU(c, 1) = std::inner_product(_Nu2y_.begin(), _Nu2y_.end(), &b[n2*c], 0.0);
          gradU(c, 2) = std::inner_product(_Nu2z_.begin(), _Nu2z_.end(), &b[n2*c], 0.0);
        }
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
    std::array<double, Tet::n_dofs_max> p_prev_block;
    std::array<double, Tet::n_dofs_max> p_next_block;
    
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
    std::array<double, Tet::n_dofs_max> phi_prev_block;
    std::array<double, Tet::n_dofs_max> phi_next_block;
    
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

void XDMFTetInterpol::enable_reflection()
{
  build_facet_neighbours(facet_neigh_, mesh, dolfin_cells_,
                         dolfin2local_.empty() ? nullptr : &dolfin2local_,
                         periodic, x_min, x_max, dim, 1e-4);
  period_ = periodic_lengths(periodic, x_min, x_max, dim);
  can_reflect = true;
}

bool XDMFTetInterpol::reflect(const Vector3d& x, Vector3d& dx, CellPos& pos)
{
  return reflect_in_cells(tets_, facet_neigh_, 4, period_, x, dx, pos,
                          [this](const Vector3d& p){ return _modx(p); });
}

#endif
