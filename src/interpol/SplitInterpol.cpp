#include "SplitInterpol.hpp"
#include "Error.hpp"
#include "loader_params.hpp"
#include "mesh_tables.hpp"
#include "near_wall.hpp"
#include "phase_timing.hpp"
#include "simplex_load.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <omp.h>

namespace {

// Net flux accepted, relative to a cell's largest facet flux floored at the stamp's
constexpr double flux_tol = 1e-9;
constexpr double flux_floor = 1e-6;

// A cell whose interior values exceed its boundary values by this much is
// warned about: the split field returns what the cell's shape forces
constexpr double sliver_warn = 50.;

// A cell's largest interior node value against its largest boundary one
inline double amplification(const double* g, const std::size_t nb,
                            const double* u_int, const std::size_t ni){
  double a = 0., b = 0.;
  for (std::size_t k = 0; k < nb; ++k) b = std::max(b, std::abs(g[k]));
  for (std::size_t k = 0; k < ni; ++k) a = std::max(a, std::abs(u_int[k]));
  return b > 0. ? a/b : 0.;
}

// The cell's affine Jacobian and its inverse, row-major; the rows of the
// inverse are the gradients of the barycentrics after the first
template<typename Cell>
inline void cell_jacobian(const Cell& cell, double* J, double* Jinv){
  constexpr int D = Cell::n_verts - 1;
  Eigen::Matrix<double, D, D> Ji;
  for (int k = 0; k < D; ++k){
    const Vector3d gl = cell.bary_grad(k + 1);
    for (int q = 0; q < D; ++q){
      Ji(k, q) = gl[q];
      Jinv[k*D + q] = gl[q];
    }
  }
  const Eigen::Matrix<double, D, D> Jm = Ji.inverse();
  for (int q = 0; q < D; ++q)
    for (int k = 0; k < D; ++k) J[q*D + k] = Jm(q, k);
}

// D times the cell's volume
template<int D>
inline double cell_measure(const double* J){
  if constexpr (D == 2) return std::abs(J[0]*J[3] - J[1]*J[2]);
  else return std::abs(J[0]*(J[4]*J[8] - J[5]*J[7]) - J[1]*(J[3]*J[8] - J[5]*J[6])
                       + J[2]*(J[3]*J[7] - J[4]*J[6]))/2.;
}

// The cell's shape: the singular value ratio of its affine Jacobian
template<int D>
double cond_J(const double* J){
  Eigen::Matrix<double, D, D> M;
  for (int i = 0; i < D; ++i)
    for (int j = 0; j < D; ++j) M(i, j) = J[i*D + j];
  const auto sv = Eigen::JacobiSVD<Eigen::Matrix<double, D, D>>(M).singularValues();
  return sv(D - 1) > 0. ? sv(0)/sv(D - 1) : std::numeric_limits<double>::infinity();
}

}  // namespace

template<typename Cell>
SplitInterpol<Cell>::SplitInterpol(const std::string& infilename)
  : MeshCore<Cell>(infilename)
{
  constexpr int nv = Cell::n_verts;
  partrac::phase_begin("load");
  dolfin_params = partrac::parse_file_or_exit(dolfin_h5_schema(mode), infilename);

  const std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  read_mesh_params();
  include_phi = dolfin_params.template get<bool>("include_phi");
  if (dolfin_params.template get<std::string>("wall_p2") == "edge")
    partrac::fail(infilename, ": wall_p2 = edge is a rule for a P1 velocity at a wall; "
                  "divfree = true has no slip from div u = 0 and replaces it");
  if (dolfin_params.template get<bool>("mesh_cache"))
    partrac::fail(infilename, ": mesh_cache is not read with divfree = true");

  ts_.initialize(get_folder() + "/" + dolfin_params.template get<std::string>("timestamps"));
  partrac::phase("params");

  simplex_load::Request req;
  req.infilename = infilename;
  // The first stamp, by the path update builds for every other one
  req.field_file = get_folder() + "/" + ts_.get(ts_.get_t_min()).prev.filename;
  req.what = "SplitInterpol";
  req.nv = nv;
  req.include_pressure = include_pressure;
  req.include_phi = include_phi;
  req.n_dofs_max = Cell::n_dofs_max;
  req.periodic = periodic;
  req.periodic_tol = periodic_tol;
  simplex_load::request_from_params(req, dolfin_params, get_folder());
  u_field = req.u_field;
  p_field = req.p_field;
  phi_field = req.phi_field;
  simplex_load::Tables t;
  simplex_load::build_tables(req, t);

  this->adopt_tables(t);
  ncoeffs_phi = t.ncoeffs_phi;
  phi_dofs_ = std::move(t.phi_dofs);
  ncells_ = t.mesh.ncells;
  nverts_ = t.mesh.nverts;

  if (ncoeffs_u != Uint(Cell::n_dofs_max))
    partrac::fail(infilename, ": divfree = true needs a P2 velocity, and this one has ",
                  ncoeffs_u, " dofs a cell; prepare the file with python/divfree/divfree_clean.py");

  auto& a = stamps_.a();
  simplex_load::read_field_by_node(req.field_file, u_field, t, t.el_u, t.node_order(t.el_u),
                                   a.u, u_map_);
  if (include_pressure)
    simplex_load::read_field_by_node(req.field_file, p_field, t, t.el_p, t.node_order(t.el_p),
                                     a.p, p_map_);
  if (include_phi)
    simplex_load::read_field_by_node(req.field_file, phi_field, t, t.el_phi,
                                     t.node_order(t.el_phi), a.phi, phi_map_);
  stamps_.hold_a(req.field_file);

  topo_ = std::move(t.mesh.topo);
  coords_ = std::move(t.mesh.coords);

  // Identify edge cells
  near_wall::label_cell_type(cell_type_, facet_neigh_, nv);
  std::cout << "Built neighbour list" << std::endl;

  // The per-cell geometry and the tree, from the topology
  mesh_tables::build_cells_and_tree<Cell>(topo_, coords_, ncells_, nverts_, dim, cells_, tree_);

  int rank = 0;
  R_ = split_eval::reference_matrix<Cell>(rank);
  if (rank != int(Split::n_int)*D)
    partrac::fail("SplitInterpol: the split's interior field is not unique: ", rank,
                  " of ", Split::n_int*std::size_t(D), " dofs pinned");
  partrac::phase("reference matrix");

  build_interior(req.field_file, a);
  partrac::phase("interior values");
  partrac::phase_total();

  std::cout << "Setting max threads: " << omp_get_max_threads() << std::endl;
}

template<typename Cell>
void SplitInterpol<Cell>::build_interior(const std::string& file, Stamp& s)
{
  constexpr std::size_t nb = Split::n_bnd*std::size_t(D);
  constexpr std::size_t ni = Split::n_int*std::size_t(D);
  s.u_int.resize(ncells_*ni);
  const double* u = s.u.data();
  double* out = s.u_int.data();
  double facet_max = 0., over = 0., loudest = 0.;
  bool nonfinite = false;
  // A throw here would terminate, so the refusal waits for the loop to end
#pragma omp parallel for schedule(static) \
  reduction(max:facet_max) reduction(max:over) reduction(max:loudest) \
  reduction(||:nonfinite)
  for (std::ptrdiff_t i = 0; i < std::ptrdiff_t(ncells_); ++i){
    const std::size_t c = std::size_t(i);
    double g[nb];
    split_eval::gather_bnd<Cell, D>(u_dofs_[c], u, g);
    const Cell& cell = cells_[c];
    double J[D*D], Jinv[D*D];
    cell_jacobian(cell, J, Jinv);
    double net = 0., big = 0.;
    split_eval::flux_parts<Cell>(cell, g, cell_measure<D>(J), net, big);
    nonfinite = nonfinite || !std::isfinite(net) || !std::isfinite(big);
    facet_max = std::max(facet_max, big);
    // The largest net flux its own scale refuses; the floor is applied after the loop
    if (!(std::abs(net) <= flux_tol*big)) over = std::max(over, std::abs(net));
    double* ui = out + c*ni;
    split_eval::interior_values<Cell>(R_, J, Jinv, g, ui);
    loudest = std::max(loudest, amplification(g, nb, ui, ni));
  }
  if (nonfinite){
    for (std::size_t c = 0; c < ncells_; ++c){
      double g[nb], J[D*D], Jinv[D*D], net = 0., big = 0.;
      split_eval::gather_bnd<Cell, D>(u_dofs_[c], u, g);
      cell_jacobian(cells_[c], J, Jinv);
      split_eval::flux_parts<Cell>(cells_[c], g, cell_measure<D>(J), net, big);
      if (!std::isfinite(net) || !std::isfinite(big))
        partrac::fail(file, ": cell ", c, " has a non-finite net or facet flux; "
                      "the velocity holds a NaN or an infinity");
    }
  }
  if (!(over <= flux_tol*flux_floor*facet_max)){
    const double floor = flux_floor*facet_max;
    std::size_t bad = 0;
    double worst = 0.;
    for (std::size_t c = 0; c < ncells_; ++c){
      double g[nb], J[D*D], Jinv[D*D], net = 0., big = 0.;
      split_eval::gather_bnd<Cell, D>(u_dofs_[c], u, g);
      cell_jacobian(cells_[c], J, Jinv);
      split_eval::flux_parts<Cell>(cells_[c], g, cell_measure<D>(J), net, big);
      const double r = std::abs(net)/std::max(big, floor);
      if (r > worst){ worst = r; bad = c; }
    }
    partrac::fail(file, ": cell ", bad, " has a net flux ", worst,
                  " of its scale, above ", flux_tol,
                  "; divfree = true reads a file prepared by python/divfree/divfree_clean.py");
  }
  // One line a stamp: a sliver returns interior values in proportion to cond(J)
  if (loudest > sliver_warn){
    std::size_t bad = 0;
    double cond = 0.;
    for (std::size_t c = 0; c < ncells_; ++c){
      double g[nb], J[D*D], Jinv[D*D], ui[ni];
      split_eval::gather_bnd<Cell, D>(u_dofs_[c], u, g);
      cell_jacobian(cells_[c], J, Jinv);
      split_eval::interior_values<Cell>(R_, J, Jinv, g, ui);
      if (amplification(g, nb, ui, ni) >= loudest){ bad = c; cond = cond_J<D>(J); break; }
    }
    std::cerr << "Warning: " << file << ": cell " << bad << " of aspect cond(J) = " << cond
              << " returns interior values " << loudest
              << " times its largest boundary value, above " << sliver_warn << std::endl;
  }
}

template<typename Cell>
void SplitInterpol<Cell>::read_stamp(const std::string& file, Stamp& s)
{
  // The node counts every stamp shares
  const auto& a = stamps_.a();
  simplex_load::read_into(file, u_field, u_map_, a.u, s.u);
  if (include_pressure) simplex_load::read_into(file, p_field, p_map_, a.p, s.p);
  if (include_phi)      simplex_load::read_into(file, phi_field, phi_map_, a.phi, s.phi);
  build_interior(file, s);
}

template<typename Cell>
void SplitInterpol<Cell>::update(const double t)
{
  const auto sp = ts_.get(t);
  const std::string prev_file = get_folder() + "/" + sp.prev.filename;
  const std::string next_file = get_folder() + "/" + sp.next.filename;

  if (partrac::stamp_reload(is_initialized, t, ts_.get_t_max(),
                            {t_prev, t_next}, {sp.prev.t, sp.next.t})){
    partrac::phase_begin("update");
    const auto fill = stamps_.load(prev_file, next_file,
                                   [&](const std::string& key, Stamp& s){ read_stamp(key, s); });
    std::cout << "Prev: Timestep = " << sp.prev.t << ", "
              << partrac::stamp_note(fill.first, sp.prev.filename) << std::endl;
    std::cout << "Next: Timestep = " << sp.next.t << ", "
              << partrac::stamp_note(fill.second, sp.next.filename) << std::endl;

    u_prev_ = stamps_.prev().u.data();
    u_next_ = stamps_.next().u.data();
    int_prev_ = stamps_.prev().u_int.data();
    int_next_ = stamps_.next().u_int.data();
    p_prev_ = stamps_.prev().p.data();
    p_next_ = stamps_.next().p.data();
    phi_prev_ = stamps_.prev().phi.data();
    phi_next_ = stamps_.next().phi.data();

    partrac::phase("vector read");
    partrac::phase_total();

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  t_update = t;
}

template class SplitInterpol<Triangle>;
template class SplitInterpol<Tet>;
