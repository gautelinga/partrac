#include "TriangleFreqInterpol.hpp"
#include "geometry.hpp"
#include "loader_params.hpp"
#include "mesh_tables.hpp"
#include "p12_eval.hpp"
#include "phase_timing.hpp"
#include "simplex_load.hpp"
#include <array>
#include <cassert>
#include <numeric>

namespace {

// The element a parameter file declares; the file's own signature decides the
// element, but a name this loader does not know is still a mistake
int declared_degree(const std::string& space, const char* what){
  if (space == "P1") return 1;
  if (space == "P2") return 2;
  partrac::fail("unrecognized ", what, " element: ", space);
  return 0;
}

// The components are one function space written many times: a later file that
// holds another element is not summable with the first
void same_element(const simplex_load::Element& a, const simplex_load::Element& b,
                  const std::string& path, const std::string& field){
  if (a.degree != b.degree || a.ncomp != b.ncomp)
    partrac::fail(path, ": '", field, "' is of degree ", b.degree, " with ", b.ncomp,
                  " components, the first component's degree ", a.degree, " with ", a.ncomp);
}

}  // namespace

TriangleFreqInterpol::TriangleFreqInterpol(const std::string& infilename)
  : MeshCore<Triangle>(infilename)
{
  constexpr int nv = Triangle::n_verts;
  constexpr int ne = nv*(nv-1)/2;
  constexpr int D = nv - 1;
  partrac::phase_begin("load");
  dolfin_params = partrac::parse_file_or_exit(triangle_freq_schema(), infilename);

  // base frequency
  const double tau = dolfin_params.get<double>("tau");
  if (tau > 0)
    omega0 = 2 * M_PI / tau;

  const std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  // Using the FreqStamps class to hold frequency data
  fs.initialize(get_folder() + "/" + dolfin_params.get<std::string>("freqstamps"));
  if (fs.size() < 1)
    partrac::fail(infilename, ": the frequency stamps list no component");

  read_mesh_params();
  const std::string u_field = dolfin_params.get<std::string>("velocity_field");
  const std::string p_field = dolfin_params.get<std::string>("pressure_field");
  const int want_u = declared_degree(dolfin_params.get<std::string>("velocity_space"), "velocity");
  const int want_p = include_pressure
    ? declared_degree(dolfin_params.get<std::string>("pressure_space"), "pressure") : 0;
  partrac::phase("params");

  const std::string mesh_file = get_folder() + "/" + dolfin_params.get<std::string>("mesh");
  simplex_load::MeshData m;
  std::vector<std::uint32_t> cell_perm;
  simplex_load::read_mesh(mesh_file, nv, m, cell_perm);
  cell_perm.clear();
  cell_perm.shrink_to_fit();
  dim = m.gdim;
  x_min = m.x_min;
  x_max = m.x_max;
  ncells_ = m.ncells;
  nverts_ = m.nverts;
  set_period();

  // The element comes from the first component; the parameter file must not claim another
  const std::string first = get_folder() + "/" + fs.get(0).filename;
  const simplex_load::Element el_u = simplex_load::read_element(first, u_field, nv);
  simplex_load::Element el_p;
  if (include_pressure)
    el_p = simplex_load::read_element(first, p_field, nv);
  if (want_u != el_u.degree)
    partrac::fail(infilename, ": velocity_space is P", want_u, ", but '", u_field,
                  "' in ", first, " is of degree ", el_u.degree);
  if (include_pressure && want_p != el_p.degree)
    partrac::fail(infilename, ": pressure_space is P", want_p, ", but '", p_field,
                  "' in ", first, " is of degree ", el_p.degree);
  if (el_u.ncomp != std::size_t(D))
    partrac::fail(first, ": '", u_field, "' has ", el_u.ncomp, " components, not ", D);
  if (include_pressure && el_p.ncomp != 1)
    partrac::fail(first, ": '", p_field, "' has ", el_p.ncomp, " components, not one");
  ncoeffs_u = Uint(simplex_load::nodes_per_cell(nv, el_u.degree));
  ncoeffs_p = include_pressure ? Uint(simplex_load::nodes_per_cell(nv, el_p.degree)) : 0;
  check_dofs_fit(ncoeffs_u, ncoeffs_p, Triangle::n_dofs_max, "TriangleFreqInterpol");

  // The edges, when either field carries midside nodes
  std::vector<std::uint32_t> edges;
  std::size_t nedges = 0;
  const bool quadratic = el_u.degree == 2 || (include_pressure && el_p.degree == 2);
  if (quadratic)
    nedges = mesh_tables::build_edge_table<nv>(m.topo, m.ncells, edges);
  partrac::phase("edge table");

  // A node and its periodic images are one node, the master's dof serving all
  simplex_load::NodePairs np;
  np.build(m, edges, nedges, nv, periodic, x_min, x_max, periodic_tol);
  const std::size_t n_u = m.nverts + (el_u.degree == 2 ? nedges : 0);
  const std::vector<std::uint32_t> nodes_u = np.masters(std::vector<std::uint32_t>(), n_u);
  u_dofs_.fill(m.topo, edges, m.ncells, nv, ne, el_u.degree == 2, m.nverts,
               nodes_u.empty() ? nullptr : nodes_u.data());
  u_dofs_.check_stride(ncoeffs_u, "TriangleFreqInterpol");
  if (include_pressure){
    const std::size_t n_p = m.nverts + (el_p.degree == 2 ? nedges : 0);
    const std::vector<std::uint32_t> nodes_p = np.masters(std::vector<std::uint32_t>(), n_p);
    p_dofs_.fill(m.topo, edges, m.ncells, nv, ne, el_p.degree == 2, m.nverts,
                 nodes_p.empty() ? nullptr : nodes_p.data());
    p_dofs_.check_stride(ncoeffs_p, "TriangleFreqInterpol");
  }
  partrac::phase("cell dofs");

  // Every component through one dof -> nodes mapping: the components are one
  // function space, so only the first has to be read with its dof table
  u_coefficients_.resize(std::size_t(fs.size()));
  if (include_pressure)
    p_coefficients_.resize(std::size_t(fs.size()));
  mesh_tables::DofNodes u_map, p_map;
  std::vector<double> u_values, p_values;
  const std::vector<std::uint32_t> no_edges;
  for (std::size_t iFreq = 0; iFreq < std::size_t(fs.size()); ++iFreq){
    const std::string path = get_folder() + "/" + fs.get(int(iFreq)).filename;
    if (iFreq == 0){
      std::vector<std::uint32_t> rows;
      std::vector<double> vec;
      simplex_load::read_field(path, u_field, m, el_u, nv, rows, vec);
      mesh_tables::scatter_dofs_to_nodes<nv>(m.topo, el_u.degree == 2 ? edges : no_edges, m.ncells,
                                             m.nverts, el_u.degree == 2 ? nedges : 0, el_u.ncomp,
                                             rows, vec, u_values, u_map);
      np.check(u_values, el_u.ncomp, path, u_field);
      if (include_pressure){
        simplex_load::read_field(path, p_field, m, el_p, nv, rows, vec);
        mesh_tables::scatter_dofs_to_nodes<nv>(m.topo, el_p.degree == 2 ? edges : no_edges, m.ncells,
                                               m.nverts, el_p.degree == 2 ? nedges : 0, el_p.ncomp,
                                               rows, vec, p_values, p_map);
        np.check(p_values, el_p.ncomp, path, p_field);
      }
    }
    else {
      same_element(el_u, simplex_load::read_element(path, u_field, nv), path, u_field);
      u_values.assign(u_values.size(), 0.);
      simplex_load::read_vector(path, u_field, u_map, u_values);
      if (include_pressure){
        same_element(el_p, simplex_load::read_element(path, p_field, nv), path, p_field);
        p_values.assign(p_values.size(), 0.);
        simplex_load::read_vector(path, p_field, p_map, p_values);
      }
    }
    // The cells' coefficients, the velocity's components blocked
    std::vector<double>& uc = u_coefficients_[iFreq];
    uc.assign(m.ncells*std::size_t(D)*ncoeffs_u, 0.);
#pragma omp parallel for schedule(static)
    for (std::ptrdiff_t i = 0; i < std::ptrdiff_t(m.ncells); ++i){
      const std::uint32_t* row = u_dofs_[std::size_t(i)];
      double* dst = uc.data() + std::size_t(i)*std::size_t(D)*ncoeffs_u;
      for (Uint k = 0; k < ncoeffs_u; ++k)
        for (int c = 0; c < D; ++c)
          dst[std::size_t(c)*ncoeffs_u + k] = u_values[std::size_t(row[k])*std::size_t(D) + std::size_t(c)];
    }
    if (include_pressure){
      std::vector<double>& pc = p_coefficients_[iFreq];
      pc.assign(m.ncells*ncoeffs_p, 0.);
#pragma omp parallel for schedule(static)
      for (std::ptrdiff_t i = 0; i < std::ptrdiff_t(m.ncells); ++i){
        const std::uint32_t* row = p_dofs_[std::size_t(i)];
        double* dst = pc.data() + std::size_t(i)*ncoeffs_p;
        for (Uint k = 0; k < ncoeffs_p; ++k) dst[k] = p_values[row[k]];
      }
    }
  }
  edges.clear();
  edges.shrink_to_fit();
  partrac::phase("component coefficients");

  mesh_tables::build_facet_neighbours<nv>(m.topo, m.ncells, m.coords, dim, periodic,
                                          x_min, x_max, periodic_tol, facet_neigh_);
  partrac::phase("facet table");

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
      cells_[i] = Triangle(c + std::size_t(row[0])*g, c + std::size_t(row[1])*g,
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

bool TriangleFreqInterpol::locate_tree(const Vector3d& xx, CellPos& pos)
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

void TriangleFreqInterpol::update(const double t)
{
  if (!is_initialized){
    is_initialized = true;
  }
  t_update = t;
}

void TriangleFreqInterpol::evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  const int id = pos.id;
  // Assume found in fluid domain

  // update frequency weights
  static thread_local std::vector<double> w_f_; w_f_.resize(fs.size());
  static thread_local std::vector<double> wt_f_; wt_f_.resize(fs.size());
  for (std::size_t iFreq=0; iFreq < static_cast<std::size_t>(fs.size()); ++iFreq){
    FreqStamp& f = fs.get(iFreq);
    double a = f.a;
    double t_shift = f.t;
    w_f_[iFreq] = a * cos(omega0 * (iFreq * t + t_shift));
    wt_f_[iFreq] = - a * omega0 * iFreq * sin(omega0 * (iFreq * t + t_shift));
  }

  // Compute Pk-Pl basis at x

  std::array<double, Triangle::n_dofs_max> Nu_;
  std::array<double, Triangle::n_dofs_max> Np_;
  std::array<double, Triangle::n_dofs_max> Nux_;
  std::array<double, Triangle::n_dofs_max> Nuy_;

  cell_basis(cells_[id], pos.bary, ncoeffs_u, Nu_.data(), "u");
  if (include_pressure)
    cell_basis(cells_[id], pos.bary, ncoeffs_p, Np_.data(), "p");

  const std::size_t u_row = std::size_t(id)*std::size_t(Triangle::n_verts - 1)*ncoeffs_u;
  const std::size_t p_row = std::size_t(id)*ncoeffs_p;

  static thread_local std::vector<double> ux_f_; ux_f_.resize(fs.size());
  static thread_local std::vector<double> uy_f_; uy_f_.resize(fs.size());
  static thread_local std::vector<double> p_f_; p_f_.resize(fs.size());

  for (std::size_t iFreq=0; iFreq < static_cast<std::size_t>(fs.size()); ++iFreq){
    const double* u_c = u_coefficients_[iFreq].data() + u_row;
    ux_f_[iFreq] = std::inner_product(Nu_.data(), Nu_.data()+ncoeffs_u, u_c, 0.0);
    uy_f_[iFreq] = std::inner_product(Nu_.data(), Nu_.data()+ncoeffs_u, u_c + ncoeffs_u, 0.0);
    if (include_pressure)
      p_f_[iFreq] = std::inner_product(Np_.data(), Np_.data()+ncoeffs_p,
                                       p_coefficients_[iFreq].data() + p_row, 0.0);
  }

  // Update
  fields.U = { std::inner_product(w_f_.begin(), w_f_.end(), ux_f_.begin(), 0.0),
               std::inner_product(w_f_.begin(), w_f_.end(), uy_f_.begin(), 0.0),
               0.};
  fields.A = { std::inner_product(wt_f_.begin(), wt_f_.end(), ux_f_.begin(), 0.0),
               std::inner_product(wt_f_.begin(), wt_f_.end(), uy_f_.begin(), 0.0),
               0.};
  if (include_pressure){
    fields.P = std::inner_product(w_f_.begin(), w_f_.end(), p_f_.begin(), 0.0);
  }

  if (wants_gradient()){
    cell_deriv(cells_[id], pos.bary, ncoeffs_u, Nux_.data(), Nuy_.data(), nullptr, "u");

    static thread_local std::vector<double> uxx_f_; uxx_f_.resize(fs.size());
    static thread_local std::vector<double> uxy_f_; uxy_f_.resize(fs.size());
    static thread_local std::vector<double> uyx_f_; uyx_f_.resize(fs.size());
    static thread_local std::vector<double> uyy_f_; uyy_f_.resize(fs.size());

    for (std::size_t iFreq=0; iFreq < static_cast<std::size_t>(fs.size()); ++iFreq){
      const double* u_c = u_coefficients_[iFreq].data() + u_row;
      uxx_f_[iFreq] = std::inner_product(Nux_.data(), Nux_.data()+ncoeffs_u, u_c, 0.0);
      uxy_f_[iFreq] = std::inner_product(Nuy_.data(), Nuy_.data()+ncoeffs_u, u_c, 0.0);
      uyx_f_[iFreq] = std::inner_product(Nux_.data(), Nux_.data()+ncoeffs_u, u_c + ncoeffs_u, 0.0);
      uyy_f_[iFreq] = std::inner_product(Nuy_.data(), Nuy_.data()+ncoeffs_u, u_c + ncoeffs_u, 0.0);
    }

    fields.gradU(0, 0) = std::inner_product(w_f_.begin(), w_f_.end(), uxx_f_.begin(), 0.0);
    fields.gradU(0, 1) = std::inner_product(w_f_.begin(), w_f_.end(), uxy_f_.begin(), 0.0);
    fields.gradU(1, 0) = std::inner_product(w_f_.begin(), w_f_.end(), uyx_f_.begin(), 0.0);
    fields.gradU(1, 1) = std::inner_product(w_f_.begin(), w_f_.end(), uyy_f_.begin(), 0.0);

    fields.gradA(0, 0) = std::inner_product(wt_f_.begin(), wt_f_.end(), uxx_f_.begin(), 0.0);
    fields.gradA(0, 1) = std::inner_product(wt_f_.begin(), wt_f_.end(), uxy_f_.begin(), 0.0);
    fields.gradA(1, 0) = std::inner_product(wt_f_.begin(), wt_f_.end(), uyx_f_.begin(), 0.0);
    fields.gradA(1, 1) = std::inner_product(wt_f_.begin(), wt_f_.end(), uyy_f_.begin(), 0.0);
  }
}
