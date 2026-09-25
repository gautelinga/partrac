// The pieces the mesh loaders are built from, on dolfin's unit meshes built in
// memory: the P1/P2 bases of Triangle and Tet, the evaluation in p12_eval.hpp
// against dolfin's own interpolation, locate, the reflecting walk, and the
// periodic tables, the wall normal. These are what a change to locate or to
// the dof gather touches, and the app tests reach them only through whole runs.
#ifdef USE_DOLFIN
#include <catch2/catch.hpp>

#include <array>
#include <memory>
#include <random>
#include <vector>
#include <dolfin.h>
#include <omp.h>

#include "typedefs.hpp"
#include "Triangle.hpp"
#include "Tet.hpp"
#include "dolfin_ref.hpp"
#include "taylor_hood.hpp"
#include "p12_eval.hpp"
#include "MeshCore.hpp"
#include "StructuredInterpol.hpp"

namespace {

// The walk, then the tree, as MeshCore::locate
template<typename Cell>
bool locate_in_cells(const std::vector<Cell>& cells, const std::vector<std::int32_t>& across,
                     const dolfin::Mesh& mesh, const Uint dim, const Vector3d& xx, CellPos& pos,
                     const std::vector<std::uint32_t>* dolfin2local, FoundCounts* count){
  return walk_to_cell(cells, across, xx, pos, count)
      || tree_to_cell(cells, mesh, dim, xx, pos, dolfin2local, count);
}

template<typename Cell> constexpr int dim_of = Cell::n_verts - 1;

// Unit square or cube with n cells a side, optionally sheared so no cell is
// axis-aligned (an affine map keeps the cells simplices)
template<typename Cell>
std::shared_ptr<dolfin::Mesh> unit_mesh(const std::size_t n, const bool sheared){
  std::shared_ptr<dolfin::Mesh> mesh;
  if constexpr (dim_of<Cell> == 2) mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitSquareMesh(n, n));
  else                             mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitCubeMesh(n, n, n));
  if (sheared){
    std::vector<double>& c = mesh->coordinates();
    const std::size_t d = dim_of<Cell>;
    for (std::size_t i = 0; i < c.size(); i += d){
      const double x = c[i], y = c[i+1];
      c[i] = x + 0.3*y;
      c[i+1] = 0.8*y + 0.1*x;
      if (d == 3) c[i+2] += 0.2*x;
    }
  }
  mesh->init();
  mesh->bounding_box_tree();
  return mesh;
}

// The loaders' per-cell tables, in dolfin's cell order
template<typename Cell>
struct Cells {
  std::shared_ptr<dolfin::Mesh> mesh;
  std::vector<Cell> cells;
  std::vector<dolfin::Cell> dolfin_cells;
  std::vector<std::int32_t> across;
  FoundCounts found;
  explicit Cells(std::shared_ptr<dolfin::Mesh> m) : mesh(m) {
    for (dolfin::CellIterator c(*mesh); !c.end(); ++c){
      dolfin_cells.push_back(*c);
      cells.push_back(Cell(*c));
    }
    const Vector3d lo(0., 0., 0.), hi(1., 1., dim_of<Cell> == 3 ? 1. : 0.);
    build_facet_neighbours(across, mesh, dolfin_cells, nullptr, {false, false, false}, lo, hi, dim_of<Cell>, 1e-12);
  }
  bool locate(const Vector3d& x, CellPos& pos){
    return locate_in_cells(cells, across, *mesh, dim_of<Cell>, x, pos, nullptr, &found);
  }
  // Vertex k of cell id, in the order the cell's barycentrics use
  Vector3d vertex(const std::size_t id, const std::size_t k) const {
    for (dolfin::VertexIterator v(dolfin_cells[id]); !v.end(); ++v)
      if (v.pos() == k){
        Vector3d x = Vector3d::Zero();
        for (int d = 0; d < dim_of<Cell>; ++d) x[d] = v->x(d);
        return x;
      }
    return Vector3d::Zero();
  }
};

// A vector field, linear or quadratic in every component
struct Polynomial {
  int dim;
  bool quadratic;
  double a(const int i) const { return 0.2 - 0.1*i; }
  double A(const int i, const int j) const { return 0.3*(i + 1) - 0.2*j; }
  double B(const int i, const int j, const int k) const { return 0.1*(i + 1)*(j + 1) - 0.05*k; }
  double value(const int i, const Vector3d& x) const {
    double v = a(i);
    for (int j = 0; j < dim; ++j){
      v += A(i, j)*x[j];
      if (quadratic)
        for (int k = 0; k < dim; ++k) v += B(i, j, k)*x[j]*x[k];
    }
    return v;
  }
  // d value(i) / dx_m
  double derivative(const int i, const int m, const Vector3d& x) const {
    double v = A(i, m);
    if (quadratic)
      for (int k = 0; k < dim; ++k) v += (B(i, m, k) + B(i, k, m))*x[k];
    return v;
  }
};

struct PolynomialExpression : public dolfin::Expression {
  Polynomial p;
  explicit PolynomialExpression(const Polynomial& p) : dolfin::Expression(p.dim), p(p) {}
  void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override {
    Vector3d y = Vector3d::Zero();
    for (int d = 0; d < p.dim; ++d) y[d] = x[d];
    for (int i = 0; i < p.dim; ++i) values[i] = p.value(i, y);
  }
};

// A point strictly inside cell id, from barycentrics
template<typename Cell>
Vector3d inside_point(const Cells<Cell>& c, const std::size_t id, std::mt19937& gen){
  std::uniform_real_distribution<> u(0.1, 1.);
  std::array<double, 4> w{};
  double sum = 0.;
  for (int k = 0; k < Cell::n_verts; ++k) sum += (w[k] = u(gen));
  Vector3d x = Vector3d::Zero();
  for (int k = 0; k < Cell::n_verts; ++k) x += (w[k]/sum)*c.vertex(id, k);
  return x;
}

// The midpoints mid_ names: edges 01, 02, 12 of a triangle; 01, 02, 03, 12, 13, 23 of a tet
template<typename Cell>
std::vector<std::array<int, 2>> edges(){
  if constexpr (dim_of<Cell> == 2) return {{0, 1}, {0, 2}, {1, 2}};
  else return {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
}

template<typename Cell>
void check_nodal_bases(){
  Cells<Cell> c(unit_mesh<Cell>(3, true));
  std::mt19937 gen(1);
  for (const std::size_t id : {std::size_t(0), c.cells.size()/2, c.cells.size() - 1}){
    const Cell& cell = c.cells[id];
    std::array<double, 4> bary;
    std::array<double, Cell::n_dofs_max> N;
    // Each basis function is 1 at its own node and 0 at the others
    for (int k = 0; k < Cell::n_verts; ++k){
      cell.contains(c.vertex(id, k), bary);
      cell_basis(cell, bary, Cell::n_verts, N.data(), "test");
      for (int j = 0; j < Cell::n_verts; ++j)
        REQUIRE(N[j] == Approx(j == k ? 1. : 0.).margin(1e-12));
      cell_basis(cell, bary, Cell::n_dofs_max, N.data(), "test");
      for (std::size_t j = 0; j < Cell::n_dofs_max; ++j)
        REQUIRE(N[j] == Approx(int(j) == k ? 1. : 0.).margin(1e-12));
    }
    const auto e = edges<Cell>();
    for (std::size_t m = 0; m < e.size(); ++m){
      cell.contains(0.5*(c.vertex(id, e[m][0]) + c.vertex(id, e[m][1])), bary);
      cell_basis(cell, bary, Cell::n_dofs_max, N.data(), "test");
      for (std::size_t j = 0; j < Cell::n_dofs_max; ++j){
        INFO("edge " << e[m][0] << e[m][1] << ", slot " << j);
        REQUIRE(N[j] == Approx(int(j) == Cell::mid_[m] ? 1. : 0.).margin(1e-12));
      }
    }
    // A partition of unity, whose derivatives sum to zero
    for (int trial = 0; trial < 5; ++trial){
      cell.contains(inside_point(c, id, gen), bary);
      for (const std::size_t n : {std::size_t(Cell::n_verts), Cell::n_dofs_max}){
        std::array<double, Cell::n_dofs_max> dNx, dNy, dNz;
        cell_basis(cell, bary, n, N.data(), "test");
        cell_deriv(cell, bary, n, dNx.data(), dNy.data(), dNz.data(), "test");
        double s = 0., sx = 0., sy = 0., sz = 0.;
        for (std::size_t j = 0; j < n; ++j){
          s += N[j]; sx += dNx[j]; sy += dNy[j];
          if (dim_of<Cell> == 3) sz += dNz[j];
        }
        REQUIRE(s == Approx(1.).margin(1e-12));
        REQUIRE(sx == Approx(0.).margin(1e-9));
        REQUIRE(sy == Approx(0.).margin(1e-9));
        REQUIRE(sz == Approx(0.).margin(1e-9));
      }
    }
  }
}

// The loaders' evaluation, from dolfin's interpolation of a polynomial: exact
// when the element holds the polynomial, and not when it does not
template<typename Cell>
double evaluation_error(const std::string& element, const bool quadratic, const bool gradient){
  constexpr int D = dim_of<Cell>;
  Cells<Cell> c(unit_mesh<Cell>(3, true));
  std::shared_ptr<dolfin::FunctionSpace> u_space, p_space;
  Uint ncoeffs_u = 0, ncoeffs_p = 0;
  taylor_hood_spaces<Cell>(element, "P1", false, c.mesh, nullptr, u_space, p_space, ncoeffs_u, ncoeffs_p);
  const Polynomial poly{D, quadratic};
  dolfin::Function u(u_space);
  u.interpolate(PolynomialExpression(poly));
  std::vector<double> values;
  u.vector()->get_local(values);
  CellDofs dofs;
  dofs.build(*u_space->dofmap(), c.dolfin_cells, "test");
  REQUIRE(dofs.stride() == std::size_t(D)*ncoeffs_u);

  std::mt19937 gen(2);
  std::uniform_real_distribution<> pick(0, c.cells.size() - 1);
  double err = 0.;
  for (int trial = 0; trial < 200; ++trial){
    const std::size_t start = std::size_t(pick(gen));
    const Vector3d x = inside_point(c, start, gen);
    CellPos pos{};
    REQUIRE(c.locate(x, pos));
    const Cell& cell = c.cells[pos.id];
    std::array<double, D*Cell::n_dofs_max> prev{}, next{};
    gather_stamps<D*Cell::n_verts, D*Cell::n_dofs_max>(dofs[pos.id], dofs.stride(), values, values,
                                                       prev.data(), next.data());
    // Both stamps read at the cell's dofs
    for (std::size_t i = 0; i < dofs.stride(); ++i){
      REQUIRE(prev[i] == values[dofs[pos.id][i]]);
      REQUIRE(next[i] == prev[i]);
    }
    if (!gradient){
      std::array<double, Cell::n_dofs_max> N;
      cell_basis(cell, pos.bary, ncoeffs_u, N.data(), "u");
      const Vector3d U = block_value<D>(N.data(), prev.data(), ncoeffs_u);
      for (int i = 0; i < D; ++i) err = std::max(err, std::abs(U[i] - poly.value(i, x)));
      for (int i = D; i < 3; ++i) REQUIRE(U[i] == 0.);
    }
    else {
      std::array<double, Cell::n_dofs_max> dNx, dNy, dNz;
      cell_deriv(cell, pos.bary, ncoeffs_u, dNx.data(), dNy.data(), dNz.data(), "u");
      const Matrix3d G = block_gradient<D>(dNx.data(), dNy.data(), dNz.data(), prev.data(), ncoeffs_u);
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j){
          const double exact = (i < D && j < D) ? poly.derivative(i, j, x) : 0.;
          err = std::max(err, std::abs(G(i, j) - exact));
        }
    }
  }
  return err;
}

template<typename Cell>
void check_evaluation(){
  // P1 holds a linear field, P2 a quadratic one: value and gradient exact
  REQUIRE(evaluation_error<Cell>("P1", false, false) < 1e-12);
  REQUIRE(evaluation_error<Cell>("P1", false, true) < 1e-10);
  REQUIRE(evaluation_error<Cell>("P2", true, false) < 1e-12);
  REQUIRE(evaluation_error<Cell>("P2", true, true) < 1e-10);
  // and P1 does not hold a quadratic one, so the checks above can fail
  REQUIRE(evaluation_error<Cell>("P1", true, false) > 1e-3);
}

template<typename Cell>
void check_locate(){
  Cells<Cell> c(unit_mesh<Cell>(4, true));
  std::mt19937 gen(3);
  const FoundCounts& n = c.found;
  for (int trial = 0; trial < 20; ++trial){
    const std::size_t id = std::size_t(gen() % c.cells.size());
    // Found in the cell it starts from, with barycentrics that rebuild the point
    const Vector3d x = inside_point(c, id, gen);
    CellPos pos{};
    pos.id = int(id);
    const auto same = n.same;
    REQUIRE(c.locate(x, pos));
    REQUIRE(pos.id == int(id));
    REQUIRE(n.same == same + 1);
    Vector3d back = Vector3d::Zero();
    for (int k = 0; k < Cell::n_verts; ++k) back += pos.bary[k]*c.vertex(id, k);
    REQUIRE((back - x).norm() < 1e-12);
    // In a face neighbour: one step of the walk, without the tree
    std::int32_t nb = -1;
    for (int k = 0; nb < 0; ++k) nb = c.across[id*Cell::n_verts + k];
    const Vector3d y = inside_point(c, nb, gen);
    const auto walk = n.walk;
    pos.id = int(id);
    REQUIRE(c.locate(y, pos));
    REQUIRE(pos.id == int(nb));
    REQUIRE(n.walk == walk + 1);
  }
  // Far away: the walk or the tree, the right cell either way
  CellPos pos{};
  pos.id = 0;
  const std::size_t far = c.cells.size() - 1;
  const Vector3d z = inside_point(c, far, gen);
  const auto other = n.tree, walk = n.walk;
  REQUIRE(c.locate(z, pos));
  REQUIRE(pos.id == int(far));
  REQUIRE(n.tree + n.walk == other + walk + 1);
  // A few cells on: the walk
  for (int trial = 0; trial < 20; ++trial){
    const std::size_t id = std::size_t(gen() % c.cells.size());
    std::size_t to = id;
    for (int s = 0; s < 3; ++s){
      const std::int32_t a = c.across[to*Cell::n_verts + gen() % Cell::n_verts];
      if (a >= 0) to = std::size_t(a);
    }
    const Vector3d y = inside_point(c, to, gen);
    const auto tree = n.tree;
    pos.id = int(id);
    REQUIRE(c.locate(y, pos));
    REQUIRE(pos.id == int(to));
    REQUIRE(n.tree == tree);
  }
  // No known cell: the tree
  const auto tree = n.tree;
  pos.id = -1;
  REQUIRE(c.locate(z, pos));
  REQUIRE(pos.id == int(far));
  REQUIRE(n.tree == tree + 1);
  // Outside the mesh: not found, and the cell is kept
  pos.id = 5;
  REQUIRE_FALSE(c.locate(Vector3d(-1., -1., dim_of<Cell> == 3 ? -1. : 0.), pos));
  REQUIRE(pos.id == 5);
}

// Walls on every side, or periodic along x
template<typename Cell>
struct Walk {
  Cells<Cell> c;
  std::vector<std::int32_t> across;
  Vector3d period;
  std::vector<bool> periodic;
  explicit Walk(const bool periodic_x) : c(unit_mesh<Cell>(4, false)), periodic{periodic_x, false, false} {
    const Vector3d lo(0., 0., 0.), hi(1., 1., dim_of<Cell> == 3 ? 1. : 0.);
    build_facet_neighbours(across, c.mesh, c.dolfin_cells, nullptr, periodic, lo, hi, dim_of<Cell>, 1e-12);
    period = Vector3d::Zero();
    if (periodic[0]) period[0] = hi[0] - lo[0];
  }
  // As the loaders' _modx
  Vector3d wrap(Vector3d x) const {
    if (periodic[0]) x[0] -= std::floor(x[0]);
    return x;
  }
  bool reflect(const Vector3d& x, Vector3d& dx, CellPos& pos){
    REQUIRE(c.locate(x, pos));
    return reflect_in_cells(c.cells, across, Cell::n_verts, period, x, dx, pos,
                            [this](const Vector3d& p){ return wrap(p); });
  }
};

template<typename Cell>
void check_reflect(){
  const double z = dim_of<Cell> == 3 ? 0.45 : 0.;
  Walk<Cell> w(false);
  CellPos pos{};
  // A move inside is unchanged
  Vector3d dx(0.1, 0.05, 0.);
  REQUIRE(w.reflect(Vector3d(0.3, 0.4, z), dx, pos));
  REQUIRE((dx - Vector3d(0.1, 0.05, 0.)).norm() < 1e-14);
  // Through the wall x = 0: the normal part mirrored, the tangential kept
  pos = CellPos();
  dx = Vector3d(-0.15, 0.02, 0.);
  REQUIRE(w.reflect(Vector3d(0.05, 0.4, z), dx, pos));
  REQUIRE((dx - Vector3d(0.05, 0.02, 0.)).norm() < 1e-12);
  std::array<double, 4> bary;
  REQUIRE(w.c.cells[pos.id].contains(Vector3d(0.1, 0.42, z), bary));
  // Into a corner: both walls
  pos = CellPos();
  dx = Vector3d(-0.15, -0.2, 0.);
  REQUIRE(w.reflect(Vector3d(0.05, 0.1, z), dx, pos));
  REQUIRE((dx - Vector3d(0.05, 0., 0.)).norm() < 1e-12);
}

template<typename Cell>
void check_periodic_walk(){
  const double z = dim_of<Cell> == 3 ? 0.45 : 0.;
  Walk<Cell> w(true);
  // Every facet on x = 0 and x = 1 has a partner
  std::size_t on_x = 0;
  for (dolfin::FacetIterator f(*w.c.mesh); !f.end(); ++f){
    const double x = f->midpoint().x();
    if (f->exterior() && (x < 1e-12 || x > 1. - 1e-12)) ++on_x;
  }
  REQUIRE(std::count_if(w.across.begin(), w.across.end(),
                        [](const std::int32_t a){ return a <= -2; }) == std::ptrdiff_t(on_x));
  // locate: a point in a periodic partner is found there, without the tree
  std::mt19937 gen(5);
  const FoundCounts& n = w.c.found;
  std::size_t partnered = 0;
  for (std::size_t id = 0; id < w.c.cells.size(); ++id)
    for (int k = 0; k < Cell::n_verts; ++k){
      const std::int32_t a = w.across[id*Cell::n_verts + k];
      if (a > -2) continue;
      const int j = facet_periodic(a);
      const Vector3d y = inside_point(w.c, std::size_t(j), gen);
      CellPos pos{};
      pos.id = int(id);
      const auto walk = n.walk, other = n.tree;
      REQUIRE(locate_in_cells(w.c.cells, w.across, *w.c.mesh, dim_of<Cell>, y, pos, nullptr, &w.c.found));
      REQUIRE(pos.id == j);
      REQUIRE(n.walk == walk + 1);
      REQUIRE(n.tree == other);
      ++partnered;
      // Beside the partner: found, by the walk or the tree
      for (int m = 0; m < Cell::n_verts; ++m){
        const std::int32_t beside = w.across[std::size_t(j)*Cell::n_verts + m];
        if (beside < 0) continue;
        const Vector3d z = inside_point(w.c, std::size_t(beside), gen);
        pos.id = int(id);
        REQUIRE(locate_in_cells(w.c.cells, w.across, *w.c.mesh, dim_of<Cell>, z, pos, nullptr, &w.c.found));
        REQUIRE(pos.id == beside);
      }
    }
  REQUIRE(partnered == on_x);
  // Through x = 1: the move is kept whole and ends in the cell of its image
  CellPos pos{};
  Vector3d dx(0.1, 0.03, 0.);
  REQUIRE(w.reflect(Vector3d(0.95, 0.4, z), dx, pos));
  REQUIRE((dx - Vector3d(0.1, 0.03, 0.)).norm() < 1e-12);
  std::array<double, 4> bary;
  REQUIRE(w.c.cells[pos.id].contains(Vector3d(0.05, 0.43, z), bary));
  // and back through x = 0
  pos = CellPos();
  dx = Vector3d(-0.1, 0., 0.);
  REQUIRE(w.reflect(Vector3d(0.05, 0.4, z), dx, pos));
  REQUIRE(w.c.cells[pos.id].contains(Vector3d(0.95, 0.4, z), bary));
}

// The cell tables on a unit mesh, built from dolfin as DolfInterpol builds them,
// without fields or files
template<typename Cell>
struct BareMesh : public MeshCore<Cell> {
  std::shared_ptr<dolfin::Mesh> mesh;
  std::vector<dolfin::Cell> dolfin_cells_;
  const Cell& cell(const int id) const { return this->cells_[id]; }
  BareMesh(std::shared_ptr<dolfin::Mesh> m, const std::vector<bool>& periodic) : MeshCore<Cell>("") {
    mesh = m;
    this->periodic = periodic;
    this->dim = m->geometry().dim();
    m->init();
    const std::vector<double> xx = m->coordinates();
    for (Uint i = 0; i < this->dim; ++i){ this->x_min[i] = xx[i]; this->x_max[i] = xx[i]; }
    for (Uint i = 0; i < xx.size(); ++i){
      const Uint d = i % this->dim;
      this->x_min[d] = std::min(this->x_min[d], xx[i]);
      this->x_max[d] = std::max(this->x_max[d], xx[i]);
    }
    this->hmin_ = m->hmin();
    for (dolfin::CellIterator c(*m); !c.end(); ++c){
      dolfin_cells_.push_back(*c);
      this->cells_.push_back(Cell(*c));
    }
    build_facet_neighbours(this->facet_neigh_, mesh, dolfin_cells_, nullptr, this->periodic,
                           this->x_min, this->x_max, this->dim, this->periodic_tol);
    this->set_period();
  }
  // The walk first; the fallback is dolfin's tree, as DolfInterpol's is
  bool locate(const Vector3d& x, const double, CellPos& pos){
    const Vector3d xx = this->_modx(x);
    return walk_to_cell(this->cells_, this->facet_neigh_, xx, pos)
        || tree_to_cell(this->cells_, *mesh, this->dim, xx, pos);
  }
  using MeshCore<Cell>::locate;
  void update(const double) {}
  void evaluate(const Vector3d&, const double, const CellPos&, PointValues&) {}
  double get_t_min() { return 0.; }
  double get_t_max() { return 1.; }
};

template<typename Cell>
void check_wall_normals(const bool periodic_x){
  auto mesh = unit_mesh<Cell>(3, false);
  BareMesh<Cell> m(mesh, {periodic_x, false, false});
  std::size_t on_wall = 0;
  for (std::size_t id = 0; id < m.dolfin_cells_.size(); ++id){
    // dolfin's outward normals of the cell's exterior facets, periodic ones left out
    Vector3d expected = Vector3d::Zero();
    for (dolfin::FacetIterator f(m.dolfin_cells_[id]); !f.end(); ++f){
      const double x = f->midpoint().x();
      if (!f->exterior() || (periodic_x && (x < 1e-12 || x > 1. - 1e-12))) continue;
      Vector3d n = Vector3d::Zero();
      for (int d = 0; d < dim_of<Cell>; ++d) n[d] = f->normal(d);
      expected += n;
    }
    if (expected.norm() > 0.){ expected.normalize(); ++on_wall; }
    int cell_id = int(id);
    const Vector3d n = m.get_boundary_normal(Vector3d::Zero(), cell_id);
    INFO("cell " << id);
    REQUIRE((n - expected).norm() < 1e-12);
  }
  REQUIRE(on_wall > 0);
  int outside = -1;
  REQUIRE(m.get_boundary_normal(Vector3d::Zero(), outside).norm() == 0.);
}

}  // namespace

TEST_CASE("The wall normal is the outward normal of the cell's wall facets", "[interpol]") {
  SECTION("triangle, walls") { check_wall_normals<Triangle>(false); }
  SECTION("tet, walls") { check_wall_normals<Tet>(false); }
  SECTION("triangle, periodic in x") { check_wall_normals<Triangle>(true); }
  SECTION("tet, periodic in x") { check_wall_normals<Tet>(true); }
}

TEST_CASE("P1 and P2 bases are nodal and a partition of unity", "[interpol]") {
  SECTION("triangle") { check_nodal_bases<Triangle>(); }
  SECTION("tet") { check_nodal_bases<Tet>(); }
}

TEST_CASE("The P1/P2 evaluation reproduces the element's polynomials", "[interpol]") {
  SECTION("triangle") { check_evaluation<Triangle>(); }
  SECTION("tet") { check_evaluation<Tet>(); }
}

TEST_CASE("locate finds the cell, walks to it, or goes to the tree", "[interpol]") {
  SECTION("triangle") { check_locate<Triangle>(); }
  SECTION("tet") { check_locate<Tet>(); }
}

TEST_CASE("The reflecting walk mirrors a move at the walls", "[interpol]") {
  SECTION("triangle") { check_reflect<Triangle>(); }
  SECTION("tet") { check_reflect<Tet>(); }
}

TEST_CASE("The walk crosses a periodic boundary into the image cell", "[interpol]") {
  SECTION("triangle") { check_periodic_walk<Triangle>(); }
  SECTION("tet") { check_periodic_walk<Tet>(); }
}

// MeshCore::locate, wrap included, on a box periodic in x and y
template<typename Cell>
void check_periodic_locate(){
  const double z = dim_of<Cell> == 3 ? 0.45 : 0.;
  BareMesh<Cell> m(unit_mesh<Cell>(4, false), {true, true, false});
  std::array<double, 4> bary;
  const auto found_at = [&](const Vector3d& from, const Vector3d& to, const Vector3d& image){
    CellPos pos{};
    REQUIRE(m.locate(from, 0., pos));
    const int start = pos.id;
    REQUIRE(m.locate(to, 0., pos));
    INFO("from cell " << start << " to cell " << pos.id);
    REQUIRE(m.cell(pos.id).contains(image, bary));
  };
  // Through a face, an edge of the box and its corner, and several boxes away
  found_at(Vector3d(0.97, 0.40, z), Vector3d(1.03, 0.41, z), Vector3d(0.03, 0.41, z));
  found_at(Vector3d(0.97, 0.97, z), Vector3d(1.03, 1.02, z), Vector3d(0.03, 0.02, z));
  found_at(Vector3d(0.02, 0.03, z), Vector3d(-0.04, -0.02, z), Vector3d(0.96, 0.98, z));
  found_at(Vector3d(0.30, 0.40, z), Vector3d(7.30, -3.60, z), Vector3d(0.30, 0.40, z));
  // Inside the box a point keeps its bits: same cell, same barycentrics
  CellPos a{}, b{};
  REQUIRE(m.locate(Vector3d(0.3125, 0.4375, z), 0., a));
  b = a;
  REQUIRE(m.locate(Vector3d(0.3125, 0.4375, z), 0., b));
  REQUIRE(a.id == b.id);
  for (int k = 0; k < Cell::n_verts; ++k) REQUIRE(a.bary[k] == b.bary[k]);
  // Not periodic in z: outside
  if (dim_of<Cell> == 3){
    CellPos pos{};
    REQUIRE_FALSE(m.locate(Vector3d(0.5, 0.5, 1.2), 0., pos));
  }
}

TEST_CASE("locate wraps a point through a face, an edge and a corner of a periodic box", "[interpol]") {
  SECTION("triangle") { check_periodic_locate<Triangle>(); }
  SECTION("tet") { check_periodic_locate<Tet>(); }
}

TEST_CASE("A mesh or element the evaluation cannot take throws partrac::Error", "[interpol][errors]") {
  // More dofs than the buffers hold
  REQUIRE_THROWS_AS(check_dofs_fit(11, 4, Tet::n_dofs_max, "test"), partrac::Error);
  REQUIRE_NOTHROW(check_dofs_fit(10, 4, Tet::n_dofs_max, "test"));
  // A dof table narrower than evaluate reads
  Cells<Tet> c(unit_mesh<Tet>(2, false));
  std::shared_ptr<dolfin::FunctionSpace> u_space, p_space;
  Uint ncoeffs_u = 0, ncoeffs_p = 0;
  taylor_hood_spaces<Tet>("P1", "P1", false, c.mesh, nullptr, u_space, p_space, ncoeffs_u, ncoeffs_p);
  CellDofs dofs;
  dofs.build(*u_space->dofmap(), c.dolfin_cells, "test");
  REQUIRE_THROWS_AS(dofs.check_stride(dofs.stride() + 1, "test"), partrac::Error);
  // An unknown element name
  REQUIRE_THROWS_AS(taylor_hood_spaces<Tet>("P7", "P1", false, c.mesh, nullptr, u_space, p_space,
                                            ncoeffs_u, ncoeffs_p), partrac::Error);
  // An unknown renumber_cells
  std::vector<std::uint32_t> map;
  REQUIRE_THROWS_WITH(cell_order(*u_space->dofmap(), c.cells.size(), "sometimes", map),
                      Catch::Contains("renumber_cells must be auto, never or always"));
}

TEST_CASE("The felbm wall normal points into the solid nodes next to a point", "[interpol]") {
  // 8^3 lattice, unit spacing, walls at x = 0 and x = 7, a solid node at (3, 3, 0) on the z = 0 face
  const Uint n[3] = {8, 8, 8};
  const Vector3d dx(1., 1., 1.);
  const auto solid = [](const Uint i, const Uint j, const Uint k){
    return i == 0 || i == 7 || (i == 3 && j == 3 && k == 0);
  };
  // Bulk: no solid neighbour
  REQUIRE(lattice_wall_normal(Vector3d(4., 4., 4.), dx, n, solid).norm() == 0.);
  // Next to x = 0 and next to x = 7, outward from the fluid
  REQUIRE((lattice_wall_normal(Vector3d(1.2, 4., 4.), dx, n, solid) - Vector3d(-1., 0., 0.)).norm() < 1e-14);
  REQUIRE((lattice_wall_normal(Vector3d(5.9, 4., 4.), dx, n, solid) - Vector3d(1., 0., 0.)).norm() < 1e-14);
  // At an edge between x = 0 and the node below, periodic in z: the mean of the two
  const Vector3d edge = lattice_wall_normal(Vector3d(1., 3., 7.), dx, n, [](const Uint i, const Uint, const Uint k){
    return i == 0 || k == 0;
  });
  REQUIRE((edge - Vector3d(-1., 0., 1.).normalized()).norm() < 1e-14);
}

TEST_CASE("PeriodicBC identifies the boundary dofs of a periodic space", "[interpol]") {
  // P1 on a 4 x 4 square: 25 vertices, 5 identified per periodic direction
  auto mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitSquareMesh(4, 4));
  const Vector3d lo(0., 0., 0.), hi(1., 1., 0.);
  std::shared_ptr<dolfin::FunctionSpace> u_space, p_space;
  Uint nu = 0, np = 0;
  const auto dofs = [&](const std::vector<bool>& periodic){
    auto bc = std::make_shared<PeriodicBC>(periodic, lo, hi, 2);
    taylor_hood_spaces<Triangle>("P1", "P1", true, mesh, bc, u_space, p_space, nu, np);
    return p_space->dim();
  };
  REQUIRE(dofs({false, false, false}) == 25);
  REQUIRE(dofs({true, false, false}) == 20);
  REQUIRE(dofs({true, true, false}) == 16);
}

#endif

// The OpenFOAM split and its stand-in node values, on meshes built here or
// read from the checked-in fixtures by a test reader: no OpenFOAM and no
// dolfin needed. The split must conform -- every interior facet in two
// simplices of opposite orientation, every boundary facet in one, every
// cyclic facet paired with its image -- or the facet table refuses the mesh
// or leaves holes; every simplex must be positive and a cell's simplices
// must fill it; and it must equal the reference split in split_expected.h5,
// simplex for simplex.
#include <catch2/catch.hpp>

#include <cmath>
#include <map>
#include <string>

#include "foam_arrays.hpp"
#include "h5direct.hpp"
#include "openfoam_nodes.hpp"
#include "openfoam_phase.hpp"
#include "openfoam_split.hpp"

namespace {

using openfoam_split::SplitData;

// Six times the signed volume (tets) or twice the signed area (triangles) of simplex t
double measure(const SplitData& s, const std::size_t t){
  const auto x = [&](const int j, const int d){
    return s.node_x[std::size_t(s.gdim)*s.cells[std::size_t(s.nv)*t + std::size_t(j)] + std::size_t(d)];
  };
  if (s.nv == 3)
    return (x(1, 0) - x(0, 0))*(x(2, 1) - x(0, 1)) - (x(1, 1) - x(0, 1))*(x(2, 0) - x(0, 0));
  double a[3][3];
  for (int j = 0; j < 3; ++j)
    for (int d = 0; d < 3; ++d) a[j][d] = x(j + 1, d) - x(0, d);
  return a[0][0]*(a[1][1]*a[2][2] - a[1][2]*a[2][1]) - a[0][1]*(a[1][0]*a[2][2] - a[1][2]*a[2][0])
       + a[0][2]*(a[1][0]*a[2][1] - a[1][1]*a[2][0]);
}

// Every simplex positive, and a 3D cell's tets summing to its volume
void check_positive(const openfoam_load::CaseData& c, const SplitData& s, const bool planar = true){
  std::vector<double> sum(std::size_t(c.ncells), 0.);
  for (std::size_t t = 0; t < s.nsimplices(); ++t){
    const double m = measure(s, t);
    REQUIRE(m > 0.);
    sum[std::size_t(s.cell_of[t])] += m/6.;
  }
  if (s.nv == 4 && planar)
    for (std::int64_t i = 0; i < c.ncells; ++i)
      REQUIRE(std::abs(sum[std::size_t(i)] - c.cell_volumes[std::size_t(i)]) <= 1e-12*c.cell_volumes[std::size_t(i)]);
}

// Interior facets twice, oppositely oriented; boundary facets once, with
// their patch; a cyclic facet's partner at its vertices' images
void check_conforming(const openfoam_load::CaseData& c, const SplitData& s){
  const int nv = s.nv;
  std::map<std::vector<std::uint32_t>, std::vector<std::pair<std::size_t, int>>> facets;
  for (std::size_t t = 0; t < s.nsimplices(); ++t)
    for (int k = 0; k < nv; ++k){
      std::vector<std::uint32_t> f;
      for (int j = 0; j < nv; ++j) if (j != k) f.push_back(s.cells[std::size_t(nv)*t + std::size_t(j)]);
      // the facet's orientation: the parity of its sort, flipped with the vertex it faces
      int parity = k % 2;
      for (std::size_t a = 0; a < f.size(); ++a)
        for (std::size_t b = a + 1; b < f.size(); ++b) parity ^= f[a] > f[b];
      std::sort(f.begin(), f.end());
      facets[f].push_back({std::size_t(nv)*t + std::size_t(k), parity});
    }
  std::size_t boundary = 0;
  for (const auto& kv : facets){
    const auto& v = kv.second;
    REQUIRE(v.size() <= 2);
    if (v.size() == 2){
      REQUIRE(v[0].second != v[1].second);
      REQUIRE(s.facet_patch[v[0].first] < 0);
      REQUIRE(s.facet_patch[v[1].first] < 0);
    }
    else {
      REQUIRE(s.facet_patch[v[0].first] >= 0);
      ++boundary;
    }
  }
  std::size_t expected = 0;
  for (const auto& p : c.patches)
    for (std::int64_t f = p.start; f < p.start + p.size; ++f)
      if (s.nv == 4) expected += std::size_t(c.face_start[std::size_t(f) + 1] - c.face_start[std::size_t(f)] - 2);
      else if (p.type != "empty") ++expected;
  REQUIRE(boundary == expected);
  // Cyclic facets: each paired, the partner at the images of its vertices
  for (std::size_t q = 0; q < s.facet_partner.size(); ++q){
    const std::int32_t patch = s.facet_patch[q];
    const bool cyclic = patch >= 0 && c.patches[std::size_t(patch)].type == "cyclic";
    REQUIRE((s.facet_partner[q] >= 0) == cyclic);
    if (!cyclic) continue;
    const auto& sep = c.patches[std::size_t(patch)].separation;
    const std::size_t o = std::size_t(s.facet_partner[q]);
    REQUIRE(s.facet_partner[o] == std::int64_t(q));
    REQUIRE(s.facet_patch[o] == c.patches[std::size_t(patch)].neighbour);
    std::vector<std::array<double, 3>> here, there;
    for (int j = 0; j < nv; ++j){
      if (j != int(q % std::size_t(nv))){
        std::array<double, 3> x{0., 0., 0.};
        for (int d = 0; d < s.gdim; ++d)
          x[std::size_t(d)] = s.node_x[std::size_t(s.gdim)*s.cells[(q/std::size_t(nv))*std::size_t(nv) + std::size_t(j)] + std::size_t(d)]
                            + sep[std::size_t(s.gdim == 3 ? d : s.inplane[std::size_t(d)])];
        here.push_back(x);
      }
      if (j != int(o % std::size_t(nv))){
        std::array<double, 3> x{0., 0., 0.};
        for (int d = 0; d < s.gdim; ++d)
          x[std::size_t(d)] = s.node_x[std::size_t(s.gdim)*s.cells[(o/std::size_t(nv))*std::size_t(nv) + std::size_t(j)] + std::size_t(d)];
        there.push_back(x);
      }
    }
    for (const auto& x : here){
      bool found = false;
      for (const auto& y : there)
        found = found || (std::abs(x[0] - y[0]) + std::abs(x[1] - y[1]) + std::abs(x[2] - y[2]) < 1e-12);
      REQUIRE(found);
    }
  }
}

// The unit cube's corners, relabelled: point perm[c] is corner c = x + 2y + 4z;
// each face listed from rot[f] onwards, which on a unit cube is its base point
std::vector<std::vector<std::int32_t>> cube_faces(const std::array<int, 8>& perm, const std::array<int, 6>& rot){
  const auto idx = [&](int x, int y, int z){ return std::int32_t(perm[std::size_t(x + 2*y + 4*z)]); };
  std::vector<std::vector<std::int32_t>> f = {
    {idx(0, 0, 0), idx(0, 0, 1), idx(0, 1, 1), idx(0, 1, 0)},
    {idx(1, 0, 0), idx(1, 1, 0), idx(1, 1, 1), idx(1, 0, 1)},
    {idx(0, 0, 0), idx(1, 0, 0), idx(1, 0, 1), idx(0, 0, 1)},
    {idx(0, 1, 0), idx(0, 1, 1), idx(1, 1, 1), idx(1, 1, 0)},
    {idx(0, 0, 0), idx(0, 1, 0), idx(1, 1, 0), idx(1, 0, 0)},
    {idx(0, 0, 1), idx(1, 0, 1), idx(1, 1, 1), idx(0, 1, 1)}};
  for (std::size_t i = 0; i < 6; ++i) std::rotate(f[i].begin(), f[i].begin() + rot[i], f[i].end());
  return f;
}

std::vector<std::array<double, 3>> cube_points(const std::array<int, 8>& perm){
  std::vector<std::array<double, 3>> p(8);
  for (int c = 0; c < 8; ++c) p[std::size_t(perm[std::size_t(c)])] = {double(c & 1), double((c >> 1) & 1), double((c >> 2) & 1)};
  return p;
}

// 64-bit FNV-1a over the int64 little-endian bytes of the values, as expected_split.py hashes
std::uint64_t fnv1a(const std::vector<std::uint32_t>& v){
  std::uint64_t h = 0xcbf29ce484222325ull;
  for (const std::uint32_t x : v){
    const std::int64_t y = x;
    for (int b = 0; b < 8; ++b){
      h ^= std::uint64_t((y >> (8*b)) & 0xff);
      h *= 0x100000001b3ull;
    }
  }
  return h;
}

template<typename T>
T h5_attribute(const hid_t loc, const char* name, const hid_t type){
  const partrac::H5Id a(H5Aopen(loc, name, H5P_DEFAULT), H5Aclose);
  REQUIRE(a.valid());
  T v{};
  REQUIRE(H5Aread(a, type, &v) >= 0);
  return v;
}

// The split of a fixture against its reference split in split_expected.h5,
// which tests/openfoam_expected.py writes from the fixture by the split's
// rules, independently of this code
void check_fixture(const std::string& name, const bool planar){
  const std::string dir = std::string(PARTRAC_SOURCE_DIR) + "/data_example/" + name;
  const openfoam_load::CaseData c = foam_arrays::read_case(dir);
  const partrac::H5Id file = partrac::h5_open_read(dir + "/split_expected.h5");
  for (const int n : {12, 6}){
    INFO(name << " split " << n);
    const SplitData s = openfoam_split::split(c, n);
    const std::string g = "w" + std::to_string(n);
    const auto same = [&](const std::string& ds, const auto& mine){
      std::vector<std::int64_t> v;
      partrac::h5_read(file, g + "/" + ds, v);
      REQUIRE(v.size() == mine.size());
      for (std::size_t i = 0; i < v.size(); ++i)
        if (v[i] != std::int64_t(mine[i])) FAIL(ds << "[" << i << "] " << mine[i] << ", expected " << v[i]);
    };
    same("cell_of", s.cell_of);
    same("node_kind", s.node_kind);
    same("node_point", s.node_point);
    same("node_cell", s.node_cell);
    same("facet_patch", s.facet_patch);
    const partrac::H5Id group(H5Gopen2(file, g.c_str(), H5P_DEFAULT), H5Gclose);
    if (H5Aexists(group, "cells_fnv1a") > 0){
      std::vector<std::int64_t> head;
      partrac::h5_read(file, g + "/cells_head", head);
      for (std::size_t i = 0; i < head.size(); ++i) REQUIRE(std::int64_t(s.cells[i]) == head[i]);
      REQUIRE(h5_attribute<std::int64_t>(group, "cells_size", H5T_NATIVE_INT64) == std::int64_t(s.cells.size()));
      REQUIRE(h5_attribute<std::uint64_t>(group, "cells_fnv1a", H5T_NATIVE_UINT64) == fnv1a(s.cells));
    }
    else {
      same("cells", s.cells);
    }
    check_positive(c, s, planar);
    check_conforming(c, s);
  }
}

}  // namespace

TEST_CASE("The OpenFOAM geometry is OpenFOAM's on a skewed cell", "[openfoam]") {
  // A unit cube with its top raised at one corner: faces not all planar
  auto p = cube_points({0, 1, 2, 3, 4, 5, 6, 7});
  p[7][2] = 1.5;
  const auto c = foam_arrays::single_cell(p, cube_faces({0, 1, 2, 3, 4, 5, 6, 7}, {0, 0, 0, 0, 0, 0}));
  // The volume under the faces' centre-point fans, as makeCellCentresAndVols has it
  REQUIRE(std::abs(c.cell_volumes[0] - 1.125) < 0.01);
  const auto u = foam_arrays::single_cell(cube_points({0, 1, 2, 3, 4, 5, 6, 7}),
                                          cube_faces({0, 1, 2, 3, 4, 5, 6, 7}, {0, 0, 0, 0, 0, 0}));
  REQUIRE(std::abs(u.cell_volumes[0] - 1.) < 1e-15);
  for (int d = 0; d < 3; ++d) REQUIRE(std::abs(u.cell_centres[std::size_t(d)] - 0.5) < 1e-15);
}

TEST_CASE("The OpenFOAM split conforms on a two-hex mesh and a 4^3 lattice", "[openfoam]") {
  for (const auto& n : {std::array<int, 3>{2, 1, 1}, std::array<int, 3>{4, 4, 4}})
    for (const unsigned rot : {0u, 7u})
      for (const int w : {12, 6}){
        INFO(n[0] << "x" << n[1] << "x" << n[2] << " rot " << rot << " W" << w);
        const auto c = foam_arrays::lattice(n, {false, false, false}, rot ? 3u : 0u, rot, rot ? 0.2 : 0.);
        const SplitData s = openfoam_split::split(c, w);
        check_positive(c, s, rot == 0);
        check_conforming(c, s);
        if (w == 12) REQUIRE(s.nsimplices() == 12*std::size_t(c.ncells));
      }
}

TEST_CASE("The OpenFOAM split pairs every cyclic facet with its image under a shuffled numbering", "[openfoam]") {
  for (const int w : {12, 6}){
    INFO("W" << w);
    const auto c = foam_arrays::lattice({5, 5, 4}, {true, true, false}, 11u, 5u);
    const SplitData s = openfoam_split::split(c, w);
    check_positive(c, s);
    check_conforming(c, s);
    std::size_t paired = 0;
    for (const auto o : s.facet_partner) paired += o >= 0;
    // two triangles a cyclic quad, on four sides of 20 quads
    REQUIRE(paired == 2*4*20);
    // the eight corners of the box share one master, the edge points along z four a column
    std::map<std::uint32_t, int> images;
    for (std::size_t n = 0; n < s.nnodes(); ++n) ++images[s.node_master[n]];
    int corners = 0;
    for (const auto& kv : images) corners += kv.second == 4;
    REQUIRE(corners == 5);
  }
}

TEST_CASE("Dompierre's split takes each of the eight cases of the far faces' diagonals", "[openfoam]") {
  // Corner 0 is point 0, the lowest; its three faces start there, so their
  // diagonals meet at it. Each far face starts at the opposite corner 7 or
  // next to it: its diagonal passes 7 or not. No far diagonal through 7: 5 tets.
  const std::array<int, 8> perm = {0, 1, 2, 3, 4, 5, 6, 7};
  for (int mask = 0; mask < 8; ++mask){
    INFO("mask " << mask);
    // far faces: x = 1 (1, starts at 1), y = 1 (3, starts at 2), z = 1 (5, starts at 4)
    std::array<int, 6> rot = {0, 0, 0, 0, 0, 0};
    const auto start_at = [&](const std::vector<std::int32_t>& f, const std::int32_t p){
      return int(std::find(f.begin(), f.end(), p) - f.begin());
    };
    const auto plain = cube_faces(perm, rot);
    rot[1] = start_at(plain[1], (mask & 1) ? 7 : 3);
    rot[3] = start_at(plain[3], (mask & 2) ? 7 : 3);
    rot[5] = start_at(plain[5], (mask & 4) ? 7 : 5);
    const auto c = foam_arrays::single_cell(cube_points(perm), cube_faces(perm, rot));
    const SplitData s = openfoam_split::split(c, 6);
    REQUIRE(s.nsimplices() == (mask == 0 ? 5u : 6u));
    REQUIRE(s.fan_cells == 0);
    REQUIRE(s.fallback == 0);
    check_positive(c, s);
    check_conforming(c, s);
  }
  // Diagonals that meet at no corner: the hex is fanned
  const std::array<int, 6> twisted = {1, 1, 0, 0, 0, 0};
  const auto c = foam_arrays::single_cell(cube_points(perm), cube_faces(perm, twisted));
  const SplitData s = openfoam_split::split(c, 6);
  REQUIRE(s.fan_cells == 1);
  REQUIRE(s.nsimplices() == 12);
  check_positive(c, s);
  check_conforming(c, s);
}

TEST_CASE("A hex whose centre fan folds takes Dompierre's split under W12", "[openfoam]") {
  // Corners 2 and 7 of the unit cube pulled in: faces so warped that the fan
  // from the cell centre has a tet under 1e-12 of the cell, while Dompierre's
  // tets from a corner are positive and fill it, so W12 falls back to them
  const std::array<int, 8> perm = {0, 1, 2, 3, 4, 5, 6, 7};
  auto p = cube_points(perm);
  p[2] = {0., 0.5, -0.5};
  p[7] = {0.25, 0.5, 0.5};
  const auto c = foam_arrays::single_cell(p, cube_faces(perm, {0, 0, 0, 0, 0, 0}));
  const SplitData s = openfoam_split::split(c, 12);
  REQUIRE(s.fallback == 1);
  REQUIRE(s.fan_cells == 0);
  REQUIRE(s.nsimplices() == 6);
  check_positive(c, s, false);
  check_conforming(c, s);
}

TEST_CASE("W12 is twelve tets a hex, the cell centre first", "[openfoam]") {
  const auto c = foam_arrays::lattice({3, 2, 2});
  const SplitData s = openfoam_split::split(c, 12);
  REQUIRE(s.nsimplices() == 12*std::size_t(c.ncells));
  REQUIRE(s.fan_cells == std::size_t(c.ncells));
  for (std::size_t t = 0; t < s.nsimplices(); ++t)
    REQUIRE(s.node_cell[s.cells[4*t]] == s.cell_of[t]);
  for (std::size_t n = 0; n < s.nnodes(); ++n)
    if (s.node_kind[n])
      for (int d = 0; d < 3; ++d)
        REQUIRE(s.node_x[3*n + std::size_t(d)] == c.cell_centres[3*std::size_t(s.node_cell[n]) + std::size_t(d)]);
  check_positive(c, s);
}

TEST_CASE("A prism and a pyramid are fanned whatever the split", "[openfoam]") {
  const auto prism = foam_arrays::single_cell({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}},
                                              {{0, 2, 1}, {3, 4, 5}, {0, 1, 4, 3}, {1, 2, 5, 4}, {0, 3, 5, 2}});
  const auto pyramid = foam_arrays::single_cell({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, 1}},
                                                {{0, 3, 2, 1}, {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}});
  REQUIRE(std::abs(prism.cell_volumes[0] - 0.5) < 1e-14);
  REQUIRE(std::abs(pyramid.cell_volumes[0] - 1./3.) < 1e-14);
  for (const int w : {12, 6}){
    const SplitData a = openfoam_split::split(prism, w);
    REQUIRE(a.nsimplices() == 8);
    REQUIRE(a.fan_cells == 1);
    check_positive(prism, a);
    check_conforming(prism, a);
    const SplitData b = openfoam_split::split(pyramid, w);
    REQUIRE(b.nsimplices() == 6);
    check_positive(pyramid, b);
    check_conforming(pyramid, b);
  }
}

TEST_CASE("A face with a hanging point is fanned from a base that leaves no flat tet", "[openfoam]") {
  // Point 8 halves the edge 4-5 on the top and the y = 0 faces, both listed
  // from point 4: their fans from point 4 hold three collinear points
  const auto c = foam_arrays::single_cell(
    {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}, {0.5, 0, 1}},
    {{0, 3, 2, 1}, {4, 8, 5, 6, 7}, {4, 0, 1, 5, 8}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}});
  const auto base = openfoam_split::face_bases(c);
  REQUIRE(base[1] != 0);
  REQUIRE(base[2] != 0);
  for (const int w : {12, 6}){
    const SplitData s = openfoam_split::split(c, w);
    REQUIRE(s.base_failures == 0);
    REQUIRE(s.invalid_cells == 0);
    check_positive(c, s);
    check_conforming(c, s);
    for (std::size_t t = 0; t < s.nsimplices(); ++t) REQUIRE(measure(s, t)/6. > 1e-3);
  }
}

TEST_CASE("The 2D split is the front plane: four triangles a quad, or two", "[openfoam]") {
  const auto c = foam_arrays::read_case(std::string(PARTRAC_SOURCE_DIR) + "/data_example/openfoam_cavity");
  REQUIRE(c.empty_axis == 2);
  for (const int w : {12, 6}){
    const SplitData s = openfoam_split::split(c, w);
    REQUIRE(s.nv == 3);
    REQUIRE(s.inplane == std::vector<int>{0, 1});
    REQUIRE(s.nsimplices() == (w == 12 ? 1600u : 800u));
    REQUIRE(s.nnodes() == (w == 12 ? 841u : 441u));
    check_positive(c, s);
    check_conforming(c, s);
  }
}

TEST_CASE("The fixtures split as their reference splits, simplex for simplex", "[openfoam]") {
  check_fixture("openfoam_cavity", true);
  // jittered: warped faces, whose cells OpenFOAM measures otherwise
  check_fixture("openfoam_channel3d", false);
}

TEST_CASE("The stand-in node values are volPointInterpolation's on a lattice", "[openfoam]") {
  // A linear field: inverse distance is exact at interior points of a uniform
  // lattice and at a wall point amid its wall faces, and the images of a
  // cyclic point read one value, exact for a field periodic across them; a
  // centre node carries its cell's value
  const auto c = foam_arrays::lattice({4, 4, 4}, {true, false, false});
  const SplitData s = openfoam_split::split(c, 12);
  // periodic in x, so constant along it
  const auto f = [](const double* x, const int d){ return 1. + 2.*x[1]*(d == 1) - 0.3*x[2]*(d + 1); };
  openfoam_load::FieldData u;
  u.name = "U";
  u.ncomp = 3;
  for (std::int64_t i = 0; i < c.ncells; ++i)
    for (int d = 0; d < 3; ++d) u.internal.push_back(f(&c.cell_centres[3*std::size_t(i)], d));
  for (const auto& p : c.patches){
    openfoam_load::FieldPatch fp;
    fp.condition = p.type == "cyclic" ? "cyclic" : "fixedValue";
    fp.fixes_value = fp.has_value = p.type != "cyclic";
    if (fp.has_value)
      for (std::int64_t q = p.start; q < p.start + p.size; ++q)
        for (int d = 0; d < 3; ++d) fp.values.push_back(f(&c.bface_centres[3*std::size_t(q - c.n_internal)], d));
    u.patches.push_back(fp);
  }
  const openfoam_nodes::Geometry g(c);
  const openfoam_nodes::InverseDistance idw(g, s);
  std::vector<double> v;
  idw.apply(u, {0, 1, 2}, v);
  int interior = 0, wall = 0;
  for (std::size_t n = 0; n < s.nnodes(); ++n){
    const double* x = &s.node_x[3*n];
    if (s.node_kind[n]){
      for (int d = 0; d < 3; ++d) REQUIRE(v[3*n + std::size_t(d)] == u.internal[3*std::size_t(s.node_cell[n]) + std::size_t(d)]);
      continue;
    }
    const auto inner = [](const double a){ return a > 1e-9 && a < 1. - 1e-9; };
    const bool on_x = !inner(x[0]);
    if (inner(x[1]) && inner(x[2])){
      // interior, or on the cyclic sides, whose images are interior too
      for (int d = 0; d < 3; ++d) REQUIRE(std::abs(v[3*n + std::size_t(d)] - f(x, d)) < 1e-14);
      ++interior;
      if (on_x){
        const std::size_t m = s.node_master[n];
        for (int d = 0; d < 3; ++d) REQUIRE(v[3*n + std::size_t(d)] == v[3*m + std::size_t(d)]);
      }
    }
    else if (inner(x[0]) && (inner(x[1]) != inner(x[2]))){
      // a wall point off the wall's edges: the mean of its four wall faces
      for (int d = 0; d < 3; ++d) REQUIRE(std::abs(v[3*n + std::size_t(d)] - f(x, d)) < 1e-14);
      ++wall;
    }
  }
  REQUIRE(interior == 5*3*3);
  REQUIRE(wall == 3*3*4);
}

// The least-squares node values W, rule by rule, on hex lattices built here
// and on the fixtures. W's promise is exactness for linear fields everywhere
// -- interior, face, edge and corner points -- whatever the conditions,
// which needs the next ring at boundary points (their cells' centres are
// coplanar), the fixed
// faces' values where the condition fixes them and never where it does not,
// mirrored cells at symmetry planes, and the cells across cyclic seams; and
// no-slip nodes exactly at rest, so the near-wall rule fires. The inverse-
// distance stand-in, cellPoint's oracle, is pinned alongside.
#include <functional>
#include <numeric>
#include <set>
#include <sstream>
#include <tuple>

namespace {

using openfoam_load::CaseData;
using openfoam_load::FieldData;
using openfoam_nodes::LeastSquares;
using openfoam_nodes::LeastSquaresReport;
using V3d = std::array<double, 3>;
using FieldFn = std::function<V3d(const double*)>;
using Conds = std::map<std::string, std::string>;

// Every mesh point a node, in point order; with centres, every cell's centre after them (W12's nodes)
SplitData point_nodes(const CaseData& c, const bool centres = false){
  SplitData s;
  const std::size_t np = c.npoints();
  for (std::size_t p = 0; p < np; ++p){
    s.node_kind.push_back(0);
    s.node_point.push_back(std::int32_t(p));
    s.node_cell.push_back(-1);
  }
  if (centres)
    for (std::int64_t i = 0; i < c.ncells; ++i){
      s.node_kind.push_back(1);
      s.node_point.push_back(-1);
      s.node_cell.push_back(std::int32_t(i));
    }
  return s;
}

// fn's first k components at the cell centres, and at the face centres of the
// patches whose condition (conds, else the patch type's or zeroGradient) fixes the value
FieldData make_field(const CaseData& c, const int k, const FieldFn& fn, const Conds& conds = {}){
  FieldData f;
  f.name = k == 3 ? "U" : "p";
  f.time = "0";
  f.ncomp = k;
  for (std::int64_t i = 0; i < c.ncells; ++i){
    const V3d v = fn(&c.cell_centres[3*std::size_t(i)]);
    f.internal.insert(f.internal.end(), v.begin(), v.begin() + k);
  }
  for (const auto& p : c.patches){
    openfoam_load::FieldPatch fp;
    const auto it = conds.find(p.name);
    const bool constraint = p.type == "cyclic" || p.type == "empty" || p.type == "symmetry" || p.type == "symmetryPlane";
    fp.condition = it != conds.end() ? it->second : constraint ? p.type : "zeroGradient";
    fp.fixes_value = fp.has_value = fp.condition == "fixedValue" || fp.condition == "noSlip";
    if (fp.has_value)
      for (std::int64_t q = p.start; q < p.start + p.size; ++q){
        const V3d v = fn(&c.bface_centres[3*std::size_t(q - c.n_internal)]);
        fp.values.insert(fp.values.end(), v.begin(), v.begin() + k);
      }
    f.patches.push_back(fp);
  }
  return f;
}

std::vector<int> all_comps(const int k){
  std::vector<int> v(static_cast<std::size_t>(k));
  std::iota(v.begin(), v.end(), 0);
  return v;
}

// W's node values, k a node
std::vector<double> ls_values(const CaseData& c, const SplitData& s, const FieldData& f,
                              LeastSquaresReport* r = nullptr, const double tol = LeastSquares::default_rank_tol,
                              LeastSquares* keep = nullptr){
  const openfoam_nodes::Geometry g(c);
  const LeastSquares w(g, s, f, r, tol);
  std::vector<double> out;
  w.apply(f, all_comps(f.ncomp), out);
  if (keep) *keep = w;
  return out;
}

std::vector<double> idw_values(const CaseData& c, const SplitData& s, const FieldData& f){
  const openfoam_nodes::Geometry g(c);
  const openfoam_nodes::InverseDistance w(g, s);
  std::vector<double> out;
  w.apply(f, all_comps(f.ncomp), out);
  return out;
}

const double* point_x(const CaseData& c, const std::size_t p){ return &c.points[3*p]; }

bool near(const double a, const double b){ return std::abs(a - b) < 1e-9; }

// 0 interior, 1 face, 2 edge, 3 corner of the box
std::vector<int> point_class(const CaseData& c){
  V3d lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
  for (std::size_t p = 0; p < c.npoints(); ++p)
    for (int d = 0; d < 3; ++d){
      lo[std::size_t(d)] = std::min(lo[std::size_t(d)], point_x(c, p)[d]);
      hi[std::size_t(d)] = std::max(hi[std::size_t(d)], point_x(c, p)[d]);
    }
  std::vector<int> cls(c.npoints(), 0);
  for (std::size_t p = 0; p < c.npoints(); ++p)
    for (int d = 0; d < 3; ++d) cls[p] += near(point_x(c, p)[d], lo[std::size_t(d)]) || near(point_x(c, p)[d], hi[std::size_t(d)]);
  return cls;
}

// The largest error of a point's node value against fn, over its k components
double point_error(const CaseData& c, const std::vector<double>& u, const int k, const FieldFn& fn, const std::size_t p){
  const V3d v = fn(point_x(c, p));
  double e = 0.;
  for (int d = 0; d < k; ++d) e = std::max(e, std::abs(u[std::size_t(k)*p + std::size_t(d)] - v[std::size_t(d)]));
  return e;
}

// Linear fields, and fields with the symmetry of the planes y = 0 and z = 0
V3d lin_vec(const double* x){
  return {1 + 2*x[0] - 3*x[1] + 0.5*x[2], -1 + x[0] + x[1] - 2*x[2], 0.3 - x[0] + 4*x[1] + x[2]};
}
V3d lin_scal(const double* x){ return {2 - x[0] + 3*x[1] - 0.7*x[2], 0., 0.}; }
V3d sym_vec(const double* x){ return {1 + 2*x[0], 3*x[1], -2*x[2]}; }
V3d sym_scal(const double* x){ return {1.5 - 2*x[0], 0., 0.}; }
V3d sym_smooth(const double* x){
  return {std::cos(x[1])*std::cos(x[2])*std::exp(x[0]), std::sin(2*x[1])*std::cos(x[2]), std::cos(x[1])*std::sin(3*x[2])};
}
V3d smooth3(const double* x){ return {std::sin(2*x[0] + 1)*std::cos(3*x[1])*std::exp(x[2]), 0., 0.}; }

std::array<std::vector<double>, 3> scaled(const std::vector<double>& a, const double sa, const std::vector<double>& b,
                                          const double sb, const std::vector<double>& c, const double sc){
  std::array<std::vector<double>, 3> ax{a, b, c};
  const double s[3] = {sa, sb, sc};
  for (int d = 0; d < 3; ++d) for (double& v : ax[std::size_t(d)]) v *= s[d];
  return ax;
}

std::array<std::vector<double>, 3> uniform_axes(){
  return {foam_arrays::linspace(0, 1, 5), foam_arrays::linspace(0, 1.2, 6), foam_arrays::linspace(0, 0.8, 4)};
}
std::array<std::vector<double>, 3> graded_axes(){
  return scaled(foam_arrays::graded(6, 4.), 1., foam_arrays::graded(7, 0.1), 1.2, foam_arrays::graded(5, 2.), 0.8);
}

// The eigenvalues of a symmetric 4 x 4 matrix, ascending, by Jacobi
std::array<double, 4> eigenvalues4(std::array<std::array<double, 4>, 4> a){
  for (int sweep = 0; sweep < 100; ++sweep)
    for (int p = 0; p < 4; ++p)
      for (int q = p + 1; q < 4; ++q){
        if (a[std::size_t(p)][std::size_t(q)] == 0.) continue;
        const double th = (a[std::size_t(q)][std::size_t(q)] - a[std::size_t(p)][std::size_t(p)])/(2*a[std::size_t(p)][std::size_t(q)]);
        const double t = (th >= 0 ? 1. : -1.)/(std::abs(th) + std::sqrt(th*th + 1.));
        const double cs = 1./std::sqrt(t*t + 1.), sn = t*cs;
        for (std::size_t k = 0; k < 4; ++k){
          const double kp = a[k][std::size_t(p)], kq = a[k][std::size_t(q)];
          a[k][std::size_t(p)] = cs*kp - sn*kq;
          a[k][std::size_t(q)] = sn*kp + cs*kq;
        }
        for (std::size_t k = 0; k < 4; ++k){
          const double pk = a[std::size_t(p)][k], qk = a[std::size_t(q)][k];
          a[std::size_t(p)][k] = cs*pk - sn*qk;
          a[std::size_t(q)][k] = sn*pk + cs*qk;
        }
      }
  std::array<double, 4> l{a[0][0], a[1][1], a[2][2], a[3][3]};
  std::sort(l.begin(), l.end());
  return l;
}

}  // namespace

TEST_CASE("W is exact for linear fields at every point with Dirichlet and zeroGradient sides", "[openfoam]") {
  // Fixed faces take the value; zeroGradient faces are not data, and the
  // next ring keeps every point exact, interior to corner; a fixed node is the
  // patch's value whatever the cells hold, a fitted one moves with them
  for (const std::string kind : {"uniform", "graded", "jittered"}){
    const auto c = foam_arrays::box(kind == "graded" ? graded_axes() : uniform_axes(), {}, kind == "jittered" ? 0.2 : 0., 4);
    const auto s = point_nodes(c);
    const Conds conds{{"xmin", "fixedValue"}, {"ymax", "fixedValue"}};
    for (const int k : {3, 1}){
      INFO(kind << " k " << k);
      const FieldFn fn = k == 3 ? FieldFn(lin_vec) : FieldFn(lin_scal);
      FieldData f = make_field(c, k, fn, conds);
      LeastSquaresReport r;
      const auto u = ls_values(c, s, f, &r);
      const auto cls = point_class(c);
      for (int cl = 0; cl < 4; ++cl){
        double e = 0.;
        int n = 0;
        for (std::size_t p = 0; p < c.npoints(); ++p)
          if (cls[p] == cl){
            e = std::max(e, point_error(c, u, k, fn, p));
            ++n;
          }
        INFO("class " << cl);
        REQUIRE(n > 0);
        REQUIRE(e < 1e-12);
      }
      REQUIRE(std::count(r.ring.begin(), r.ring.end(), 0) > 0);
      REQUIRE(std::count(r.ring.begin(), r.ring.end(), 2) > 0);
      std::mt19937 rng(0);
      std::normal_distribution<double> nd;
      for (double& v : f.internal) v += nd(rng);
      const auto u2 = ls_values(c, s, f);
      for (std::size_t p = 0; p < c.npoints(); ++p)
        for (int d = 0; d < k; ++d){
          const std::size_t i = std::size_t(k)*p + std::size_t(d);
          if (r.ring[p] == 0) REQUIRE(u2[i] == u[i]);
          else REQUIRE(u2[i] != u[i]);
        }
    }
  }
}

TEST_CASE("W with mirrored cells gives a symmetric field no normal component on the planes", "[openfoam]") {
  // Each cell and its mirror carry opposite normal components at equal
  // weights, edges of two planes included; without the mirror it is O(h^2)
  for (const std::string type : {"symmetry", "symmetryPlane"}){
    INFO(type);
    const auto c = foam_arrays::box(graded_axes(), {{"ymin", type}, {"zmin", type}});
    const auto f = make_field(c, 3, sym_smooth, {{"xmax", "fixedValue"}});
    LeastSquaresReport r;
    const auto u = ls_values(c, point_nodes(c), f, &r);
    int ony = 0, onz = 0, both = 0;
    for (std::size_t p = 0; p < c.npoints(); ++p){
      if (r.ring[p] == 0) continue;
      const bool y0 = near(point_x(c, p)[1], 0.), z0 = near(point_x(c, p)[2], 0.);
      if (y0) REQUIRE(std::abs(u[3*p + 1]) < 1e-14);
      if (z0) REQUIRE(std::abs(u[3*p + 2]) < 1e-14);
      ony += y0;
      onz += z0;
      both += y0 && z0;
    }
    REQUIRE(ony > 20);
    REQUIRE(onz > 20);
    REQUIRE(both > 3);
  }
}

TEST_CASE("W is exact at symmetry planes for linear fields with the symmetry", "[openfoam]") {
  // Symmetry faces are not data; the mirrored cells enter the fit with the
  // mirrored value, so the face, edge (two planes) and corner (two planes and
  // a zeroGradient side) points stay exact, a vector's normal component zero
  for (const std::string kind : {"uniform", "graded"})
    for (const std::string type : {"symmetry", "symmetryPlane"}){
      const auto c = foam_arrays::box(kind == "graded" ? graded_axes() : uniform_axes(), {{"ymin", type}, {"zmin", type}});
      const auto cls = point_class(c);
      for (const int k : {3, 1}){
        INFO(kind << " " << type << " k " << k);
        const FieldFn fn = k == 3 ? FieldFn(sym_vec) : FieldFn(sym_scal);
        const auto f = make_field(c, k, fn, {{"xmax", "fixedValue"}});
        LeastSquares w;
        const auto u = ls_values(c, point_nodes(c), f, nullptr, LeastSquares::default_rank_tol, &w);
        for (int cl = 0; cl < 4; ++cl)
          for (std::size_t p = 0; p < c.npoints(); ++p)
            if (cls[p] == cl) REQUIRE(point_error(c, u, k, fn, p) < 1e-12);
        REQUIRE(w.mirrored > 0);
        if (k == 3)
          for (std::size_t p = 0; p < c.npoints(); ++p)
            if (near(point_x(c, p)[1], 0.)) REQUIRE(std::abs(u[3*p + 1]) < 1e-12);
      }
    }
}

TEST_CASE("W's rank test fires at exactly the boundary points of a lattice", "[openfoam]") {
  // A boundary point's cells have coplanar centres, an interior point's do not
  const auto c = foam_arrays::box({foam_arrays::linspace(0, 1, 6), foam_arrays::linspace(0, 1, 6), foam_arrays::linspace(0, 1, 6)});
  LeastSquaresReport r;
  ls_values(c, point_nodes(c), make_field(c, 1, lin_scal), &r);
  const auto cls = point_class(c);
  for (std::size_t p = 0; p < c.npoints(); ++p) REQUIRE(int(r.ring[p]) == (cls[p] > 0 ? 2 : 1));
}

TEST_CASE("W's rank test sends the nearly coplanar rings of a jittered boundary to the next ring", "[openfoam]") {
  // Jittered interior points leave a boundary point's cells off their plane
  // by a little: scaled to unit diagonal the test at 1e-3 still takes the next
  // ring, and the boundary error of a smooth field stays the uniform
  // lattice's; at 1e-8 the first ring passes and extrapolates across its thin
  // spread, a hundred times worse
  const std::array<std::vector<double>, 3> ax{foam_arrays::linspace(0, 1, 16), foam_arrays::linspace(0, 1, 16),
                                              foam_arrays::linspace(0, 1, 16)};
  const auto uni = foam_arrays::box(ax);
  const auto jit = foam_arrays::box(ax, {}, 0.2, 1);
  const auto cls = point_class(jit);
  std::map<std::string, double> e;
  for (const auto& run : {std::make_tuple("uniform", &uni, 1e-3), std::make_tuple("jitter", &jit, 1e-3),
                          std::make_tuple("jitter 1e-8", &jit, 1e-8)}){
    const std::string name = std::get<0>(run);
    const CaseData& c = *std::get<1>(run);
    const double tol = std::get<2>(run);
    INFO(name);
    LeastSquaresReport r;
    const auto u = ls_values(c, point_nodes(c), make_field(c, 1, smooth3), &r, tol);
    const auto lin = ls_values(c, point_nodes(c), make_field(c, 1, lin_scal), nullptr, tol);
    double eb = 0., el = 0.;
    int b2 = 0, b1 = 0, nb = 0;
    for (std::size_t p = 0; p < c.npoints(); ++p){
      el = std::max(el, point_error(c, lin, 1, lin_scal, p));
      if (cls[p] == 0){
        if (name == "jitter") REQUIRE(r.ring[p] == 1);
        continue;
      }
      eb = std::max(eb, point_error(c, u, 1, smooth3, p));
      ++nb;
      b2 += r.ring[p] == 2;
      b1 += r.ring[p] == 1;
    }
    e[name] = eb;
    REQUIRE(el < (tol > 1e-8 ? 1e-12 : 1e-7));
    if (name == "jitter") REQUIRE(double(b2) > 0.99*nb);
    if (name == "jitter 1e-8") REQUIRE(b1 > 1000);
  }
  INFO("uniform " << e["uniform"] << " jitter " << e["jitter"] << " jitter 1e-8 " << e["jitter 1e-8"]);
  REQUIRE(e["jitter"] < 1.5*e["uniform"]);
  REQUIRE(e["jitter 1e-8"] > 100*e["uniform"]);
}

TEST_CASE("W's rank test never fires inside an aspect-1000 boundary layer", "[openfoam]") {
  // Scaled to unit diagonal the normal matrix does not see the cells'
  // aspect; scaled by the mean distance alone its condition is the aspect squared
  const auto c = foam_arrays::box(scaled(foam_arrays::linspace(0, 1, 12), 1., foam_arrays::linspace(0, 1, 12), 1.,
                                         foam_arrays::graded(12, 1000.), 0.01));
  LeastSquaresReport r;
  const auto u = ls_values(c, point_nodes(c), make_field(c, 1, lin_scal), &r);
  const auto cls = point_class(c);
  for (std::size_t p = 0; p < c.npoints(); ++p){
    REQUIRE(r.ring[p] > 0);
    REQUIRE(point_error(c, u, 1, lin_scal, p) < 1e-12);
    if (cls[p] == 0){
      REQUIRE(r.ring[p] == 1);
      REQUIRE(r.cond[p] < 10.);
    }
  }
}

namespace {

// The lid-driven cavity: n x n, 2D (one cell between empty sides) or n^3
CaseData cavity_box(const int n, const bool two_d){
  Conds t{{"xmin", "wall"}, {"xmax", "wall"}, {"ymin", "wall"}, {"ymax", "wall"}};
  t["zmin"] = t["zmax"] = two_d ? "empty" : "wall";
  return foam_arrays::box({foam_arrays::linspace(0, 1, n), foam_arrays::linspace(0, 1, n),
                           two_d ? foam_arrays::linspace(0, 0.1, 1) : foam_arrays::linspace(0, 1, n)}, t);
}

// The lid (ymax) at speed U, the other walls no-slip, random cell values
FieldData cavity_field(const CaseData& c, const std::string& lid = "ymax", const double U = 1.){
  FieldData f;
  f.name = "U";
  f.ncomp = 3;
  std::mt19937 rng(3);
  std::uniform_real_distribution<double> ud(-1., 1.);
  for (std::int64_t i = 0; i < 3*c.ncells; ++i) f.internal.push_back(ud(rng));
  for (const auto& p : c.patches){
    openfoam_load::FieldPatch fp;
    if (p.type == "empty") fp.condition = "empty";
    else {
      fp.condition = p.name == lid ? "fixedValue" : "noSlip";
      fp.fixes_value = fp.has_value = true;
      for (std::int64_t q = 0; q < p.size; ++q){
        fp.values.push_back(p.name == lid ? U : 0.);
        fp.values.push_back(0.);
        fp.values.push_back(0.);
      }
    }
    f.patches.push_back(fp);
  }
  return f;
}

}  // namespace

TEST_CASE("W's zero wins where a moving lid meets the no-slip walls", "[openfoam]") {
  // The lid's corners are at rest (cellPoint gives half the lid speed there),
  // the lid's other points exactly at its speed, the walls exactly at rest
  for (const bool two_d : {true, false}){
    INFO((two_d ? "2D" : "3D"));
    const auto c = cavity_box(8, two_d);
    LeastSquaresReport r;
    const auto u = ls_values(c, point_nodes(c), cavity_field(c), &r);
    std::set<std::int32_t> corners;
    for (std::size_t p = 0; p < c.npoints(); ++p){
      const double* x = point_x(c, p);
      const bool top = near(x[1], 1.);
      bool side = near(x[0], 0.) || near(x[0], 1.);
      if (!two_d) side = side || near(x[2], 0.) || near(x[2], 1.);
      const bool wall = side || near(x[1], 0.);
      const V3d v{u[3*p], u[3*p + 1], u[3*p + 2]};
      if (top && side){
        REQUIRE(v == V3d{0., 0., 0.});
        corners.insert(std::int32_t(p));
      }
      else if (top) REQUIRE(v == V3d{1., 0., 0.});
      else if (wall) REQUIRE(v == V3d{0., 0., 0.});
    }
    REQUIRE(std::set<std::int32_t>(r.zero_wins.begin(), r.zero_wins.end()) == corners);
  }
}

TEST_CASE("The inverse-distance nodes give the lid's corners half its speed", "[openfoam]") {
  // cellPoint's blending of the lid and the wall on a uniform mesh, for contrast
  const auto c = cavity_box(8, true);
  const auto u = idw_values(c, point_nodes(c), cavity_field(c));
  int n = 0;
  for (std::size_t p = 0; p < c.npoints(); ++p){
    const double* x = point_x(c, p);
    if (!near(x[1], 1.) || !(near(x[0], 0.) || near(x[0], 1.))) continue;
    REQUIRE(std::abs(u[3*p] - 0.5) < 1e-15);
    REQUIRE(std::abs(u[3*p + 1]) < 1e-15);
    ++n;
  }
  REQUIRE(n == 4);
}

TEST_CASE("W reproduces a moving wall's linear profile along it", "[openfoam]") {
  // A wall whose value varies linearly: its nodes, the patch's edges and
  // corners included, take a fit over the patch's faces, their next ring
  // where the first is one-sided
  for (const bool two_d : {false, true}){
    INFO((two_d ? "2D" : "3D"));
    const auto c = foam_arrays::box({foam_arrays::graded(6, 3.), foam_arrays::linspace(0, 1, 4),
                                     two_d ? foam_arrays::linspace(0, 0.1, 1) : foam_arrays::graded(5, 3.)},
                                    two_d ? Conds{{"zmin", "empty"}, {"zmax", "empty"}} : Conds{});
    const FieldFn prof = [two_d](const double* x){
      return V3d{0.5 + 2*x[0] - (two_d ? 0. : 1.5*x[2]), 0., 0.2 + x[0]};
    };
    const auto u = ls_values(c, point_nodes(c), make_field(c, 3, prof, {{"ymax", "fixedValue"}}));
    int n = 0;
    for (std::size_t p = 0; p < c.npoints(); ++p){
      if (!near(point_x(c, p)[1], 1.)) continue;
      REQUIRE(point_error(c, u, 3, prof, p) < 1e-12);
      ++n;
    }
    REQUIRE(n == (two_d ? 14 : 42));
  }
}

TEST_CASE("Cyclic images chain across two cyclic pairs", "[openfoam]") {
  // A point on an x-y edge has four images, x + shift its master, the pairs
  // separated along their axes
  const auto c = foam_arrays::box({foam_arrays::linspace(0, 1, 4), foam_arrays::linspace(0, 1, 4), foam_arrays::linspace(0, 1, 4)},
                                  {{"xmin", "cyclic"}, {"xmax", "cyclic"}, {"ymin", "cyclic"}, {"ymax", "cyclic"}});
  const openfoam_nodes::Geometry g(c);
  std::map<std::int32_t, int> images;
  for (std::size_t p = 0; p < c.npoints(); ++p){
    REQUIRE(g.master[p] <= std::int32_t(p));
    for (int d = 0; d < 3; ++d)
      REQUIRE(std::abs(point_x(c, p)[d] + g.shift[3*p + std::size_t(d)] - point_x(c, std::size_t(g.master[p]))[d]) < 1e-15);
    ++images[g.master[p]];
  }
  for (std::size_t p = 0; p < c.npoints(); ++p)
    if (near(point_x(c, p)[0], 0.) && near(point_x(c, p)[1], 0.)) REQUIRE(images[g.master[p]] == 4);
  REQUIRE(c.patches[0].separation == std::array<double, 3>{1., 0., 0.});
  REQUIRE(c.patches[2].separation == std::array<double, 3>{0., 1., 0.});
}

namespace {

V3d periodic_field(const double* x){
  const double pi = std::acos(-1.);
  return {1 + std::cos(2*pi*x[0]) + 2*x[2], 0.5*std::sin(2*pi*x[1]) - x[2], 3.};
}

}  // namespace

TEST_CASE("W gives a point's cyclic images one value and a seam point an interior point's accuracy", "[openfoam]") {
  // Periodic in x and y, walls in z: the images read their master's row, whose
  // ring takes the cells on every side, so a seam point is its translate
  // mid-channel; with zeroGradient walls a wall point beside the seam takes
  // its next ring across it; linear in z is exact on the wall-seam edges
  for (const bool graded : {false, true}){
    INFO((graded ? "graded" : "uniform"));
    const auto xs = foam_arrays::linspace(0, 1, 8);
    const auto c = foam_arrays::box({xs, foam_arrays::linspace(0, 1, 8), graded ? foam_arrays::graded(6, 3.) : foam_arrays::linspace(0, 1, 6)},
                                    {{"xmin", "cyclic"}, {"xmax", "cyclic"}, {"ymin", "cyclic"}, {"ymax", "cyclic"}});
    const auto s = point_nodes(c);
    const openfoam_nodes::Geometry g(c);
    LeastSquaresReport r;
    const auto u = ls_values(c, s, make_field(c, 3, periodic_field, {{"zmin", "noSlip"}, {"zmax", "fixedValue"}}), &r);
    std::vector<double> seam, mid;
    for (std::size_t p = 0; p < c.npoints(); ++p){
      for (int d = 0; d < 3; ++d) REQUIRE(u[3*p + std::size_t(d)] == u[3*std::size_t(g.master[p]) + std::size_t(d)]);
      const double* x = point_x(c, p);
      if (x[2] <= 1e-9 || x[2] >= 1 - 1e-9 || !near(x[1], 0.5)) continue;
      if (near(x[0], 0.)){
        seam.push_back(point_error(c, u, 3, periodic_field, p));
        REQUIRE(r.ring[p] == 1);
      }
      if (near(x[0], 0.5)) mid.push_back(point_error(c, u, 3, periodic_field, p));
    }
    REQUIRE(seam.size() == mid.size());
    REQUIRE(!seam.empty());
    std::sort(seam.begin(), seam.end());
    std::sort(mid.begin(), mid.end());
    for (std::size_t i = 0; i < seam.size(); ++i) REQUIRE(std::abs(seam[i] - mid[i]) <= 1e-14 + 1e-8*std::abs(mid[i]));
    LeastSquaresReport rz;
    const auto uz = ls_values(c, s, make_field(c, 3, periodic_field), &rz);
    const double h = xs[1];
    double e_near = -1., e_far = -1.;
    for (std::size_t p = 0; p < c.npoints(); ++p){
      const double* x = point_x(c, p);
      if (!near(x[2], 0.) || !near(x[1], 0.5)) continue;
      if (near(x[0], h)){
        REQUIRE(rz.ring[p] == 2);
        e_near = point_error(c, uz, 3, periodic_field, p);
      }
      if (near(x[0], 0.5 + h)) e_far = point_error(c, uz, 3, periodic_field, p);
    }
    REQUIRE(e_near >= 0.);
    REQUIRE(std::abs(e_near - e_far) <= 1e-14 + 1e-8*e_far);
    const FieldFn lin = [](const double* x){ return V3d{1 + 2*x[2], -x[2], 3.}; };
    const auto u2 = ls_values(c, s, make_field(c, 3, lin, {{"zmin", "fixedValue"}, {"zmax", "zeroGradient"}}));
    for (std::size_t p = 0; p < c.npoints(); ++p) REQUIRE(point_error(c, u2, 3, lin, p) < 1e-12);
  }
}

TEST_CASE("The inverse-distance nodes are exact inside a uniform lattice", "[openfoam]") {
  const auto c = foam_arrays::box({foam_arrays::linspace(0, 1, 6), foam_arrays::linspace(0, 1, 6), foam_arrays::linspace(0, 1, 6)});
  const auto u = idw_values(c, point_nodes(c), make_field(c, 3, lin_vec));
  const auto cls = point_class(c);
  for (std::size_t p = 0; p < c.npoints(); ++p)
    if (cls[p] == 0) REQUIRE(point_error(c, u, 3, lin_vec, p) < 1e-13);
}

namespace {

std::vector<double> alternating(const int n){
  std::vector<double> x(std::size_t(n) + 1, 0.);
  for (int i = 0; i < n; ++i) x[std::size_t(i) + 1] = x[std::size_t(i)] + (i % 2 ? 2. : 1.);
  for (double& v : x) v /= x.back();
  return x;
}

// The inverse-distance nodes' largest error inside a one-cell-thick lattice, so thin that its distances are the plane's
double idw_interior_error(const std::vector<double>& xs, const std::vector<double>& ys, const FieldFn& fn){
  const auto c = foam_arrays::box({xs, ys, {0., 1e-7}}, {{"zmin", "empty"}, {"zmax", "empty"}});
  const auto u = idw_values(c, point_nodes(c), make_field(c, 1, fn));
  double e = 0.;
  for (std::size_t p = 0; p < c.npoints(); ++p){
    const double* x = point_x(c, p);
    if (near(x[0], 0.) || near(x[0], 1.) || near(x[1], 0.) || near(x[1], 1.)) continue;
    e = std::max(e, point_error(c, u, 1, fn, p));
  }
  return e;
}

}  // namespace

TEST_CASE("The inverse-distance nodes converge at second order on a graded lattice, first on an alternating one", "[openfoam]") {
  // A linear field on the graded lattice (ratios 4 and 10): 1.7e-2, 4.2e-3,
  // 1.1e-3, 2.6e-4 at n = 8 .. 64; a smooth one on h, 2h, h, 2h at rate 1,
  // where least squares keeps 2
  const FieldFn lin2 = [](const double* x){ return V3d{1 + 2*x[0] - 3*x[1], 0., 0.}; };
  const double expect[4] = {1.7e-2, 4.2e-3, 1.1e-3, 2.6e-4};
  int i = 0;
  for (const int n : {8, 16, 32, 64}){
    const double e = idw_interior_error(foam_arrays::graded(n, 4.), foam_arrays::graded(n, 10.), lin2);
    INFO("n " << n << ": " << e);
    REQUIRE(std::abs(e - expect[i]) < 0.06*expect[i]);
    ++i;
  }
  const FieldFn smooth2 = [](const double* x){ return V3d{std::sin(2*x[0] + 1)*std::cos(3*x[1]), 0., 0.}; };
  std::vector<double> e;
  for (const int n : {16, 32, 64}) e.push_back(idw_interior_error(alternating(n), alternating(n), smooth2));
  for (std::size_t j = 0; j + 1 < e.size(); ++j){
    const double rate = std::log2(e[j]/e[j + 1]);
    INFO("rate " << rate);
    REQUIRE(std::abs(rate - 1.) < 0.15);
  }
}

TEST_CASE("The inverse-distance nodes blend wall faces and sum a seam point's images", "[openfoam]") {
  // A wall point between two faces of different values takes their
  // 1/|x - x_f| mean and ignores the cells; zeroGradient faces carry their
  // cells' values; an empty side does not blend; a cyclic point sums its cells
  // on both sides
  const std::vector<double> xs{0., 0.3, 1.}, ys{0., 0.5, 1.2}, zs{0., 0.1};
  const auto c = foam_arrays::box({xs, ys, zs}, {{"zmin", "empty"}, {"zmax", "empty"}});
  const FieldFn seven = [](const double*){ return V3d{7., 0., 0.}; };
  FieldData f = make_field(c, 1, seven);
  f.patches[2].condition = "fixedValue";
  f.patches[2].fixes_value = f.patches[2].has_value = true;
  f.patches[2].values = {2., 5.};
  f.internal = {1., 3., 4., 8.};
  const auto u = idw_values(c, point_nodes(c), f);
  const auto find = [](const CaseData& m, const double x, const double y, const double z){
    for (std::size_t p = 0; p < m.npoints(); ++p)
      if (near(point_x(m, p)[0], x) && near(point_x(m, p)[1], y) && near(point_x(m, p)[2], z)) return p;
    FAIL("no point");
    return std::size_t(0);
  };
  const auto dist = [](const CaseData& m, const std::size_t p, const double* y){
    const double* x = point_x(m, p);
    return std::sqrt((x[0] - y[0])*(x[0] - y[0]) + (x[1] - y[1])*(x[1] - y[1]) + (x[2] - y[2])*(x[2] - y[2]));
  };
  const auto fc = [&](const CaseData& m, const std::int64_t f){ return &m.bface_centres[3*std::size_t(f - m.n_internal)]; };
  const std::size_t p = find(c, 0.3, 0., 0.);
  const auto& ymin = c.patches[2];
  const double w0 = 1/dist(c, p, fc(c, ymin.start)), w1 = 1/dist(c, p, fc(c, ymin.start + 1));
  REQUIRE(std::abs(u[p] - (2*w0 + 5*w1)/(w0 + w1)) < 1e-14*u[p]);
  const std::size_t q = find(c, 0., 0.5, 0.1);
  const auto& xmin = c.patches[0];
  const double v0 = 1/dist(c, q, fc(c, xmin.start)), v1 = 1/dist(c, q, fc(c, xmin.start + 1));
  REQUIRE(std::abs(u[q] - (1*v0 + 4*v1)/(v0 + v1)) < 1e-14*u[q]);
  const std::size_t r = find(c, 0.3, 0.5, 0.);
  double num = 0., den = 0.;
  for (std::size_t i = 0; i < 4; ++i){
    const double w = 1/dist(c, r, &c.cell_centres[3*i]);
    num += w*f.internal[i];
    den += w;
  }
  REQUIRE(std::abs(u[r] - num/den) < 1e-14*u[r]);
  const auto m = foam_arrays::box({xs, ys, zs}, {{"xmin", "cyclic"}, {"xmax", "cyclic"}, {"zmin", "empty"}, {"zmax", "empty"}});
  FieldData g = make_field(m, 1, seven);
  g.internal = {1., 3., 4., 8.};
  const auto uc = idw_values(m, point_nodes(m), g);
  const std::size_t a = find(m, 0., 0.5, 0.), b = find(m, 1., 0.5, 0.);
  const double wa0 = 1/dist(m, a, &m.cell_centres[0]), wa2 = 1/dist(m, a, &m.cell_centres[6]);
  const double wb1 = 1/dist(m, b, &m.cell_centres[3]), wb3 = 1/dist(m, b, &m.cell_centres[9]);
  const double want = (wa0*1 + wa2*4 + wb1*3 + wb3*8)/(wa0 + wa2 + wb1 + wb3);
  REQUIRE(std::abs(uc[a] - want) < 1e-14*want);
  REQUIRE(uc[b] == uc[a]);
}

TEST_CASE("The inverse-distance nodes lose the normal component at a symmetry plane", "[openfoam]") {
  const auto c = foam_arrays::box({foam_arrays::linspace(0, 1, 3), foam_arrays::linspace(0, 1, 3), foam_arrays::linspace(0, 1, 3)},
                                  {{"ymin", "symmetryPlane"}});
  const auto u = idw_values(c, point_nodes(c), make_field(c, 3, lin_vec));
  double uy = 0., ux = 0.;
  for (std::size_t p = 0; p < c.npoints(); ++p)
    if (near(point_x(c, p)[1], 0.)){
      uy = std::max(uy, std::abs(u[3*p + 1]));
      ux = std::max(ux, std::abs(u[3*p]));
    }
  REQUIRE(uy == 0.);
  REQUIRE(ux > 0.);
}

TEST_CASE("W's normal matrix is conditioned independently of the cell size", "[openfoam]") {
  // Scaled to unit diagonal, the interior condition number of a graded,
  // jittered lattice stays put under refinement; unscaled it grows as 1/h^2
  std::vector<double> conds, raw;
  for (const int n : {8, 16, 32}){
    const auto c = foam_arrays::box({foam_arrays::graded(n, 4.), foam_arrays::graded(n, 10.), foam_arrays::graded(n, 2.)},
                                    {}, 0.2, unsigned(n));
    LeastSquaresReport r;
    ls_values(c, point_nodes(c), make_field(c, 1, lin_scal), &r);
    const openfoam_nodes::Geometry g(c);
    const auto cls = point_class(c);
    double cmax = 0., rmax = 0.;
    for (std::size_t p = 0; p < c.npoints(); ++p){
      if (cls[p]) continue;
      cmax = std::max(cmax, r.cond[p]);
      std::array<std::array<double, 4>, 4> N{};
      for (std::size_t j = g.pc_start[p]; j < g.pc_start[p + 1]; ++j){
        const double* y = &c.cell_centres[3*g.pc[j]];
        const double A[4] = {1., y[0] - point_x(c, p)[0], y[1] - point_x(c, p)[1], y[2] - point_x(c, p)[2]};
        for (std::size_t a = 0; a < 4; ++a)
          for (std::size_t b = 0; b < 4; ++b) N[a][b] += A[a]*A[b];
      }
      const auto l = eigenvalues4(N);
      rmax = std::max(rmax, l[3]/l[0]);
    }
    conds.push_back(cmax);
    raw.push_back(rmax);
  }
  INFO(conds[0] << " " << conds[1] << " " << conds[2] << "; raw " << raw[0] << " " << raw[2]);
  REQUIRE(*std::max_element(conds.begin(), conds.end()) < 1.5*(*std::min_element(conds.begin(), conds.end())));
  REQUIRE(*std::max_element(conds.begin(), conds.end()) < 10.);
  REQUIRE(raw[2] > 10*raw[0]);
}

TEST_CASE("W on a 64^3 lattice: exact for a linear field, a centre its cell's value", "[openfoam]") {
  // Cyclic in x, walls in z, W12's nodes (the points, then every cell's
  // centre): the operator built on a lattice of 262,144 cells
  const int n = 64;
  const auto c = foam_arrays::box({foam_arrays::linspace(0, 1, n), foam_arrays::linspace(0, 1, n), foam_arrays::linspace(0, 1, n)},
                                  {{"xmin", "cyclic"}, {"xmax", "cyclic"}});
  const FieldFn lin = [](const double* x){ return V3d{1 + 2*x[1] - x[2], x[2], 3 - x[1]}; };
  const auto f = make_field(c, 3, lin, {{"zmin", "noSlip"}, {"zmax", "fixedValue"}});
  const auto s = point_nodes(c, true);
  const openfoam_nodes::Geometry g(c);
  const LeastSquares w(g, s, f);
  std::vector<double> u;
  w.apply(f, {0, 1, 2}, u);
  REQUIRE(u.size() == 3*(c.npoints() + std::size_t(n)*n*n));
  for (std::size_t p = 0; p < c.npoints(); ++p) REQUIRE(point_error(c, u, 3, lin, p) < 1e-12);
  for (std::size_t i = c.npoints(); i < s.nnodes(); ++i)
    for (int d = 0; d < 3; ++d) REQUIRE(u[3*i + std::size_t(d)] == f.internal[3*std::size_t(s.node_cell[i]) + std::size_t(d)]);
}

TEST_CASE("W refuses a stamp whose conditions differ from the first stamp's", "[openfoam]") {
  // The rows are fixed at load: a patch that fixes the value at one stamp and
  // not at another is refused, naming the patch and the stamp
  const auto c = cavity_box(4, false);
  const auto f = cavity_field(c);
  const openfoam_nodes::Geometry g(c);
  const LeastSquares w(g, point_nodes(c), f);
  FieldData later = f;
  later.time = "7";
  later.patches[3].condition = "zeroGradient";
  std::vector<double> u;
  REQUIRE_THROWS_WITH(w.apply(later, {0, 1, 2}, u), Catch::Contains("U at 7: patch ymax is zeroGradient"));
  later = f;
  later.time = "8";
  later.patches[0].condition = "inletOutlet";
  REQUIRE_THROWS_WITH(w.apply(later, {0, 1, 2}, u), Catch::Contains("U at 8: patch xmin is inletOutlet"));
}

TEST_CASE("W on the checked-in channel: linear in z exact at every point of the split", "[openfoam]") {
  // The jittered channel (cyclic in x and y, walls in z) through the test
  // reader and the W12 split: exact with fixed walls or zeroGradient ones
  // (the next ring), the images one value, a centre its cell's value
  const auto c = foam_arrays::read_case(std::string(PARTRAC_SOURCE_DIR) + "/data_example/openfoam_channel3d");
  const SplitData s = openfoam_split::split(c, 12);
  const FieldFn lin = [](const double* x){ return V3d{1 + 2*x[2], -3*x[2], 0.5}; };
  for (const Conds& conds : {Conds{{"bottom", "noSlip"}, {"top", "fixedValue"}}, Conds{}}){
    const auto f = make_field(c, 3, lin, conds);
    const auto u = ls_values(c, s, f);
    for (std::size_t n = 0; n < s.nnodes(); ++n){
      if (s.node_kind[n]){
        for (int d = 0; d < 3; ++d) REQUIRE(u[3*n + std::size_t(d)] == f.internal[3*std::size_t(s.node_cell[n]) + std::size_t(d)]);
        continue;
      }
      REQUIRE(point_error(c, u, 3, lin, std::size_t(s.node_point[n])) < 1e-12);
      for (int d = 0; d < 3; ++d) REQUIRE(u[3*n + std::size_t(d)] == u[3*std::size_t(s.node_master[n]) + std::size_t(d)]);
    }
  }
}

TEST_CASE("W on the checked-in cavity: the lid's corners at rest, the pressure fitted at the walls", "[openfoam]") {
  // Through the test reader and the 2D split: zero wins at the lid's corners
  // (the inverse-distance nodes give half the lid speed), the lid at its
  // speed; the pressure's walls are zeroGradient, so no node is fixed and the
  // wall nodes take the next ring
  const auto c = foam_arrays::read_case(std::string(PARTRAC_SOURCE_DIR) + "/data_example/openfoam_cavity");
  const SplitData s = openfoam_split::split(c, 12);
  const auto u = ls_values(c, s, cavity_field(c, "movingWall"));
  const auto ui = idw_values(c, s, cavity_field(c, "movingWall"));
  double lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300};
  for (std::size_t p = 0; p < c.npoints(); ++p)
    for (int d = 0; d < 2; ++d){
      lo[d] = std::min(lo[d], point_x(c, p)[d]);
      hi[d] = std::max(hi[d], point_x(c, p)[d]);
    }
  int corners = 0, lid = 0;
  for (std::size_t n = 0; n < s.nnodes(); ++n){
    if (s.node_kind[n]) continue;
    const double* x = point_x(c, std::size_t(s.node_point[n]));
    const bool top = near(x[1], hi[1]), side = near(x[0], lo[0]) || near(x[0], hi[0]);
    if (top && side){
      REQUIRE(u[3*n] == 0.);
      REQUIRE(std::abs(ui[3*n] - 0.5) < 1e-15);
      ++corners;
    }
    else if (top){
      REQUIRE(u[3*n] == 1.);
      ++lid;
    }
  }
  REQUIRE(corners == 2);
  REQUIRE(lid == 19);
  const FieldFn zero = [](const double*){ return V3d{0., 0., 0.}; };
  LeastSquaresReport r;
  ls_values(c, s, make_field(c, 1, zero, {{"frontAndBack", "empty"}}), &r);
  for (std::size_t n = 0; n < s.nnodes(); ++n){
    if (s.node_kind[n]) continue;
    const double* x = point_x(c, std::size_t(s.node_point[n]));
    const bool wall = near(x[1], hi[1]) || near(x[1], lo[1]) || near(x[0], lo[0]) || near(x[0], hi[0]);
    REQUIRE(int(r.ring[std::size_t(s.node_point[n])]) == (wall ? 2 : 1));
  }
}

// The phase field: the volume above 1/2 of a P1 field, closed form per
// simplex; each fanned cell's centre moved so its cell's volume above 1/2 is
// its volume fraction's, which W6's hexes, without a centre, cannot; and a
// gradient of its own, the least-squares fit's, exact for linear fields
namespace {

using openfoam_phase::above_half;

// The fraction of a simplex above 1/2 by Monte Carlo, uniform points by
// exponential spacings, and its standard error
std::pair<double, double> above_half_mc(const double* v, const int nv, std::mt19937& rng, const int n = 400000){
  std::exponential_distribution<double> ex;
  int hits = 0;
  for (int i = 0; i < n; ++i){
    double l[4], sum = 0., f = 0.;
    for (int j = 0; j < nv; ++j) sum += l[j] = ex(rng);
    for (int j = 0; j < nv; ++j) f += l[j]/sum*v[j];
    hits += f > 0.5;
  }
  const double p = double(hits)/n;
  return {p, std::sqrt(std::max(p*(1 - p), 1./n)/n)};
}

// Each OpenFOAM cell's centre node on the split, -1 for none
openfoam_phase::Cells phase_cells(const CaseData& c, const SplitData& s){
  std::vector<std::int32_t> centre(std::size_t(c.ncells), -1);
  for (std::size_t n = 0; n < s.nnodes(); ++n)
    if (s.node_kind[n]) centre[std::size_t(s.node_cell[n])] = std::int32_t(n);
  return openfoam_phase::Cells(std::size_t(c.ncells), s.cell_of, centre);
}

// Each cell's volume and its volume above 1/2 from node values v
std::pair<std::vector<double>, std::vector<double>> cell_volumes(const CaseData& c, const SplitData& s,
                                                                 const std::vector<double>& v){
  std::vector<double> vol(std::size_t(c.ncells), 0.), above(vol);
  const std::size_t k = std::size_t(s.nv);
  for (std::size_t t = 0; t < s.nsimplices(); ++t){
    double w[4];
    for (std::size_t j = 0; j < k; ++j) w[j] = v[s.cells[k*t + j]];
    const double m = measure(s, t)/(s.nv == 3 ? 2. : 6.);
    vol[std::size_t(s.cell_of[t])] += m;
    above[std::size_t(s.cell_of[t])] += m*above_half(w, s.nv);
  }
  return {vol, above};
}

// A unit lattice with walls, jittered, and a 2D box one cell thick between an empty pair
CaseData phase_lattice(const bool two_d){
  if (!two_d) return foam_arrays::lattice({7, 8, 6}, {false, false, false}, 0, 0, 0.25);
  return foam_arrays::box({foam_arrays::linspace(0, 1, 13), foam_arrays::linspace(0, 1, 11), {0., 0.1}},
                          {{"zmin", "empty"}, {"zmax", "empty"}}, 0.25);
}

// A tanh step 0.05 wide across a tilted plane, and across a sphere (a circle in 2D)
V3d tilted_step(const double* x){
  return {0.5*(1 + std::tanh((0.36*x[0] + 0.6*x[1] + 0.71*x[2] - 0.83)/0.05)), 0., 0.};
}
V3d round_step(const double* x){
  const double r = std::sqrt((x[0] - 0.45)*(x[0] - 0.45) + (x[1] - 0.52)*(x[1] - 0.52) + (x[2] - 0.55)*(x[2] - 0.55));
  return {0.5*(1 - std::tanh((r - 0.3)/0.05)), 0., 0.};
}
V3d round_step_2d(const double* x){
  const double y[3] = {x[0], x[1], 0.55};
  return round_step(y);
}

}  // namespace

TEST_CASE("The volume above 1/2 of a linear field on a simplex is exact for every sign pattern", "[openfoam]") {
  // Against the closed form for distinct values, sum_i max(s_i, 0)^d over
  // prod_{j != i} (s_i - s_j) with s = v - 1/2, the values kept 0.1 apart so
  // it does not cancel; against Monte Carlo for every pattern and where values
  // repeat or sit at 1/2, which does not count as above
  std::mt19937 rng(3);
  std::uniform_real_distribution<double> mag(0.05, 1.);
  for (const int nv : {3, 4})
    for (int pattern = 0; pattern < (1 << nv); ++pattern){
      INFO("nv " << nv << " pattern " << pattern);
      double v[4], s[4];
      for (int trial = 0; trial < 50; ++trial){
        bool apart = false;
        while (!apart){
          for (int i = 0; i < nv; ++i) s[i] = ((pattern >> i) & 1 ? 1. : -1.)*mag(rng);
          apart = true;
          for (int i = 0; i < nv; ++i)
            for (int j = i + 1; j < nv; ++j) apart = apart && std::abs(s[i] - s[j]) > 0.1;
        }
        double exact = 0.;
        for (int i = 0; i < nv; ++i){
          double den = 1.;
          for (int j = 0; j < nv; ++j) if (j != i) den *= s[i] - s[j];
          exact += std::pow(std::max(s[i], 0.), nv - 1)/den;
          v[i] = 0.5 + s[i];
        }
        REQUIRE(std::abs(above_half(v, nv) - exact) < 1e-12);
      }
      const auto mc = above_half_mc(v, nv, rng);
      REQUIRE(std::abs(above_half(v, nv) - mc.first) < 5*mc.second);
    }
  const std::vector<std::vector<double>> repeated = {
    {0.9, 0.9, 0.2, 0.2}, {0.7, 0.7, 0.7, 0.1}, {0.5, 0.9, 0.1, 0.3}, {0.5, 0.5, 0.5, 0.9}, {0.5, 0.5, 0.2, 0.2},
    {0.8, 0.8, 0.1}, {0.5, 0.9, 0.2}, {0.5, 0.5, 0.1}};
  for (const auto& v : repeated){
    const int nv = int(v.size());
    const auto mc = above_half_mc(v.data(), nv, rng);
    INFO(v[0] << " " << v[1] << " " << v[2]);
    REQUIRE(std::abs(above_half(v.data(), nv) - mc.first) < 5*mc.second);
  }
  // All at 1/2: nothing above, where Monte Carlo reads round-off
  const double half[4] = {0.5, 0.5, 0.5, 0.5};
  REQUIRE(above_half(half, 4) == 0.);
  REQUIRE(above_half(half, 3) == 0.);
}

TEST_CASE("A mixed fanned cell's centre brings its volume above 1/2 to its volume fraction's within [0, 1]", "[openfoam]") {
  // W12 on a jittered lattice and a 2D box, cellPoint's node values of a tanh
  // step across a tilted plane and across a sphere: a cell within purity of 0
  // or 1 keeps its centre; a corrected one is V_c alpha_c to 1e-12 of V_c;
  // a clipped one, which the 3D cases have, sits at 0 or 1 with its target
  // past what that end gives; the sum's error is no more than the pure and
  // clipped cells' errors; only centres move, and none leaves [0, 1]
  for (const bool two_d : {false, true})
    for (const FieldFn& fn : {FieldFn(tilted_step), FieldFn(two_d ? round_step_2d : round_step)}){
      INFO((two_d ? "2D" : "3D"));
      const CaseData c = phase_lattice(two_d);
      const SplitData s = openfoam_split::split(c, 12);
      REQUIRE(s.fan_cells == std::size_t(c.ncells));
      const FieldData f = make_field(c, 1, fn);
      const std::vector<double> v0 = idw_values(c, s, f);
      std::vector<double> v = v0;
      const auto r = openfoam_phase::conserve(phase_cells(c, s), s.cells, s.node_x, s.nv, f.internal, v);
      REQUIRE(r.pure + r.corrected + r.clipped + r.unreachable == std::size_t(c.ncells));
      REQUIRE(r.unreachable == 0);
      REQUIRE(r.unfanned == 0);
      INFO("pure " << r.pure << " corrected " << r.corrected << " clipped " << r.clipped);
      REQUIRE(r.pure > 0);
      REQUIRE(r.corrected > std::size_t(c.ncells)/10);
      if (!two_d) REQUIRE(r.clipped > 0);
      REQUIRE(r.worst_before > 1e-2);
      const auto [vol, above] = cell_volumes(c, s, v);
      std::vector<std::size_t> centre(vol.size());
      for (std::size_t n = 0; n < s.nnodes(); ++n){
        if (!s.node_kind[n]){
          REQUIRE(v[n] == v0[n]);
          continue;
        }
        centre[std::size_t(s.node_cell[n])] = n;
        REQUIRE(v[n] >= 0.);
        REQUIRE(v[n] <= 1.);
      }
      double target = 0., kept = 0.;
      std::size_t pure = 0, clipped = 0;
      for (std::size_t i = 0; i < vol.size(); ++i){
        const double a = f.internal[i], want = vol[i]*a, err = std::abs(above[i] - want);
        target += want;
        const std::size_t n = centre[i];
        if (a < openfoam_phase::purity || a > 1 - openfoam_phase::purity){
          REQUIRE(v[n] == v0[n]);
          kept += err;
          ++pure;
        }
        else if (err > 1e-12*vol[i]){
          REQUIRE((v[n] == 0. || v[n] == 1.));
          REQUIRE((v[n] == 0.) == (above[i] > want));
          kept += err;
          ++clipped;
        }
      }
      REQUIRE(pure == r.pure);
      REQUIRE(clipped == r.clipped);
      REQUIRE(std::abs(r.target - target) < 1e-13*target);
      REQUIRE(std::abs(r.after - r.target) <= kept + 1e-12*target);
      REQUIRE(std::abs(r.after - r.target) < std::abs(r.before - r.target));
    }
}

TEST_CASE("A cell of volume fraction 0 against a point past 1/2 keeps its centre, counted unreachable", "[openfoam]") {
  // Water but for a quarter of an unjittered lattice: the points on the
  // quarter's inner edge read 3/4, so the dry cells there cannot be dry for
  // any centre value; they keep theirs, and every other cell is pure
  const CaseData c = foam_arrays::lattice({6, 6, 3});
  const SplitData s = openfoam_split::split(c, 12);
  const FieldFn quarter = [](const double* x){ return V3d{x[0] > 0.5 && x[1] > 0.5 ? 0. : 1., 0., 0.}; };
  const FieldData f = make_field(c, 1, quarter);
  const std::vector<double> v0 = idw_values(c, s, f);
  std::vector<double> v = v0;
  const auto r = openfoam_phase::conserve(phase_cells(c, s), s.cells, s.node_x, s.nv, f.internal, v);
  REQUIRE(v == v0);
  REQUIRE(r.unreachable == 3);
  REQUIRE(r.pure == std::size_t(c.ncells) - 3);
  REQUIRE(r.corrected + r.clipped == 0);
  REQUIRE(r.worst_after > 0.);
}

TEST_CASE("Under W6 a hex has no centre and its nodes are left as they are", "[openfoam]") {
  // Dompierre's tets on every hex: no node moves, the error is reported, and
  // every cell counts as off without a centre where it is
  const CaseData c = phase_lattice(false);
  const SplitData s = openfoam_split::split(c, 6);
  REQUIRE(s.fan_cells == 0);
  const FieldData f = make_field(c, 1, round_step);
  const std::vector<double> v0 = idw_values(c, s, f);
  std::vector<double> v = v0;
  const auto r = openfoam_phase::conserve(phase_cells(c, s), s.cells, s.node_x, s.nv, f.internal, v);
  REQUIRE(v == v0);
  REQUIRE(r.pure + r.corrected + r.clipped + r.unreachable == 0);
  REQUIRE(r.unfanned > 0);
  REQUIRE(r.before == r.after);
  REQUIRE(r.worst_before == r.worst_after);
  REQUIRE(r.worst_before > 1e-2);
  const auto [vol, above] = cell_volumes(c, s, v);
  double sum = 0.;
  for (const double a : above) sum += a;
  REQUIRE(std::abs(r.after - sum) < 1e-12*sum);
}

TEST_CASE("The phase field's gradient is exact for a linear field at every node", "[openfoam]") {
  // At the points the fit's G, the next ring at the walls, across a cyclic
  // pair (the field constant along it); at a centre the mean of its points';
  // on the jittered lattice and the checked-in cavity (two components)
  const auto c3 = foam_arrays::box({foam_arrays::linspace(0, 1, 6), foam_arrays::linspace(0, 1, 7), foam_arrays::linspace(0, 1, 5)},
                                   {{"xmin", "cyclic"}, {"xmax", "cyclic"}}, 0.2);
  const auto c2 = foam_arrays::read_case(std::string(PARTRAC_SOURCE_DIR) + "/data_example/openfoam_cavity");
  for (const CaseData* c : {&c3, &c2}){
    const bool two_d = c->empty_axis >= 0;
    INFO((two_d ? "cavity" : "lattice"));
    const V3d G = two_d ? V3d{3., 5., 0.} : V3d{0., 2., -1.5};
    const FieldFn lin = [&](const double* x){ return V3d{0.2 + G[0]*x[0] + G[1]*x[1] + G[2]*x[2], 0., 0.}; };
    const SplitData s = openfoam_split::split(*c, 12);
    const openfoam_nodes::Geometry g(*c);
    const openfoam_nodes::Gradient grad(g, s);
    std::vector<double> out;
    grad.apply(make_field(*c, 1, lin), out);
    const int d = two_d ? 2 : 3;
    REQUIRE(out.size() == std::size_t(d)*s.nnodes());
    REQUIRE(grad.second_ring > 0);
    openfoam_phase::centre_means(phase_cells(*c, s), s.cells, s.nv, d, out);
    for (std::size_t n = 0; n < s.nnodes(); ++n)
      for (int a = 0; a < d; ++a) REQUIRE(std::abs(out[std::size_t(d)*n + std::size_t(a)] - G[std::size_t(a)]) < 1e-10);
  }
}

TEST_CASE("The patches are classed by the velocity's condition, and the split keeps each boundary facet's patch", "[openfoam]") {
  // A wall at rest is no-slip, one fixed along itself a moving wall, one
  // fixed across itself (an inlet), an outlet or a symmetry plane other; the
  // walk treats every class but cyclic as a wall, the classes are for the log
  // and for the open boundaries, so they must come out of the conditions alone
  const auto c = foam_arrays::box({foam_arrays::linspace(0, 1, 4), foam_arrays::linspace(0, 1, 4),
                                   foam_arrays::graded(4, 3.)},
                                  {{"xmin", "cyclic"}, {"xmax", "cyclic"}, {"zmin", "wall"}, {"zmax", "wall"}});
  const FieldFn inflow = [](const double* x){ return V3d{1e-3*x[0], 1. + x[2], 0.}; };
  auto f = make_field(c, 3, inflow, {{"ymin", "fixedValue"}, {"zmin", "noSlip"}, {"zmax", "fixedValue"}});
  for (std::size_t i = 0; i < f.patches[4].values.size(); ++i) f.patches[4].values[i] = i % 3 ? 0. : 1e-13;
  for (std::size_t i = 0; i < f.patches[5].values.size(); ++i) f.patches[5].values[i] = i % 3 ? 0. : 1.;
  using openfoam_nodes::PatchClass;
  const std::vector<PatchClass> k = openfoam_nodes::classify_patches(c, f);
  // xmin xmax ymin ymax zmin zmax
  REQUIRE(k == std::vector<PatchClass>{PatchClass::cyclic, PatchClass::cyclic, PatchClass::other,
                                       PatchClass::other, PatchClass::wall, PatchClass::moving_wall});
  // above round-off the bottom moves
  f.patches[4].values[0] = 1e-9;
  REQUIRE(openfoam_nodes::classify_patches(c, f)[4] == PatchClass::moving_wall);
  // a wall moving through itself is an inlet
  f.patches[5].values[2] = 1e-3;
  REQUIRE(openfoam_nodes::classify_patches(c, f)[5] == PatchClass::other);

  for (const int w : {12, 6}){
    INFO("W" << w);
    const SplitData s = openfoam_split::split(c, w);
    std::vector<std::size_t> per_patch(c.patches.size(), 0);
    for (const std::int32_t p : s.facet_patch)
      if (p >= 0) ++per_patch[std::size_t(p)];
    // every boundary face's two triangles from its base point, on 16 faces a side
    for (const std::size_t m : per_patch) REQUIRE(m == 2*16u);
  }
}

TEST_CASE("A case with every point on a no-slip wall is refused under W6 and warned of under W12", "[openfoam]") {
  // One cell thick between the no-slip walls zmin and zmax, cyclic in x and
  // y: every point is fixed at zero, so W6's tets, which have no other node,
  // give zero velocity throughout, and W12's only the cell centres carry the
  // flow. The case is not read as 2D, its walls are the case's. Two cells
  // thick, every cell has a free point and nothing is said
  const FieldFn poiseuille = [](const double* x){ return V3d{4*x[2]*(1 - x[2]), 0., 0.}; };
  for (const int nz : {1, 2}){
    INFO(nz << " cells thick");
    const auto c = foam_arrays::lattice({4, 4, nz}, {true, true, false});
    const auto f = make_field(c, 3, poiseuille, {{"zmin", "noSlip"}, {"zmax", "noSlip"}});
    const auto k = openfoam_nodes::classify_patches(c, f);
    REQUIRE(openfoam_nodes::walled_cells(c, k) == (nz == 1 ? 16u : 0u));
    std::ostringstream log;
    if (nz == 1)
      REQUIRE_THROWS_WITH(openfoam_nodes::check_walled(c, k, 6, log),
                          Catch::Contains("every mesh point lies on a no-slip wall")
                          && Catch::Contains("under split=6 every node is fixed at zero"));
    else
      openfoam_nodes::check_walled(c, k, 6, log);
    openfoam_nodes::check_walled(c, k, 12, log);
    REQUIRE(log.str().find("Warning: every mesh point lies on a no-slip wall") == (nz == 1 ? 0 : std::string::npos));
  }
}

TEST_CASE("Cells with every point on a no-slip wall are counted and warned of under either split", "[openfoam]") {
  // 4x4x1, cyclic in x, zmin a no-slip wall and zmax one only over its y = 0
  // row: that row's four cells have every point on a wall, the others a free
  // point at the top. A warning gives the count and the fraction, and what it
  // means under each split, with no refusal; zmax walled throughout is the
  // all-walls case, still refused under W6
  auto c = foam_arrays::lattice({4, 4, 1}, {true, false, false});
  openfoam_load::Patch strip = c.patches.back();
  strip.name = "zmax_wall";
  strip.size = 4;
  c.patches.back().start += 4;
  c.patches.back().size -= 4;
  c.patches.insert(c.patches.end() - 1, strip);
  using K = openfoam_nodes::PatchClass;
  std::vector<K> k{K::cyclic, K::cyclic, K::other, K::other, K::wall, K::wall, K::other};
  REQUIRE(openfoam_nodes::walled_cells(c, k) == 4u);
  for (const int w : {6, 12}){
    INFO("W" << w);
    std::ostringstream log;
    openfoam_nodes::check_walled(c, k, w, log);
    REQUIRE_THAT(log.str(), Catch::StartsWith("Warning: 4 of 16 cells (25%) have every point on a no-slip wall")
                 && Catch::Contains(w == 6 ? "under split=6 their velocity is zero throughout"
                                           : "under split=12 only their centre values are free"));
  }
  k[6] = K::wall;
  REQUIRE(openfoam_nodes::walled_cells(c, k) == 16u);
  std::ostringstream log;
  REQUIRE_THROWS_WITH(openfoam_nodes::check_walled(c, k, 6, log),
                      Catch::Contains("every mesh point lies on a no-slip wall"));
  openfoam_nodes::check_walled(c, k, 12, log);
  REQUIRE(log.str().find("Warning: every mesh point lies on a no-slip wall") == 0);
  k[5] = K::other;
  k[6] = K::other;
  REQUIRE(openfoam_nodes::walled_cells(c, k) == 0u);
}

TEST_CASE("Cyclic facets pair by the addressing, whatever the positions of their points", "[openfoam]") {
  // The two sides of a cyclic pair written with round-off between them, and
  // with a displacement no position match would accept: the facets pair all
  // the same, since only the cyclic faces' point order says which is whose image
  for (const double eps : {1e-15, 1e-6}){
    INFO("eps " << eps);
    auto c = foam_arrays::box({foam_arrays::linspace(0, 1, 4), foam_arrays::linspace(0, 1, 5),
                               foam_arrays::graded(3, 2.)},
                              {{"xmin", "cyclic"}, {"xmax", "cyclic"}, {"ymin", "cyclic"}, {"ymax", "cyclic"}});
    std::mt19937 rng(5);
    std::uniform_real_distribution<double> ud(-1., 1.);
    for (std::size_t p = 0; p < c.npoints(); ++p)
      if (c.points[3*p] == 1. || c.points[3*p + 1] == 1.)
        for (int d = 0; d < 3; ++d) c.points[3*p + std::size_t(d)] += eps*ud(rng);
    foam_arrays::derive_geometry(c);
    for (const int w : {12, 6}){
      const SplitData s = openfoam_split::split(c, w);
      std::size_t cyclic = 0, paired = 0;
      for (std::size_t q = 0; q < s.facet_patch.size(); ++q){
        if (s.facet_patch[q] < 0 || c.patches[std::size_t(s.facet_patch[q])].type != "cyclic") continue;
        ++cyclic;
        paired += s.facet_partner[q] >= 0;
      }
      REQUIRE(cyclic == 2*2*(5*3 + 4*3));
      REQUIRE(paired == cyclic);
    }
  }
}
