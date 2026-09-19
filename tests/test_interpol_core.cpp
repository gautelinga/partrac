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
#include "cell_locate.hpp"
#include "dolfin_helpers.hpp"
#include "dolfin_spaces.hpp"
#include "p12_eval.hpp"
#include "PeriodicBC.hpp"
#include "MeshInterpol.hpp"
#include "StructuredInterpol.hpp"

namespace {

// The walk, then the tree, as MeshInterpol::locate
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
    period = periodic_lengths(periodic, lo, hi, dim_of<Cell>);
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

// MeshInterpol's cell tables on a unit mesh, without fields or files
template<typename Cell>
struct BareMesh : public MeshInterpol<Cell> {
  using MeshInterpol<Cell>::dolfin_cells_;
  BareMesh(std::shared_ptr<dolfin::Mesh> m, const std::vector<bool>& periodic) : MeshInterpol<Cell>("") {
    this->mesh = m;
    this->periodic = periodic;
    this->init_mesh_geometry();
    for (dolfin::CellIterator c(*m); !c.end(); ++c){
      this->dolfin_cells_.push_back(*c);
      this->cells_.push_back(Cell(*c));
    }
    this->build_facet_table();
  }
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
