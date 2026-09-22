// The spatial kernel of the divergence-free split field, on cells built here,
// so it is tested without a mesh, a file or dolfin.
//
// The fixture that needs no cleaner is a quadratic divergence-free polynomial:
// its P2 interpolant on a cell has exactly zero net flux, and on each sub-cell
// of the split it is itself a divergence-free P2 field, so uniqueness says the
// construction must give the polynomial back to round-off, value and gradient.
// Around that: the reference matrix against a solve on the physical cell (which
// is what the Piola map claims to replace), div u = 0 for arbitrary zero-flux
// data, continuity across the macro facets and the split's own facets, and no
// slip with u.n ~ delta^2 above a facet whose data is at rest.
#include <catch2/catch.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include <eigen3/Eigen/Dense>

#include "Tet.hpp"
#include "Triangle.hpp"
#include "p12_eval.hpp"
#include "split_eval.hpp"

namespace {

using split_eval::Split;

// A quadratic divergence-free velocity: the perpendicular gradient of a cubic
// stream function in 2D, the curl of a cubic potential in 3D. The coefficients
// are off every symmetry, so a wrong node order does not pass by luck.
Vector3d poly_u(const Vector3d& p, const int D){
  const double x = p[0], y = p[1], z = p[2];
  if (D == 2)
    return {-0.7*x*x + 1.8*x*y + 0.75*y*y - 0.3*x + 1.0*y - 0.2,
            -1.2*x*x + 1.4*x*y - 0.9*y*y - 1.2*x + 0.3*y - 0.8, 0.};
  return {-0.8*x*x + 0.75*y*y - 1.8*z*z + 1.6*x*y + 1.8*x*z + 0.4*x + 0.4,
           1.5*x*x - 1.5*y*y + 1.2*z*z + 1.2*x*y - 0.7*x + 0.4*y + 0.4*z,
          -0.9*y*y - 0.9*z*z + 0.4*x*z + 1.4*y*z + 1.6*x - 0.8*z - 0.5};
}

// grad(i, j) is du_i/dx_j
Matrix3d poly_grad(const Vector3d& p, const int D){
  const double x = p[0], y = p[1], z = p[2];
  Matrix3d g = Matrix3d::Zero();
  if (D == 2){
    g(0, 0) = -1.4*x + 1.8*y - 0.3;   g(0, 1) = 1.8*x + 1.5*y + 1.0;
    g(1, 0) = -2.4*x + 1.4*y - 1.2;   g(1, 1) = 1.4*x - 1.8*y + 0.3;
    return g;
  }
  g(0, 0) = -1.6*x + 1.6*y + 1.8*z + 0.4;  g(0, 1) = 1.6*x + 1.5*y;      g(0, 2) = 1.8*x - 3.6*z;
  g(1, 0) = 3.0*x + 1.2*y - 0.7;           g(1, 1) = 1.2*x - 3.0*y + 0.4; g(1, 2) = 2.4*z + 0.4;
  g(2, 0) = 0.4*z + 1.6;                   g(2, 1) = -1.8*y + 1.4*z;      g(2, 2) = 0.4*x + 1.4*y - 1.8*z - 0.8;
  return g;
}

// A cell of jittered vertices, and the same vertices as a flat array. Cells
// well away from degenerate: the local problem is unique whatever the shape,
// but a sliver's interior values are large and its round-off with them, which
// the checks below would read as an error of the construction.
template<typename Cell>
struct RandCell {
  static constexpr int D = Cell::n_verts - 1;
  std::array<Vector3d, Cell::n_verts> v;
  Cell cell;
  RandCell(std::mt19937& rng, const double jitter = 0.35){
    std::uniform_real_distribution<double> uni(-jitter, jitter);
    double x[Cell::n_verts][3] = {};
    for (int tries = 0; tries < 1000; ++tries){
      for (int k = 0; k < Cell::n_verts; ++k)
        for (int c = 0; c < D; ++c)
          x[k][c] = (k == c + 1 ? 1. : 0.) + uni(rng);
      Eigen::Matrix<double, D, D> J;
      for (int k = 0; k < D; ++k)
        for (int q = 0; q < D; ++q) J(q, k) = x[k + 1][q] - x[0][q];
      if (std::abs(J.determinant()) > 0.4) break;
    }
    for (int k = 0; k < Cell::n_verts; ++k) v[k] = {x[k][0], x[k][1], x[k][2]};
    if constexpr (D == 2) cell = Cell(x[0], x[1], x[2]);
    else                  cell = Cell(x[0], x[1], x[2], x[3]);
  }
  // The cell's affine Jacobian and its inverse, row-major
  void jacobian(double* J, double* Jinv) const {
    for (int k = 0; k < D; ++k){
      const Vector3d gl = cell.bary_grad(k + 1);
      for (int q = 0; q < D; ++q){
        J[q*D + k] = v[k + 1][q] - v[0][q];
        Jinv[k*D + q] = gl[q];
      }
    }
  }
  // The position of P2 slot k: a vertex, or the midpoint of the edge it holds
  Vector3d node(const std::size_t k) const {
    constexpr auto se = split_eval::slot_edges<Cell>();
    return 0.5*(v[se[k][0]] + v[se[k][1]]);
  }
};

// The interior values of one cell, by the construction under test
template<typename Cell>
std::vector<double> interior(const RandCell<Cell>& c,
                             const typename Split<Cell>::RefMatrix& R,
                             const std::vector<double>& g){
  constexpr int D = Split<Cell>::D;
  double J[D*D], Jinv[D*D];
  c.jacobian(J, Jinv);
  std::vector<double> u_int(Split<Cell>::n_int*std::size_t(D));
  split_eval::interior_values<Cell>(R, J, Jinv, g.data(), u_int.data());
  return u_int;
}

// The split field at a point of the cell: the sub-cell of the smallest
// barycentric, the P2 basis in mu, and the node values gathered through the
// tables. The node ids are the slots themselves here, so g serves as the field.
template<typename Cell>
Vector3d eval_at(const RandCell<Cell>& c, const std::vector<double>& g,
                 const std::vector<double>& u_int, const Vector3d& p, Matrix3d* grad = nullptr){
  constexpr int D = Split<Cell>::D;
  constexpr std::size_t N = Split<Cell>::n_sub;
  std::array<double, 4> bary{};
  c.cell.contains(p, bary);
  double mu[Cell::n_verts];
  const int i = split_eval::sub_cell<Cell::n_verts>(bary, mu);
  std::uint32_t row[Cell::n_dofs_max];
  for (std::size_t k = 0; k < Cell::n_dofs_max; ++k) row[k] = std::uint32_t(k);
  std::array<double, N*3> block;
  split_eval::gather_sub<Cell, D>(i, row, g.data(), u_int.data(), block.data());
  std::array<double, N> Nv;
  split_eval::sub_basis<Cell>(mu, Nv.data());
  if (grad){
    std::array<double, N> Nx, Ny, Nz;
    split_eval::sub_deriv<Cell>(c.cell, i, mu, Nx.data(), Ny.data(), Nz.data());
    *grad = block_gradient<D>(Nx.data(), Ny.data(), Nz.data(), block.data(), N);
  }
  return block_value<D>(Nv.data(), block.data(), N);
}

// A point drawn inside the cell, from barycentrics of an exponential draw
template<typename Cell>
Vector3d inside(const RandCell<Cell>& c, std::mt19937& rng){
  std::exponential_distribution<double> e(1.);
  double w[Cell::n_verts], s = 0.;
  for (int k = 0; k < Cell::n_verts; ++k){ w[k] = e(rng) + 1e-6; s += w[k]; }
  Vector3d p = Vector3d::Zero();
  for (int k = 0; k < Cell::n_verts; ++k) p += (w[k]/s)*c.v[k];
  return p;
}

// One cell's net flux, divided by the cell measure
template<typename Cell>
double net_flux(const Cell& cell, const std::vector<double>& g){
  double net = 0., big = 0.;
  split_eval::flux_parts<Cell>(cell, g.data(), 1., net, big);
  return net;
}

// Its net flux against its largest facet flux
template<typename Cell>
double flux_ratio(const Cell& cell, const std::vector<double>& g){
  double net = 0., big = 0.;
  split_eval::flux_parts<Cell>(cell, g.data(), 1., net, big);
  return big > 0. ? std::abs(net)/big : 0.;
}

// Move slot s until the cell's net flux vanishes; the flux is linear in it, so
// one solve of the D x 1 dependence does it
template<typename Cell>
void zero_flux(const Cell& cell, std::vector<double>& g, const std::size_t s){
  constexpr int D = Split<Cell>::D;
  const double f0 = net_flux<Cell>(cell, g);
  Vector3d w = Vector3d::Zero();
  for (int c = 0; c < D; ++c){
    g[s*D + std::size_t(c)] += 1.;
    w[c] = net_flux<Cell>(cell, g) - f0;
    g[s*D + std::size_t(c)] -= 1.;
  }
  const double n2 = w.squaredNorm();
  REQUIRE(n2 > 1e-12);
  for (int c = 0; c < D; ++c) g[s*D + std::size_t(c)] -= f0*w[c]/n2;
}

// The reference matrix's claim, solved on the physical cell instead: the P2
// field on its split with this boundary data and div u = 0 at the vertices of
// every sub-cell. Independent of the Piola map, so it is what checks it.
template<typename Cell>
std::vector<double> direct_interior(const RandCell<Cell>& c, const std::vector<double>& g){
  constexpr int D = Split<Cell>::D;
  constexpr int nv = Cell::n_verts;
  constexpr int nrows = nv*nv;
  constexpr int na = int(Split<Cell>::n_int)*D;
  constexpr int nb = int(Split<Cell>::n_bnd)*D;
  constexpr auto sn = split_eval::sub_nodes<Cell>();
  Vector3d z = Vector3d::Zero();
  for (int k = 0; k < nv; ++k) z += c.v[k]/double(nv);

  Eigen::Matrix<double, nrows, na> A = Eigen::Matrix<double, nrows, na>::Zero();
  Eigen::Matrix<double, nrows, nb> B = Eigen::Matrix<double, nrows, nb>::Zero();
  for (int o = 0; o < nv; ++o){
    double s[nv][3] = {};
    int m = 0;
    for (int j = 0; j < nv; ++j)
      if (j != o){ for (int q = 0; q < 3; ++q) s[m][q] = c.v[j][q]; ++m; }
    for (int q = 0; q < 3; ++q) s[nv - 1][q] = z[q];
    const Cell sub = [&]{
      if constexpr (nv == 3) return Cell(s[0], s[1], s[2]);
      else                   return Cell(s[0], s[1], s[2], s[3]);
    }();
    for (int q = 0; q < nv; ++q){
      double bary[4] = {0., 0., 0., 0.};
      bary[q] = 1.;
      std::array<double, Cell::n_dofs_max> Nx, Ny, Nz;
      if constexpr (nv == 3) sub.quadderiv(bary[0], bary[1], bary[2], Nx.data(), Ny.data());
      else                   sub.quadderiv(bary[0], bary[1], bary[2], bary[3],
                                           Nx.data(), Ny.data(), Nz.data());
      const double* dN[3] = {Nx.data(), Ny.data(), Nz.data()};
      for (std::size_t k = 0; k < Split<Cell>::n_sub; ++k){
        const int t = sn[o][k];
        for (int cc = 0; cc < D; ++cc){
          if (t >= 0) B(o*nv + q, t*D + cc) += dN[cc][k];
          else        A(o*nv + q, (-1 - t)*D + cc) += dN[cc][k];
        }
      }
    }
  }
  Eigen::Matrix<double, nb, 1> gv;
  for (int j = 0; j < nb; ++j) gv[j] = g[std::size_t(j)];
  const Eigen::Matrix<double, na, 1> x = A.colPivHouseholderQr().solve(-(B*gv));
  return std::vector<double>(x.data(), x.data() + na);
}

// The P2 interpolant of the polynomial on a cell, by slot
template<typename Cell>
std::vector<double> poly_nodes(const RandCell<Cell>& c){
  constexpr int D = Split<Cell>::D;
  std::vector<double> g(Split<Cell>::n_bnd*std::size_t(D));
  for (std::size_t k = 0; k < Split<Cell>::n_bnd; ++k){
    const Vector3d u = poly_u(c.node(k), D);
    for (int cc = 0; cc < D; ++cc) g[k*std::size_t(D) + std::size_t(cc)] = u[cc];
  }
  return g;
}

template<typename Cell>
std::vector<double> random_nodes(std::mt19937& rng){
  constexpr int D = Split<Cell>::D;
  std::uniform_real_distribution<double> uni(-1., 1.);
  std::vector<double> g(Split<Cell>::n_bnd*std::size_t(D));
  for (double& v : g) v = uni(rng);
  return g;
}


template<typename Cell>
void check_rank(){
  int rank = 0;
  const auto R = split_eval::reference_matrix<Cell>(rank);
  REQUIRE(rank == int(Split<Cell>::n_int)*Split<Cell>::D);
  REQUIRE(R.size() == Split<Cell>::n_int*std::size_t(Split<Cell>::D)
                     *Split<Cell>::n_bnd*std::size_t(Split<Cell>::D));
}

template<typename Cell>
void check_polynomial(){
  constexpr int D = Split<Cell>::D;
  int rank = 0;
  const auto R = split_eval::reference_matrix<Cell>(rank);
  std::mt19937 rng(11);
  for (int trial = 0; trial < 40; ++trial){
    const RandCell<Cell> c(rng);
    const auto g = poly_nodes(c);
    // Its P2 interpolant has no net flux, which is what the loader demands
    REQUIRE(flux_ratio<Cell>(c.cell, g) < 1e-12);
    const auto u_int = interior(c, R, g);
    // The interior nodes hold the polynomial too: the barycenter first, then
    // the midpoints of the segments from each vertex to it
    Vector3d z = Vector3d::Zero();
    for (int k = 0; k < Cell::n_verts; ++k) z += c.v[k]/double(Cell::n_verts);
    for (std::size_t n = 0; n < Split<Cell>::n_int; ++n){
      const Vector3d p = n == 0 ? z : 0.5*(c.v[n - 1] + z);
      const Vector3d w = poly_u(p, D);
      for (int cc = 0; cc < D; ++cc)
        REQUIRE(u_int[n*std::size_t(D) + std::size_t(cc)] == Approx(w[cc]).margin(1e-12));
    }
    for (int q = 0; q < 20; ++q){
      const Vector3d p = inside(c, rng);
      Matrix3d grad;
      const Vector3d u = eval_at(c, g, u_int, p, &grad);
      const Vector3d w = poly_u(p, D);
      const Matrix3d gw = poly_grad(p, D);
      for (int cc = 0; cc < D; ++cc){
        REQUIRE(u[cc] == Approx(w[cc]).margin(1e-11));
        for (int j = 0; j < D; ++j)
          REQUIRE(grad(cc, j) == Approx(gw(cc, j)).margin(1e-10));
      }
    }
  }
}

template<typename Cell>
void check_against_direct_solve(){
  int rank = 0;
  const auto R = split_eval::reference_matrix<Cell>(rank);
  std::mt19937 rng(23);
  for (int trial = 0; trial < 20; ++trial){
    const RandCell<Cell> c(rng);
    std::vector<double> g = random_nodes<Cell>(rng);
    zero_flux(c.cell, g, Split<Cell>::n_bnd - 1);
    const auto a = interior(c, R, g);
    const auto b = direct_interior(c, g);
    REQUIRE(a.size() == b.size());
    double scale = 1e-6;
    for (const double v : b) scale = std::max(scale, std::abs(v));
    for (std::size_t i = 0; i < a.size(); ++i)
      REQUIRE(a[i] == Approx(b[i]).margin(1e-9*scale));
  }
}

template<typename Cell>
void check_divergence_free(){
  constexpr int D = Split<Cell>::D;
  int rank = 0;
  const auto R = split_eval::reference_matrix<Cell>(rank);
  std::mt19937 rng(37);
  for (int trial = 0; trial < 30; ++trial){
    const RandCell<Cell> c(rng);
    std::vector<double> g = random_nodes<Cell>(rng);
    // Arbitrary data is made zero-flux by moving one midpoint, as the tool does
    zero_flux(c.cell, g, Split<Cell>::n_bnd - 1);
    REQUIRE(flux_ratio<Cell>(c.cell, g) < 1e-12);
    const auto u_int = interior(c, R, g);
    for (int q = 0; q < 30; ++q){
      Matrix3d grad;
      eval_at(c, g, u_int, inside(c, rng), &grad);
      double div = 0., size = 0.;
      for (int cc = 0; cc < D; ++cc){
        div += grad(cc, cc);
        for (int j = 0; j < D; ++j) size = std::max(size, std::abs(grad(cc, j)));
      }
      REQUIRE(std::abs(div) < 1e-11*std::max(size, 1.));
    }
  }
}

template<typename Cell>
void check_internal_facets(){
  constexpr int D = Split<Cell>::D;
  int rank = 0;
  const auto R = split_eval::reference_matrix<Cell>(rank);
  std::mt19937 rng(41);
  std::uniform_real_distribution<double> uni(0.05, 0.9);
  for (int trial = 0; trial < 20; ++trial){
    const RandCell<Cell> c(rng);
    std::vector<double> g = random_nodes<Cell>(rng);
    zero_flux(c.cell, g, Split<Cell>::n_bnd - 1);
    const auto u_int = interior(c, R, g);
    // Where two barycentrics are the smallest and equal the point lies on an
    // internal facet of the split; step a hair to either side of it
    for (int a = 0; a < Cell::n_verts; ++a){
      for (int b = a + 1; b < Cell::n_verts; ++b){
        const double s = 0.02 + 0.03*uni(rng);
        const double rest = (1. - 2.*s)/double(Cell::n_verts - 2);
        std::array<double, 4> bary{};
        for (int k = 0; k < Cell::n_verts; ++k)
          bary[k] = (k == a || k == b) ? s : rest;
        Vector3d p = Vector3d::Zero();
        for (int k = 0; k < Cell::n_verts; ++k) p += bary[k]*c.v[k];
        const Vector3d d = 1e-10*(c.v[a] - c.v[b]);
        const Vector3d u1 = eval_at(c, g, u_int, p + d);
        const Vector3d u2 = eval_at(c, g, u_int, p - d);
        for (int cc = 0; cc < D; ++cc)
          REQUIRE(u1[cc] == Approx(u2[cc]).margin(1e-8));
      }
    }
  }
}

template<typename Cell>
void check_macro_facet(){
  constexpr int D = Split<Cell>::D;
  constexpr int nv = Cell::n_verts;
  constexpr auto se = split_eval::slot_edges<Cell>();
  int rank = 0;
  const auto R = split_eval::reference_matrix<Cell>(rank);
  std::mt19937 rng(53);
  // Two cells sharing the facet opposite vertex 0: the second's apex is the
  // first's reflected through the facet's centroid
  for (int trial = 0; trial < 15; ++trial){
    RandCell<Cell> a(rng), b(rng);
    for (int k = 1; k < nv; ++k) b.v[k] = a.v[k];
    Vector3d mid = Vector3d::Zero();
    for (int k = 1; k < nv; ++k) mid += a.v[k]/double(nv - 1);
    b.v[0] = 2.*mid - a.v[0];
    double x[nv][3] = {};
    for (int k = 0; k < nv; ++k) for (int q = 0; q < 3; ++q) x[k][q] = b.v[k][q];
    if constexpr (nv == 3) b.cell = Cell(x[0], x[1], x[2]);
    else                   b.cell = Cell(x[0], x[1], x[2], x[3]);

    std::vector<double> ga = random_nodes<Cell>(rng), gb = random_nodes<Cell>(rng);
    // The P2 nodes of the shared facet are the same nodes, so the same values
    for (std::size_t k = 0; k < Split<Cell>::n_bnd; ++k)
      if (se[k][0] != 0 && se[k][1] != 0)
        for (int cc = 0; cc < D; ++cc)
          gb[k*std::size_t(D) + std::size_t(cc)] = ga[k*std::size_t(D) + std::size_t(cc)];
    // Each cell made zero-flux through a midpoint of its own, off the facet
    std::size_t off = 0;
    for (std::size_t k = std::size_t(nv); k < Split<Cell>::n_bnd; ++k)
      if (se[k][0] == 0 || se[k][1] == 0){ off = k; break; }
    zero_flux(a.cell, ga, off);
    zero_flux(b.cell, gb, off);
    const auto ia = interior(a, R, ga), ib = interior(b, R, gb);

    std::exponential_distribution<double> e(1.);
    for (int q = 0; q < 20; ++q){
      double w[nv] = {}, s = 0.;
      for (int k = 1; k < nv; ++k){ w[k] = e(rng) + 1e-6; s += w[k]; }
      Vector3d p = Vector3d::Zero();
      for (int k = 1; k < nv; ++k) p += (w[k]/s)*a.v[k];
      const Vector3d ua = eval_at(a, ga, ia, p);
      const Vector3d ub = eval_at(b, gb, ib, p);
      for (int cc = 0; cc < D; ++cc)
        REQUIRE(ua[cc] == Approx(ub[cc]).margin(1e-9));
    }
  }
}

template<typename Cell>
void check_no_slip(){
  constexpr int D = Split<Cell>::D;
  constexpr int nv = Cell::n_verts;
  constexpr auto se = split_eval::slot_edges<Cell>();
  int rank = 0;
  const auto R = split_eval::reference_matrix<Cell>(rank);
  std::mt19937 rng(67);
  for (int trial = 0; trial < 15; ++trial){
    const RandCell<Cell> c(rng);
    // The facet opposite vertex 0 is the wall: its every P2 node is at rest
    std::vector<double> g = random_nodes<Cell>(rng);
    for (std::size_t k = 0; k < Split<Cell>::n_bnd; ++k)
      if (se[k][0] != 0 && se[k][1] != 0)
        for (int cc = 0; cc < D; ++cc) g[k*std::size_t(D) + std::size_t(cc)] = 0.;
    std::size_t off = 0;
    for (std::size_t k = std::size_t(nv); k < Split<Cell>::n_bnd; ++k)
      if (se[k][0] == 0 || se[k][1] == 0){ off = k; break; }
    zero_flux(c.cell, g, off);
    const auto u_int = interior(c, R, g);

    const Vector3d n = c.cell.bary_grad(0).normalized();
    std::exponential_distribution<double> e(1.);
    for (int q = 0; q < 10; ++q){
      // Well inside the facet: nearer its boundary than delta the point above
      // it lies in another sub-cell, where the quadratic is another one
      double w[nv] = {}, s = 0.;
      for (int k = 1; k < nv; ++k){ w[k] = e(rng) + 1e-6; s += w[k]; }
      Vector3d foot = Vector3d::Zero();
      for (int k = 1; k < nv; ++k) foot += (0.5/double(nv - 1) + 0.5*w[k]/s)*c.v[k];
      // On the facet the field is at rest, up to the round-off of the point's
      // own barycentric coordinate
      const Vector3d u0 = eval_at(c, g, u_int, foot);
      for (int cc = 0; cc < D; ++cc) REQUIRE(std::abs(u0[cc]) < 1e-10);
      // Above it the normal component falls a hundredfold a decade of delta;
      // the foot's own value is subtracted, since it is zero only to round-off
      const double u0n = n.dot(u0);
      double un[2];
      for (int j = 0; j < 2; ++j){
        const double delta = (j == 0 ? 1e-4 : 1e-3);
        un[j] = std::abs(n.dot(eval_at(c, g, u_int, foot + delta*n)) - u0n);
      }
      CAPTURE(trial, q, un[0], un[1], u0n, flux_ratio<Cell>(c.cell, g));
      REQUIRE(un[1]/std::max(un[0], 1e-300) == Approx(100.).epsilon(0.02));
    }
  }
}

}  // namespace

TEST_CASE("The split's reference matrix pins every interior dof", "[split]"){
  check_rank<Triangle>();
  check_rank<Tet>();
}

TEST_CASE("The split field returns a quadratic divergence-free polynomial", "[split]"){
  check_polynomial<Triangle>();
  check_polynomial<Tet>();
}

TEST_CASE("The reference matrix and a solve on the cell itself agree", "[split]"){
  check_against_direct_solve<Triangle>();
  check_against_direct_solve<Tet>();
}

TEST_CASE("The split field is divergence-free for any zero-flux data", "[split]"){
  check_divergence_free<Triangle>();
  check_divergence_free<Tet>();
}

TEST_CASE("The split field is continuous across the split's own facets", "[split]"){
  check_internal_facets<Triangle>();
  check_internal_facets<Tet>();
}

TEST_CASE("The split fields of two cells agree on their shared facet", "[split]"){
  check_macro_facet<Triangle>();
  check_macro_facet<Tet>();
}

TEST_CASE("No slip on a facet at rest, with u.n growing as delta squared", "[split]"){
  check_no_slip<Triangle>();
  check_no_slip<Tet>();
}
