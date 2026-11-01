// The wall correction (wall_split.hpp) on single cells built here, without
// dolfin: random P2 blocks at rest on a cell's wall entities, the wall terms
// with random normals per entity, the correction from its pattern's table.
// The correction vanishes on the cell's faces, is continuous across the
// split's planes and its gradient is its derivative, the values alone those
// with the gradients; with it div u vanishes on the wall facets and on the
// wall edges (3D) or vertices (2D) from every sub-cell holding them, n.grad
// div u on the wall facets; the field does not depend on how the cell's
// vertices are numbered; with the wall terms left out the table is the least
// reference H1 field under the rows alone; and a pattern's table factored by
// its rank, as each sub-cell's polynomial, holds the dense one's nodal field.
#include <catch2/catch.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

#include "Tet.hpp"
#include "Triangle.hpp"
#include "wall_split.hpp"

namespace {

template<typename Cell>
using Coords = std::array<std::array<double, 3>, Cell::n_verts>;

template<typename Cell>
Cell cell_of(const Coords<Cell>& x){
  if constexpr (Cell::n_verts == 3) return Cell(x[0].data(), x[1].data(), x[2].data());
  else                              return Cell(x[0].data(), x[1].data(), x[2].data(), x[3].data());
}

template<typename Cell>
Coords<Cell> coords(){
  if constexpr (Cell::n_verts == 3) return {{{0.1, -0.2, 0.}, {1.3, 0.15, 0.}, {0.35, 0.9, 0.}}};
  else return {{{0.1, -0.2, 0.05}, {1.2, 0.1, -0.1}, {0.3, 1.1, 0.2}, {0.25, 0.3, 0.95}}};
}

// Rest bits closed: a facet's edges and vertices, an edge's vertices (3D)
template<int NV>
unsigned closed(unsigned rest){
  using S = wall_split::Shape<NV>;
  constexpr auto ev = wall_split::entity_verts<NV>();
  for (int k = 0; k < NV; ++k){
    if (!(rest >> k & 1)) continue;
    for (std::size_t q = 0; q < ev.size(); ++q)
      if (ev[q][0] != k && ev[q][1] != k) rest |= 1u << (NV + int(q));
  }
  if constexpr (NV == 4)
    for (std::size_t q = 0; q < ev.size(); ++q)
      if (rest >> (NV + int(q)) & 1) rest |= (1u << (S::n_bits + ev[q][0])) | (1u << (S::n_bits + ev[q][1]));
  return rest;
}

// A random unit vector fixed by a key
Vector3d keyed(const unsigned key, const int D){
  std::mt19937 g(key*2654435761u + 17u);
  std::normal_distribution<double> nd;
  Vector3d n(nd(g), nd(g), D == 3 ? nd(g) : 0.);
  return n.normalized();
}

// The wall terms of a cell at rest on rest, by bit, and their normals, keyed by
// the labels of the vertices they join: edge ends from a vertex at rest but
// along an edge at rest, side faces of an edge at rest but wall facets, and
// in 3D the bubble's facets
template<int NV>
std::uint32_t wall_terms(const unsigned rest, const std::array<int, NV>& label,
                         std::vector<double>& normals){
  constexpr int D = NV - 1;
  using S = wall_split::Shape<NV>;
  constexpr auto edges = split_eval::lex_edges<NV>();
  constexpr auto apex = wall_split::edge_apices<NV>();
  const auto vertex = [&](const int v){ return (rest >> (NV == 4 ? S::n_bits + v : NV + v) & 1) != 0; };
  const auto edge = [&](const int e){
    return (rest >> (NV == 4 ? NV + e : 3 - edges[std::size_t(e)][0] - edges[std::size_t(e)][1]) & 1) != 0;
  };
  const auto push = [&](const Vector3d& n){ for (int d = 0; d < D; ++d) normals.push_back(n[d]); };
  std::uint32_t terms = 0;
  for (int e = 0; e < S::n_edges; ++e)
    for (int o = 0; o < 2; ++o){
      const int w = edges[std::size_t(e)][std::size_t(o)], v = edges[std::size_t(e)][std::size_t(1 - o)];
      if (!vertex(w) || edge(e)) continue;
      terms |= 1u << (2*e + o);
      push(keyed(unsigned(64*label[std::size_t(w)] + label[std::size_t(v)]), D));
    }
  for (int e = 0; e < S::n_edges; ++e)
    for (int j = 0; j < NV - 2; ++j){
      const int a = edges[std::size_t(e)][0], b = edges[std::size_t(e)][1], v = apex[std::size_t(e)][std::size_t(j)];
      if (!edge(e) || (NV == 4 && (rest >> (6 - a - b - v) & 1))) continue;
      terms |= 1u << (S::triple0 + e*(NV - 2) + j);
      const int la = std::min(label[std::size_t(a)], label[std::size_t(b)]);
      const int lb = std::max(label[std::size_t(a)], label[std::size_t(b)]);
      push(keyed(unsigned(4096 + 64*64*la + 64*lb + label[std::size_t(v)]), D));
    }
  if (NV == 4)
    for (int k = 0; k < NV; ++k){
      if (!(rest >> k & 1)) continue;
      terms |= 1u << (S::bubble0 + k);
      push(keyed(unsigned(900000 + label[std::size_t(k)]), D));
      push(0.7*keyed(unsigned(910000 + label[std::size_t(k)]), D));
    }
  return terms;
}

// A cell's P2 field plus the wall terms and the correction, from two blocks u
// and v at rest on the rest bits' entities, made ready as a loader does
template<typename Cell>
struct Corrected {
  static constexpr int nv = Cell::n_verts, D = nv - 1;
  static constexpr std::size_t n2 = Cell::n_dofs_max;
  using S = wall_split::Shape<nv>;
  Cell cell;
  wall_split::Frame<Cell> frame;
  unsigned rest;
  std::array<int, nv> label;
  int r = 0;
  std::vector<std::int16_t> idx;
  std::vector<double> W, normals, F;
  std::uint32_t terms;
  std::array<double, D*n2> u{}, v{};
  double du[nv], dv[nv];
  wall_split::Tables<nv> tables;
  wall_split::Ready<Cell> rdy;
  bool terms_on = true;

  Corrected(const Cell& c, const unsigned rest, std::mt19937& rng,
            const std::array<int, nv>& label = identity())
    : cell(c), frame(c), rest(rest), label(label),
      W(wall_split::pattern_tables<nv>({std::uint16_t(rest & S::pattern_mask)})),
      F(wall_split::factor<nv>(W.data(), r, idx)), terms(wall_terms<nv>(rest, label, normals)){
    std::uniform_real_distribution<double> uni(-1., 1.);
    static constexpr auto slots = wall_split::entity_slots<Cell>();
    unsigned zero = 0;
    for (unsigned m = rest; m; m &= m - 1) zero |= slots[std::size_t(__builtin_ctz(m))];
    for (int d = 0; d < D; ++d)
      for (std::size_t q = 0; q < n2; ++q){
        const bool z = zero >> q & 1;
        u[std::size_t(d)*n2 + q] = z ? 0. : uni(rng);
        v[std::size_t(d)*n2 + q] = z ? 0. : uni(rng);
      }
    tables.cell = {wall_split::WallCell{0, terms, 0, std::uint16_t(rest)}};
    tables.normals = normals;
    tables.pattern = {std::uint16_t(rest & S::pattern_mask)};
    tables.W = F;
    tables.at = {0};
    tables.rank = {r};
    tables.idx = idx;
    update();
  }
  static std::array<int, nv> identity(){
    std::array<int, nv> l;
    std::iota(l.begin(), l.end(), 0);
    return l;
  }
  void update(){
    wall_split::vertex_divergence<Cell>(frame, u.data(), du);
    wall_split::vertex_divergence<Cell>(frame, v.data(), dv);
    tables.cell[0].terms = terms_on ? terms : 0u;
    wall_split::ready<Cell>(cell, tables, tables.cell[0], u.data(), v.data(), rdy);
  }
  wall_split::Table table() const { return {F.data(), r, idx.data()}; }
  // P2 alone at b
  void p2(const std::array<double, 4>& b, Vector3d& U, Matrix3d& J, const bool first) const {
    std::array<double, n2> N, Nx, Ny, Nz{};
    if constexpr (nv == 3){
      cell.quadbasis(b[0], b[1], b[2], N.data());
      cell.quadderiv(b[0], b[1], b[2], Nx.data(), Ny.data());
    }
    else {
      cell.quadbasis(b[0], b[1], b[2], b[3], N.data());
      cell.quadderiv(b[0], b[1], b[2], b[3], Nx.data(), Ny.data(), Nz.data());
    }
    const double* dN[3] = {Nx.data(), Ny.data(), Nz.data()};
    const double* w = first ? u.data() : v.data();
    U.setZero();
    J.setZero();
    for (int d = 0; d < D; ++d)
      for (std::size_t q = 0; q < n2; ++q){
        U[d] += N[q]*w[std::size_t(d)*n2 + q];
        for (int x = 0; x < D; ++x) J(d, x) += dN[x][q]*w[std::size_t(d)*n2 + q];
      }
  }
  // The whole field, or the correction alone (whole false), at b in sub-cell
  // j (-1: the least barycentric's)
  void at(const std::array<double, 4>& b, Vector3d& U, Matrix3d& J, const bool first = true,
          const bool whole = true, const int j = -1) const {
    Vector3d Ua = Vector3d::Zero(), Ub = Vector3d::Zero();
    Matrix3d Ja = Matrix3d::Zero(), Jb = Matrix3d::Zero();
    double mu[nv];
    const int s = j < 0 ? split_eval::sub_cell<nv>(b, mu) : j;
    if (whole) wall_split::add_ready_in<Cell, true>(tables, tables.cell[0], rdy, s, b, Ua, Ub, Ja, Jb);
    else wall_split::add_in<Cell, true>(frame, table(), s, b, du, dv, rdy.terms, Ua, Ub, Ja, Jb);
    U = first ? Ua : Ub;
    J = first ? Ja : Jb;
  }
  // The values alone of the whole field at b, both fields
  void values(const std::array<double, 4>& b, Vector3d& U, Vector3d& V) const {
    U = V = Vector3d::Zero();
    Matrix3d unused;
    wall_split::add_ready<Cell, false>(tables, tables.cell[0], rdy, b, U, V, unused, unused);
  }
  double div(const std::array<double, 4>& b, const bool first, const int j = -1) const {
    Vector3d U;
    Matrix3d J;
    at(b, U, J, first, true, j);
    return J.trace();
  }
  // The barycentrics of the point t along the inward normal of facet k from b
  std::array<double, 4> along(const std::array<double, 4>& b, const int k, const double t) const {
    const Vector3d n = frame.G[k].normalized();
    std::array<double, 4> o = b;
    for (int i = 0; i < nv; ++i) o[std::size_t(i)] += t*frame.G[i].dot(n);
    return o;
  }
};

template<typename M>
double largest(const Eigen::MatrixBase<M>& m){ return m.cwiseAbs().maxCoeff(); }

// Random barycentrics, each at least lo
template<int NV>
std::array<double, 4> random_point(std::mt19937& rng, const double lo = 0.05){
  std::uniform_real_distribution<double> uni(0., 1.);
  std::array<double, 4> b{};
  double s = 0.;
  for (int i = 0; i < NV; ++i) s += b[std::size_t(i)] = uni(rng);
  for (int i = 0; i < NV; ++i) b[std::size_t(i)] = lo + (1. - double(NV)*lo)*b[std::size_t(i)]/s;
  return b;
}

// Random barycentrics on facet k, each other one at least lo
template<int NV>
std::array<double, 4> on_facet(std::mt19937& rng, const int k, const double lo = 0.15){
  std::uniform_real_distribution<double> uni(0., 1.);
  std::array<double, 4> b{};
  double s = 0.;
  for (int i = 0; i < NV; ++i)
    if (i != k) s += b[std::size_t(i)] = uni(rng);
  for (int i = 0; i < NV; ++i)
    if (i != k) b[std::size_t(i)] = lo + (1. - double(NV - 1)*lo)*b[std::size_t(i)]/s;
  return b;
}

// Zero on the faces, the same from both sides of a split plane, the gradient
// the value's derivative by central differences
template<typename Cell>
void shape(const unsigned rest, const int seed){
  constexpr int nv = Cell::n_verts, D = nv - 1;
  std::mt19937 rng{std::uint32_t(seed)};
  const Corrected<Cell> f(cell_of<Cell>(coords<Cell>()), rest, rng);
  double gscale = 0., vscale = 0.;
  for (int p = 0; p < 40; ++p){
    const std::array<double, 4> b = random_point<nv>(rng);
    Vector3d U;
    Matrix3d J;
    f.at(b, U, J, p % 2, false);
    gscale = std::max(gscale, largest(J));
    vscale = std::max(vscale, largest(U));
  }
  INFO("rest " << rest << ", scale " << vscale << ", gradient scale " << gscale);
  REQUIRE(gscale > 0.);
  for (int p = 0; p < 40; ++p){
    const bool first = p % 2;
    const std::array<double, 4> b = random_point<nv>(rng);
    Vector3d U;
    Matrix3d J;
    f.at(b, U, J, first, false);
    // A face
    for (int k = 0; k < nv; ++k){
      std::array<double, 4> c = b;
      const double r = 1. - c[std::size_t(k)];
      for (int i = 0; i < nv; ++i) c[std::size_t(i)] = i == k ? 0. : c[std::size_t(i)]/r;
      Vector3d Uf;
      Matrix3d Jf;
      f.at(c, Uf, Jf, first, false);
      CHECK(largest(Uf) == 0.);
    }
    // A plane: the two least barycentrics tied, from either sub-cell
    int j = 0;
    for (int i = 1; i < nv; ++i) if (b[std::size_t(i)] < b[std::size_t(j)]) j = i;
    const int a = (j + 1) % nv;
    std::array<double, 4> t = b;
    t[std::size_t(a)] = t[std::size_t(j)] = 0.5*(b[std::size_t(a)] + b[std::size_t(j)]);
    Vector3d Uj, Ua;
    Matrix3d Jj, Ja;
    f.at(t, Uj, Jj, first, false, j);
    f.at(t, Ua, Ja, first, false, a);
    CHECK(largest(Uj - Ua) <= 1e-13*vscale);
    // The gradient, of the whole field, and the values alone
    f.at(b, U, J, first);
    {
      Vector3d A, B;
      f.values(b, A, B);
      CHECK(largest(Vector3d((first ? A : B) - U)) <= 1e-13*std::max(vscale, largest(U)));
    }
    for (int x = 0; x < D; ++x){
      Vector3d e = Vector3d::Zero();
      e[x] = 1e-6;
      std::array<double, 4> bp = b, bm = b;
      for (int i = 0; i < nv; ++i){
        const double dl = f.frame.G[i].dot(e);
        bp[std::size_t(i)] += dl;
        bm[std::size_t(i)] -= dl;
      }
      Vector3d Up, Um;
      Matrix3d Jp, Jm;
      f.at(bp, Up, Jp, first);
      f.at(bm, Um, Jm, first);
      CHECK(largest(Vector3d((Up - Um)/2e-6 - J.col(x))) <= 1e-6*std::max(gscale, largest(J)));
    }
  }
}

// div u on the pattern's facets (in their sub-cell) and codimension-2
// entities (in every sub-cell holding them), against P2's own; n.grad div u
// on the facets, from div u at three heights in the facet's sub-cell
template<typename Cell>
void wall_divergence(const unsigned rest, const int seed){
  constexpr int nv = Cell::n_verts;
  constexpr auto ev = wall_split::entity_verts<nv>();
  std::mt19937 rng{std::uint32_t(seed)};
  std::uniform_real_distribution<double> uni(0.1, 0.9);
  const Corrected<Cell> f(cell_of<Cell>(coords<Cell>()), rest, rng);
  REQUIRE(f.terms != 0u);
  std::size_t rows = 0;
  double div = 0., div_p2 = 0., slope = 0., scale = 0.;
  const auto point = [&](const std::array<double, 4>& b, const int j){
    for (const bool first : {true, false}){
      Vector3d U;
      Matrix3d J, J2;
      f.at(b, U, J, first, true, j);
      f.p2(b, U, J2, first);
      scale = std::max(scale, largest(J));
      div = std::max(div, std::abs(J.trace()));
      div_p2 = std::max(div_p2, std::abs(J2.trace()));
      ++rows;
    }
  };
  for (int k = 0; k < nv; ++k){
    if (!(rest >> k & 1)) continue;
    for (int p = 0; p < 6; ++p){
      const std::array<double, 4> b = on_facet<nv>(rng, k);
      point(b, k);
      // div at heights h, 2h, 4h: its first-order term, (8 q1 - 6 q2 + q4)/3 with q = div/height
      const double gk = f.frame.G[k].norm(), h = 1e-2/gk;
      for (const bool first : {true, false}){
        double q[3];
        for (int s = 0; s < 3; ++s){
          const double t = h*double(1 << s);
          q[s] = f.div(f.along(b, k, t), first, k)/t;
        }
        slope = std::max(slope, std::abs(8.*q[0] - 6.*q[1] + q[2])/(3.*gk));
      }
    }
  }
  for (std::size_t e = 0; e < ev.size(); ++e){
    if (!(rest >> (nv + int(e)) & 1)) continue;
    const int a = ev[e][0], c = ev[e][1];
    for (int p = 0; p < (a == c ? 1 : 4); ++p){
      std::array<double, 4> b{};
      const double t = a == c ? 0. : uni(rng);
      b[std::size_t(a)] += 1. - t;
      b[std::size_t(c)] += t;
      for (int j = 0; j < nv; ++j)
        if (j != a && j != c) point(b, j);
    }
  }
  INFO("rest " << rest << ", rows " << rows << ", |grad u| " << scale << ", div " << div
       << ", P2's " << div_p2 << ", n.grad div / |grad l_k| " << slope);
  REQUIRE(rows > 0);
  CHECK(div <= 1e-12*scale);
  CHECK(div_p2 > 1e-3*scale);
  CHECK(slope <= 1e-8*scale);
}

// The cell with its vertices renumbered by every permutation, the data, the
// rest bits and the wall terms' normals with them: the same field at the same
// points
template<typename Cell>
void numbering(const unsigned rest, const int seed){
  constexpr int nv = Cell::n_verts, D = nv - 1;
  constexpr std::size_t n2 = Cell::n_dofs_max;
  constexpr auto ev = wall_split::entity_verts<nv>();
  constexpr auto edges = split_eval::lex_edges<nv>();
  using S = wall_split::Shape<nv>;
  std::mt19937 rng{std::uint32_t(seed)};
  const Coords<Cell> X = coords<Cell>();
  const Corrected<Cell> ref(cell_of<Cell>(X), rest, rng);
  std::vector<std::array<double, 4>> pts;
  for (int p = 0; p < 20; ++p) pts.push_back(random_point<nv>(rng));
  std::array<int, nv> perm;
  std::iota(perm.begin(), perm.end(), 0);
  // The P2 slot of the edge between two vertices
  const auto slot = [&](const int a, const int b){
    for (std::size_t e = 0; e < edges.size(); ++e)
      if ((edges[e][0] == a && edges[e][1] == b) || (edges[e][0] == b && edges[e][1] == a))
        return std::size_t(Cell::mid_[e]);
    return std::size_t(0);
  };
  int count = 0;
  double gap = 0., scale = 0.;
  do {
    // local vertex i of the new cell is perm[i] of the old
    Coords<Cell> Y;
    for (int i = 0; i < nv; ++i) Y[std::size_t(i)] = X[std::size_t(perm[std::size_t(i)])];
    unsigned r = 0;
    for (int k = 0; k < nv; ++k){
      if (rest >> perm[std::size_t(k)] & 1) r |= 1u << k;
      if (nv == 4 && (rest >> (S::n_bits + perm[std::size_t(k)]) & 1)) r |= 1u << (S::n_bits + k);
    }
    for (std::size_t e = 0; e < ev.size(); ++e){
      const int a = perm[std::size_t(ev[e][0])], b = perm[std::size_t(ev[e][1])];
      for (std::size_t q = 0; q < ev.size(); ++q)
        if ((ev[q][0] == a && ev[q][1] == b) || (ev[q][0] == b && ev[q][1] == a))
          if (rest >> (nv + int(q)) & 1) r |= 1u << (nv + int(e));
    }
    Corrected<Cell> f(cell_of<Cell>(Y), r, rng, perm);
    for (int d = 0; d < D; ++d){
      const std::size_t o = std::size_t(d)*n2;
      for (int i = 0; i < nv; ++i){
        f.u[o + std::size_t(i)] = ref.u[o + std::size_t(perm[std::size_t(i)])];
        f.v[o + std::size_t(i)] = ref.v[o + std::size_t(perm[std::size_t(i)])];
      }
      for (std::size_t e = 0; e < edges.size(); ++e){
        const std::size_t s = slot(perm[std::size_t(edges[e][0])], perm[std::size_t(edges[e][1])]);
        f.u[o + std::size_t(Cell::mid_[e])] = ref.u[o + s];
        f.v[o + std::size_t(Cell::mid_[e])] = ref.v[o + s];
      }
    }
    f.update();
    for (const auto& b : pts){
      std::array<double, 4> c{};
      for (int i = 0; i < nv; ++i) c[std::size_t(i)] = b[std::size_t(perm[std::size_t(i)])];
      for (const bool first : {true, false}){
        Vector3d U, V;
        Matrix3d J, K;
        ref.at(b, U, J, first);
        f.at(c, V, K, first);
        scale = std::max(scale, largest(J));
        gap = std::max({gap, largest(U - V)/largest(J), largest(J - K)});
      }
    }
    ++count;
  } while (std::next_permutation(perm.begin(), perm.end()));
  INFO("rest " << rest << ", permutations " << count << ", gap " << gap << ", scale "
       << scale);
  REQUIRE(scale > 0.);
  CHECK(gap <= 1e-12*scale);
}

// With the wall terms left out, the correction against the
// least reference H1 field under the pattern's rows, solved here through the
// Schur complement of the rows: the tables' vertex divergence parts and the
// fields at random points
template<typename Cell>
void without_wall_terms(const unsigned rest, const int seed){
  constexpr int nv = Cell::n_verts, D = nv - 1;
  using S = wall_split::Shape<nv>;
  using Mat = Eigen::MatrixXd;
  const wall_split::Reference ref = wall_split::reference<nv>();
  const int ni = ref.n_int;
  const unsigned pattern = rest & S::pattern_mask;
  std::vector<int> act;
  for (int r = 0; r < int(ref.bit.size()); ++r)
    if (pattern >> ref.bit[std::size_t(r)] & 1) act.push_back(r);
  const int na = int(act.size());
  Mat B(na, ni*D), T(na, nv);
  for (int r = 0; r < na; ++r){
    B.row(r) = ref.B.row(act[std::size_t(r)]).head(ni*D);
    T.row(r) = ref.T.row(act[std::size_t(r)]);
  }
  Mat Mi = Mat::Zero(ni*D, ni*D);
  const Mat Si = ref.S.topLeftCorner(ni, ni).inverse();
  for (int a = 0; a < ni; ++a)
    for (int c = 0; c < ni; ++c)
      for (int x = 0; x < D; ++x) Mi(a*D + x, c*D + x) = Si(a, c);
  const Mat Sch = B*Mi*B.transpose();
  Eigen::JacobiSVD<Mat> svd(Sch, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd sv = svd.singularValues();
  const double cut = 1e-11*sv[0];
  for (int i = 0; i < sv.size(); ++i) sv[i] = sv[i] > cut ? 1./sv[i] : 0.;
  const Mat W = Mi*B.transpose()*(svd.matrixV()*sv.asDiagonal()*svd.matrixU().transpose())*T;
  const Mat Ehi = ref.Eh.inverse();
  double gap = 0., scale = 0.;
  std::mt19937 rng{std::uint32_t(seed)};
  Corrected<Cell> f(cell_of<Cell>(coords<Cell>()), rest, rng);
  for (int a = 0; a < ni; ++a){
    const Mat w = Ehi*W.middleRows(a*D, D);
    for (int x = 0; x < D; ++x)
      for (int i = 0; i < nv; ++i){
        gap = std::max(gap, std::abs(f.W[std::size_t(a)*S::node + std::size_t(x*nv + i)] - w(x, i)));
        scale = std::max(scale, std::abs(w(x, i)));
      }
  }
  // The fields: the solve's node values through the split's basis
  f.terms_on = false;
  f.update();
  std::vector<double> Wd(f.W.size(), 0.);
  for (int a = 0; a < ni; ++a){
    const Mat w = Ehi*W.middleRows(a*D, D);
    for (int x = 0; x < D; ++x)
      for (int i = 0; i < nv; ++i) Wd[std::size_t(a)*S::node + std::size_t(x*nv + i)] = w(x, i);
  }
  int rd = 0;
  std::vector<std::int16_t> idxd;
  const std::vector<double> Fd = wall_split::factor<nv>(Wd.data(), rd, idxd);
  for (int p = 0; p < 20; ++p){
    const std::array<double, 4> b = random_point<nv>(rng);
    Vector3d U1 = Vector3d::Zero(), V1 = Vector3d::Zero(), U2 = Vector3d::Zero(), V2 = Vector3d::Zero();
    Matrix3d J1 = Matrix3d::Zero(), K1 = Matrix3d::Zero(), J2 = Matrix3d::Zero(), K2 = Matrix3d::Zero();
    wall_split::add<Cell, true>(f.frame, f.table(), b, f.du, f.dv, f.rdy.terms, U1, V1, J1, K1);
    wall_split::add<Cell, true>(f.frame, {Fd.data(), rd, idxd.data()}, b, f.du, f.dv, f.rdy.terms, U2, V2, J2, K2);
    scale = std::max(scale, largest(J2));
    gap = std::max({gap, largest(U1 - U2), largest(V1 - V2), largest(J1 - J2), largest(K1 - K2)});
  }
  INFO("rest " << rest << ", rows " << na << ", gap " << gap << ", scale " << scale);
  REQUIRE(scale > 0.);
  CHECK(gap <= 1e-9*scale);
}

// The factored table against the pattern's dense one, per input column: in
// each sub-cell j at random points, lambda_j Q_j from the factored table
// against the dense table's interior node values through the split's degree-K
// Lagrange basis; its rank padded to a block
template<int NV>
void factored(const unsigned rest, const int seed){
  constexpr int D = NV - 1;
  using S = wall_split::Shape<NV>;
  constexpr wall_split::Degree<NV, 3> d3{};
  const std::vector<double> W = wall_split::pattern_tables<NV>({std::uint16_t(rest & S::pattern_mask)});
  int r = 0;
  std::vector<std::int16_t> idx;
  const std::vector<double> F = wall_split::factor<NV>(W.data(), r, idx);
  const wall_split::Table t{F.data(), r, idx.data()};
  REQUIRE(r % wall_split::block == 0);
  REQUIRE(F.size() == wall_split::table_size<NV>(r));
  constexpr int nc = NV + S::n_slots*D;
  std::mt19937 rng{std::uint32_t(seed)};
  double gap = 0., scale = 0.;
  for (int k = 0; k < nc; ++k){
    // The dense table's node values of input k
    std::vector<double> node(std::size_t(S::n_int*D));
    for (int a = 0; a < S::n_int; ++a)
      for (int x = 0; x < D; ++x){
        const double* w = W.data() + std::size_t(a)*S::node;
        node[std::size_t(a*D + x)] = k < NV ? w[x*NV + k]
                                            : w[S::n_d + std::size_t((k - NV)/D*D*D + x*D + (k - NV) % D)];
      }
    for (int j = 0; j < NV; ++j){
      // Q_j of input k: its reduced input is row k of the sub-cell's inputs table
      const double* z = wall_split::sub_inputs<NV>(t, j).p + std::size_t(k*r);
      wall_split::pair y[d3.n*D];
      wall_split::sub_coefficients<NV>(t, j, z, z, y);
      for (int p = 0; p < 10; ++p){
        std::array<double, 4> b = random_point<NV>(rng, 0.);
        double mu[NV];
        split_eval::sub_levels<NV>(j, b, mu);
        bool in = true;
        for (int m = 0; m < NV; ++m) in = in && mu[m] >= 0.;
        if (!in){ --p; continue; }
        for (int x = 0; x < D; ++x){
          double nodal = 0., poly = 0.;
          for (int l = 0; l < S::n_act; ++l){
            double phi = 1.;
            for (int i = 0; i < NV; ++i)
              for (int m = 0; m < S::lat.act[std::size_t(l)][std::size_t(i)]; ++m)
                phi *= (double(wall_split::K)*mu[i] - double(m))/double(m + 1);
            nodal += phi*node[std::size_t(S::in.idx[std::size_t(j)][std::size_t(l)]*D + x)];
          }
          for (int q = 0; q < d3.n; ++q){
            double mono = 1.;
            for (int i = 0; i < NV; ++i) mono *= std::pow(b[std::size_t(i)], d3.e[std::size_t(q)][std::size_t(i)]);
            poly += y[q*D + x][0]*mono;
          }
          poly *= b[std::size_t(j)];
          gap = std::max(gap, std::abs(poly - nodal));
          scale = std::max(scale, std::abs(nodal));
        }
      }
    }
  }
  INFO("rest " << rest << ", rank " << r << ", gap " << gap << ", scale " << scale);
  REQUIRE(scale > 0.);
  CHECK(r < S::n_int*D);
  CHECK(gap <= 1e-12*scale);
}

}  // namespace

// Rest bits before closure. Triangles: facets in bits 0-2, vertices in bits
// 3-5; tets: facets in bits 0-3, edges 01, 02, 03, 12, 13, 23 in bits 4-9,
// vertices in bits 10-13
const unsigned tri_facet = closed<3>(0b001u), tri_two = closed<3>(0b011u), tri_vertex = 0b001'000u,
               tri_vertices = 0b011'000u;
const unsigned tet_facet = closed<4>(0b0001u), tet_two = closed<4>(0b0011u), tet_edge = closed<4>(0b000001'0000u),
               tet_edges = closed<4>(0b100001'0000u), tet_apex = closed<4>(0b0001u) | 1u << 10,
               tet_vertex = 1u << 12;

TEST_CASE("The wall correction vanishes on the cell's faces, is continuous on the split's planes, "
          "and its gradient is its derivative", "[wall_split]"){
  SECTION("triangles"){
    shape<Triangle>(tri_facet, 1); shape<Triangle>(tri_two, 2);
    shape<Triangle>(tri_vertex, 3); shape<Triangle>(tri_vertices, 4);
  }
  SECTION("tets"){
    shape<Tet>(tet_facet, 5); shape<Tet>(tet_two, 6); shape<Tet>(tet_edge, 7);
    shape<Tet>(tet_edges, 8); shape<Tet>(tet_apex, 9);
  }
}

TEST_CASE("With the wall correction div u vanishes on the wall entities at rest, and n.grad div u on the "
          "wall facets", "[wall_split]"){
  SECTION("triangles"){
    wall_divergence<Triangle>(tri_facet, 11); wall_divergence<Triangle>(tri_two, 12);
    wall_divergence<Triangle>(tri_vertex, 13); wall_divergence<Triangle>(tri_vertices, 14);
  }
  SECTION("tets"){
    wall_divergence<Tet>(tet_facet, 15); wall_divergence<Tet>(tet_two, 16);
    wall_divergence<Tet>(tet_edge, 17); wall_divergence<Tet>(tet_edges, 18);
    wall_divergence<Tet>(tet_apex, 19);
  }
}

TEST_CASE("The wall correction does not depend on the numbering of the cell's vertices", "[wall_split]"){
  SECTION("triangles"){
    numbering<Triangle>(tri_facet, 21); numbering<Triangle>(tri_vertices, 22);
  }
  SECTION("tets"){
    numbering<Tet>(tet_two, 23); numbering<Tet>(tet_edges, 24); numbering<Tet>(tet_apex, 25);
  }
}

TEST_CASE("A pattern's factored table, as each sub-cell's polynomial, holds its dense one's nodal field",
          "[wall_split]"){
  SECTION("triangles"){
    factored<3>(tri_facet, 41); factored<3>(tri_two, 42); factored<3>(tri_vertex, 43);
    factored<3>(tri_vertices, 44);
  }
  SECTION("tets"){
    factored<4>(tet_facet, 45); factored<4>(tet_two, 46); factored<4>(tet_edge, 47); factored<4>(tet_edges, 48);
    factored<4>(tet_apex, 49);
  }
}

TEST_CASE("Without the wall terms the wall correction is the least reference H1 field under the rows",
          "[wall_split]"){
  SECTION("triangles"){
    without_wall_terms<Triangle>(tri_facet, 31); without_wall_terms<Triangle>(tri_two, 32);
    without_wall_terms<Triangle>(tri_vertices, 33);
  }
  SECTION("tets"){
    without_wall_terms<Tet>(tet_facet, 34); without_wall_terms<Tet>(tet_two, 35);
    without_wall_terms<Tet>(tet_edges, 36);
  }
}

// wall_upgrade on a P2 velocity, read by the stamped and the frequency loader
// on dolfin's unit square and cube with every side a wall, whose corner and
// edge cells have two wall facets: continuity across the faces of every cell
// with a wall vertex and across the planes of its split,
// div u = 0 on the walls, the correction only on the wall entities whose P2
// nodes are all at rest, and a later stamp that moves a wall edge or vertex at
// rest in the first refused.
#ifdef USE_DOLFIN
#include <catch2/catch.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "Error.hpp"
#include "SimplexInterpol.hpp"
#include "interpolators.hpp"
#include "held_eval.hpp"
#include "case_dir.hpp"
#include "region_cases.hpp"

namespace {

constexpr double T0 = 0.4;
const std::vector<bool> walls = {false, false, false};

// At rest on every side of the unit box, other at each stamp
template<int D>
Vector3d at_rest(const Vector3d& x, const int k){
  double b = 1.;
  for (int d = 0; d < D; ++d) b *= x[d]*(1. - x[d]);
  return b*Vector3d(4.*(1. + x[1]) + double(k)*x[0], 1. + x[0] - double(k)*x[1]*x[1],
                    D == 3 ? 1. + x[0]*x[2] + 0.5*double(k) : 0.);
}

// A wall-upgraded P2 case of two stamps, read by the stamped loader
// (h5_params.dat) and as components by the frequency loader (freq_params.dat)
template<typename Cell>
void write_upgraded(const CaseDir& c, const std::size_t n, const Field& f, const bool upgrade = true){
  const std::string key = upgrade ? "wall_upgrade=true\n" : "wall_upgrade=false\n";
  write_stamped<Cell>(c, n, "P2", walls, true, {}, f);
  std::ofstream(c.file("h5_params.dat"), std::ios::app) << key;
  write_freq_of_stamps(c, "P2", walls, "0 0 0.7 up_0.h5\n2 0.3 -1.3 up_1.h5\n", key);
}

template<typename I>
PointValues held_at(I& intp, const Region& R, const Vector3d& x, const double t){
  PointValues v(1.);
  HeldEval<I> ev(intp, v);
  ev.bind(R);
  ev.at(x, t);
  return v;
}

// The point of cell R at barycentrics lam, from a point x0 in it
template<typename I>
Vector3d point_at(I& intp, const Region& R, const Vector3d& x0, const std::array<double, 4>& lam, const int D){
  std::array<double, 4> l0;
  intp.region_point(R, x0, l0);
  Matrix3d M = Matrix3d::Identity();
  Vector3d r = Vector3d::Zero();
  for (int i = 0; i < D; ++i){
    const Vector3d g = intp.level_grad(R, i + 1);
    for (int d = 0; d < D; ++d) M(i, d) = g[d];
    r[i] = lam[std::size_t(i + 1)] - l0[std::size_t(i + 1)];
  }
  return x0 + M.partialPivLu().solve(r);
}

// Points located in cells with a vertex on the box's sides, and their regions
template<typename I>
std::vector<std::pair<Region, Vector3d>> wall_cells(I& intp, const int D, const std::size_t count){
  std::vector<std::pair<Region, Vector3d>> out;
  for (const Vector3d& x : box_points(intp, D, walls, 20*count)){
    CellPos pos;
    if (!intp.locate(x, T0, pos)) continue;
    const Region R = intp.region_in(pos.id, x);
    bool wall = false;
    for (int k = 0; k <= D; ++k){
      std::array<double, 4> e{};
      e[std::size_t(k)] = 1.;
      const Vector3d v = point_at(intp, R, x, e, D);
      for (int d = 0; d < D; ++d) wall = wall || std::abs(v[d]) < 1e-9 || std::abs(v[d] - 1.) < 1e-9;
    }
    if (wall) out.push_back({R, x});
    if (out.size() == count) break;
  }
  return out;
}

// Across each inner face of a cell with a wall vertex, the field from either
// side, and across a plane of its split (two least barycentrics equal), the
// field a step 1e-10 to either side: the same to round-off, the gradient not
// across the plane
template<typename Cell>
void continuous(const std::size_t n){
  constexpr int D = Cell::n_verts - 1;
  CaseDir c("wsplit_cont");
  write_upgraded<Cell>(c, n, at_rest<D>);
  SimplexInterpol<Cell> intp(c.file("h5_params.dat"));
  intp.set_int_order(2);
  intp.set_needs_gradient(true);
  intp.update(T0);
  std::mt19937 rng(5);
  double scale = 0., gscale = 0., face = 0., plane = 0., kink = 0.;
  std::size_t faces = 0, planes = 0;
  for (const auto& [R, x] : wall_cells(intp, D, 200)){
    for (int k = 0; k <= D; ++k){
      if (intp.is_wall(R, k)) continue;
      std::array<double, 4> lam = random_weights<Cell>(rng);
      const double s = 1. - lam[std::size_t(k)];
      for (int i = 0; i <= D; ++i) lam[std::size_t(i)] = i == k ? 0. : lam[std::size_t(i)]/s;
      const Vector3d y = point_at(intp, R, x, lam, D);
      Region S;
      REQUIRE(intp.across(R, k, y, S) >= 0);
      const PointValues a = held_at(intp, R, y, T0), b = held_at(intp, S, y, T0);
      scale = std::max(scale, largest(a.U));
      gscale = std::max(gscale, largest(a.gradU));
      face = std::max(face, largest(a.U - b.U));
      ++faces;
    }
    // The plane between the sub-cells of the two least barycentrics
    std::array<double, 4> lam = random_weights<Cell>(rng);
    int j = 0;
    for (int i = 1; i <= D; ++i) if (lam[std::size_t(i)] < lam[std::size_t(j)]) j = i;
    const int a = (j + 1) % (D + 1);
    lam[std::size_t(a)] = lam[std::size_t(j)];
    double s = 0.;
    for (int i = 0; i <= D; ++i) s += lam[std::size_t(i)];
    for (int i = 0; i <= D; ++i) lam[std::size_t(i)] /= s;
    const Vector3d y = point_at(intp, R, x, lam, D);
    const Vector3d e = 1e-10*(intp.level_grad(R, a) - intp.level_grad(R, j)).normalized();
    const PointValues p = held_at(intp, R, Vector3d(y + e), T0), q = held_at(intp, R, Vector3d(y - e), T0);
    plane = std::max(plane, largest(p.U - q.U));
    kink = std::max(kink, largest(p.gradU - q.gradU));
    ++planes;
  }
  INFO("faces " << faces << ", planes " << planes << ", scale " << scale << ", gradient scale "
       << gscale);
  REQUIRE(faces > 100);
  CHECK(face <= 1e-13*scale);
  CHECK(plane <= 1e-8*scale);
  // The correction is there: P2's gradient is smooth in the cell, the field's jumps
  CHECK(kink > 1e-4*gscale);
}

// div u on the walls, at random points of the box's sides: zero to round-off
// with the rule, P2's own O(1) without it; the stamped and the frequency
// loader
template<typename Cell>
void divergence_free_on_walls(const std::size_t n){
  constexpr int D = Cell::n_verts - 1;
  CaseDir c("wsplit_div"), p("wsplit_div_plain");
  write_upgraded<Cell>(c, n, at_rest<D>);
  write_upgraded<Cell>(p, n, at_rest<D>, false);
  SimplexInterpol<Cell> st(c.file("h5_params.dat")), st_plain(p.file("h5_params.dat"));
  SimplexFreqInterpol<Cell> fq(c.file("freq_params.dat")), fq_plain(p.file("freq_params.dat"));
  std::mt19937 rng(9);
  std::uniform_real_distribution<double> uni(0.01, 0.99);
  std::vector<Vector3d> pts;
  for (int i = 0; i < 400; ++i){
    Vector3d x = Vector3d::Zero();
    for (int d = 0; d < D; ++d) x[d] = uni(rng);
    x[i % D] = (i/D) % 2 ? 1. : 0.;
    pts.push_back(x);
  }
  const auto check = [&](auto& up, auto& plain){
    double div = 0., div_plain = 0., g = 0.;
    for (auto* i : {&up, &plain}){
      i->set_int_order(2);
      i->set_needs_gradient(true);
      i->update(T0);
    }
    for (const Vector3d& x : pts){
      CellPos pos;
      REQUIRE(up.locate(x, T0, pos));
      PointValues a(1.), b(1.);
      up.evaluate_motion(x, T0, pos, a);
      plain.evaluate_motion(x, T0, pos, b);
      g = std::max(g, largest(a.gradU));
      div = std::max({div, std::abs(a.gradU.trace()), std::abs(a.gradA.trace())});
      div_plain = std::max(div_plain, std::abs(b.gradU.trace()));
    }
    INFO("gradient scale " << g << ", div " << div << ", P2's " << div_plain);
    CHECK(div <= 1e-12*g);
    CHECK(div_plain > 1e-2*g);
  };
  check(st, st_plain);
  check(fq, fq_plain);
}

// A field at rest at every wall vertex, with every wall edge off y = 0 and
// y = 1 moving: the facets on y = 0 and 1 alone are at rest, so div u
// vanishes on them and not on the other sides; the cells whose wall facets
// lie on the other sides keep the wall terms of their wall vertices at rest
// (in 2D also a correction), so they differ from the plain P2 field
template<typename Cell>
void midpoint_off_rest(const std::size_t n){
  constexpr int D = Cell::n_verts - 1;
  const double m = M_PI*double(n);
  const Field f = [m](const Vector3d& x, const int k){
    double s = std::pow(std::sin(m*x[1]), 2);
    for (int d = 0; d < D; ++d)
      if (d != 1) s += x[1]*(1. - x[1])*std::pow(std::sin(m*x[d]), 2);
    return Vector3d(at_rest<D>(x, k) + 0.3*s*Vector3d(1., 0.5, D == 3 ? 0.2 : 0.));
  };
  CaseDir c("wsplit_mid"), p("wsplit_mid_plain");
  write_upgraded<Cell>(c, n, f);
  write_upgraded<Cell>(p, n, f, false);
  SimplexInterpol<Cell> st(c.file("h5_params.dat")), st_plain(p.file("h5_params.dat"));
  SimplexFreqInterpol<Cell> fq(c.file("freq_params.dat")), fq_plain(p.file("freq_params.dat"));
  const double reach = 1./double(n);
  std::mt19937 rng(13);
  std::uniform_real_distribution<double> uni(0.02, 0.98);
  const auto check = [&](auto& up, auto& plain){
    for (auto* i : {&up, &plain}){
      i->set_int_order(2);
      i->set_needs_gradient(true);
      i->update(T0);
    }
    std::size_t moving = 0, corrected = 0;
    for (const Vector3d& x : box_points(up, D, walls, 4000)){
      CellPos pos;
      REQUIRE(up.locate(x, T0, pos));
      PointValues a(1.), b(1.);
      up.evaluate_motion(x, T0, pos, a);
      plain.evaluate_motion(x, T0, pos, b);
      const bool same = a.U == b.U && a.gradU == b.gradU && a.A == b.A && a.gradA == b.gradA;
      double others = 1.;
      for (int d = 0; d < D; ++d)
        if (d != 1) others = std::min({others, x[d], 1. - x[d]});
      const double ys = std::min(x[1], 1. - x[1]);
      if (others < reach && ys > reach && !same) ++moving;
      else if (ys < reach && others > reach && !same) ++corrected;
    }
    CHECK(moving > 50);
    CHECK(corrected > 50);
    // div u on the sides, away from their edges
    double rest = 0., off = 0., g = 0.;
    for (int i = 0; i < 400; ++i){
      Vector3d x = Vector3d::Zero();
      for (int d = 0; d < D; ++d) x[d] = 0.1 + 0.8*uni(rng);
      const int side = i % D;
      x[side] = (i/D) % 2 ? 1. : 0.;
      CellPos pos;
      REQUIRE(up.locate(x, T0, pos));
      PointValues a(1.);
      up.evaluate_motion(x, T0, pos, a);
      g = std::max(g, a.gradU.cwiseAbs().maxCoeff());
      if (side == 1) rest = std::max(rest, std::abs(a.gradU.trace()));
      else off = std::max(off, std::abs(a.gradU.trace()));
    }
    INFO("gradient scale " << g << ", div at rest " << rest << ", moving " << off);
    CHECK(rest <= 1e-12*g);
    CHECK(off > 1e-3*g);
  };
  check(st, st_plain);
  check(fq, fq_plain);
}

// The first stamp at rest on the facets of z = 0 and 1 (y = 0 and 1 in 2D)
// and on some wall edges (every wall vertex in 2D) of the other sides, whose
// facets move; the second moving those edges (vertices) and keeping the
// facets: the stamped loader refuses it, naming the edges (vertices) and the
// file. With vertices, in 3D, the second moves the other sides' wall
// vertices too, which the cells touching a side at a vertex alone name.
template<typename Cell>
void entity_moved(const std::size_t n, const bool vertices){
  constexpr int D = Cell::n_verts - 1;
  const double m = M_PI*double(n);
  const Field f = [m, vertices](const Vector3d& x, const int k){
    const double t = x[D - 1]*(1. - x[D - 1]);
    double s;
    if constexpr (D == 3)
      s = std::pow(std::sin(m*x[2]), 2)
          + double(k)*t*std::pow(vertices ? std::cos(m*x[1]) : std::sin(m*x[1]), 2);
    else
      s = std::pow(std::sin(m*x[1]), 2) + double(k)*t*std::pow(std::cos(m*x[0])*std::cos(m*x[1]), 2);
    return Vector3d(at_rest<D>(x, k) + 0.3*s*Vector3d(1., 0.5, D == 3 ? 0.2 : 0.));
  };
  CaseDir c("wsplit_entity");
  write_upgraded<Cell>(c, n, f);
  SimplexInterpol<Cell> st(c.file("h5_params.dat"));
  std::string refused;
  try { st.update(T0); }
  catch (const partrac::Error& e){ refused = e.what(); }
  INFO(refused);
  CHECK(refused.find("up_1.h5") != std::string::npos);
  CHECK(refused.find(D == 3 ? "wall edges at rest" : "wall vertices at rest") != std::string::npos);
  CHECK(refused.find("wall facets") == std::string::npos);
  if (D == 3) CHECK((refused.find("wall vertices at rest") != std::string::npos) == vertices);
}

}  // namespace

TEST_CASE("The wall correction is continuous across the faces of the cells with a wall vertex and across "
          "its split's planes", "[wall_split]"){
  SECTION("triangles"){ continuous<Triangle>(6); }
  SECTION("tets"){ continuous<Tet>(3); }
}

TEST_CASE("With the wall correction div u vanishes on the walls", "[wall_split]"){
  SECTION("triangles"){ divergence_free_on_walls<Triangle>(6); }
  SECTION("tets"){ divergence_free_on_walls<Tet>(3); }
}

TEST_CASE("A wall entity with a node off rest is not at rest", "[wall_split]"){
  SECTION("triangles"){ midpoint_off_rest<Triangle>(6); }
  SECTION("tets"){ midpoint_off_rest<Tet>(3); }
}

TEST_CASE("A later stamp that moves a wall edge or vertex at rest in the first is refused", "[wall_split]"){
  SECTION("triangles"){ entity_moved<Triangle>(6, false); }
  SECTION("tets"){ entity_moved<Tet>(3, false); }
  SECTION("tets, wall vertices"){ entity_moved<Tet>(3, true); }
}

#endif
