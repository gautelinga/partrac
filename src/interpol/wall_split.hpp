#ifndef __WALL_SPLIT_HPP
#define __WALL_SPLIT_HPP

// The wall correction of a P2 velocity in a cell with a wall vertex at rest:
// u = P2 + w + g, w continuous quartic wall terms vanishing at the P2 nodes, g
// zero on the cell's faces, the least |grad g|^2 with div u = 0 at the walls at
// rest and n.grad div u = 0 on the wall facets; on sub-cell j, g = lambda_j E Q_j

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

#include "split_eval.hpp"
#include "typedefs.hpp"

namespace wall_split {

// The degree of the correction on each sub-cell
constexpr int K = 4;
// A factored table's rank is padded to a block
constexpr int block = 4;

// Two doubles: the two fields' lanes, and the tables' kernels' unit; a vector
// where the target has them
typedef double pair __attribute__((vector_size(16)));
inline pair load(const double* p){ pair a; std::memcpy(&a, p, sizeof a); return a; }
inline void store(double* p, const pair a){ std::memcpy(p, &a, sizeof a); }
inline pair splat(const double x){ return pair{x, x}; }
// The two lanes' sums of a and b
inline pair sums(const pair a, const pair b){ return pair{a[0], b[0]} + pair{a[1], b[1]}; }

constexpr int binom(const int n, const int k){
  int r = 1;
  for (int i = 1; i <= k; ++i) r = r*(n - k + i)/i;
  return r;
}

// A sub-cell's degree-K lattice over its barycentrics (the macro vertices
// but j in order, then the centroid), lexicographic, and the nodes off the
// macro faces: the centroid's index at least 1
template<int NV>
struct Lattice {
  static constexpr int n_loc = binom(K + NV - 1, NV - 1);
  static constexpr int n_act = n_loc - binom(K + NV - 2, NV - 2);
  std::array<std::array<int, NV>, n_loc> all{};
  std::array<std::array<int, NV>, n_act> act{};
  constexpr Lattice(){
    int c = 0, a = 0;
    int pw = 1;
    for (int i = 0; i < NV; ++i) pw *= K + 1;
    for (int code = 0; code < pw; ++code){
      std::array<int, NV> m{};
      int r = code, s = 0;
      for (int i = NV - 1; i >= 0; --i){ m[std::size_t(i)] = r % (K + 1); r /= K + 1; s += m[std::size_t(i)]; }
      if (s != K) continue;
      all[std::size_t(c++)] = m;
      if (m[NV - 1] >= 1) act[std::size_t(a++)] = m;
    }
  }
};

// The split's interior nodes: sub-cell j's active node l is interior node
// idx[j][l], numbered by sub-cell then node, first seen first
template<int NV>
struct Interior {
  static constexpr Lattice<NV> lat{};
  static constexpr int n_act = Lattice<NV>::n_act;
  static constexpr int n_max = NV*n_act;
  std::array<std::array<int, n_act>, NV> idx{};
  int n = 0;
  constexpr Interior(){
    std::array<std::array<int, NV>, n_max> keys{};
    for (int j = 0; j < NV; ++j)
      for (int l = 0; l < n_act; ++l){
        const auto& a = lat.act[std::size_t(l)];
        // NV K lambda, an integer per macro vertex
        std::array<int, NV> key{};
        int s = 0;
        for (int m = 0; m < NV; ++m){
          if (m == j){ key[std::size_t(m)] = a[NV - 1]; continue; }
          key[std::size_t(m)] = NV*a[std::size_t(s++)] + a[NV - 1];
        }
        int found = -1;
        for (int q = 0; q < n && found < 0; ++q){
          bool same = true;
          for (int m = 0; m < NV; ++m) same = same && keys[std::size_t(q)][std::size_t(m)] == key[std::size_t(m)];
          if (same) found = q;
        }
        if (found < 0){ keys[std::size_t(n)] = key; found = n++; }
        idx[std::size_t(j)][std::size_t(l)] = found;
      }
  }
};

// The vertices off each lexicographic edge, ascending
template<int NV>
constexpr std::array<std::array<int, NV - 2>, NV*(NV - 1)/2> edge_apices(){
  constexpr auto edges = split_eval::lex_edges<NV>();
  std::array<std::array<int, NV - 2>, NV*(NV - 1)/2> t{};
  for (std::size_t e = 0; e < edges.size(); ++e){
    int j = 0;
    for (int v = 0; v < NV; ++v)
      if (v != edges[e][0] && v != edges[e][1]) t[e][std::size_t(j++)] = v;
  }
  return t;
}

template<int NV>
struct Shape {
  static constexpr int nv = NV, D = NV - 1;
  static constexpr Lattice<NV> lat{};
  static constexpr Interior<NV> in{};
  static constexpr int n_act = Lattice<NV>::n_act;
  static constexpr int n_int = in.n;
  // Pattern bits: the facets (bit k faces vertex k), then the edges 01, 02,
  // 03, 12, 13, 23 in 3D, the vertices in 2D; the rest bits add the
  // vertices in 3D
  static constexpr int n_ent = NV == 4 ? 6 : 3;
  static constexpr int n_bits = NV + n_ent;
  static constexpr int n_rest = n_bits + (NV == 4 ? NV : 0);
  static constexpr unsigned pattern_mask = (1u << n_bits) - 1;
  // The wall terms by bit: edge ends 2e + o (o: the end at the edge's second
  // vertex), side faces e (NV - 2) + apex, then the bubble's facets (3D)
  static constexpr int n_edges = NV*(NV - 1)/2;
  static constexpr int n_ends = 2*n_edges;
  static constexpr int n_triples = n_edges*(NV - 2);
  static constexpr int n_terms = n_ends + n_triples + (NV == 4 ? NV : 0);
  static constexpr int triple0 = n_ends, bubble0 = n_ends + n_triples;
  // Their monomials, the slots: the edge ends, the face products times a
  // vertex of the face (3D: face k's three; 2D: the cell's), the bubble (3D)
  static constexpr int n_faces = NV == 4 ? 12 : 3;
  static constexpr int n_slots = n_ends + n_faces + (NV == 4 ? 1 : 0);
  static constexpr int face0 = n_ends, bubble = n_ends + n_faces;
  // A table: per interior node Wd (D x NV), then Wc (D x D) a slot
  static constexpr std::size_t n_d = std::size_t(D*NV);
  static constexpr std::size_t node = n_d + std::size_t(n_slots*D*D);
  static constexpr std::size_t stride = std::size_t(n_int)*node;
  // A factored table's rank, padded to a block, at most
  static constexpr int r_max = (n_int*D + block - 1)/block*block;
};

// The homogeneous monomials of degree d in NV barycentrics, lexicographic
// from the highest power of the first
template<int NV, int d>
struct Degree {
  static constexpr int n = binom(d + NV - 1, NV - 1);
  std::array<std::array<int, NV>, n> e{};
  constexpr Degree(){
    int c = 0, pw = 1;
    for (int i = 0; i < NV; ++i) pw *= d + 1;
    for (int code = pw - 1; code >= 0; --code){
      std::array<int, NV> a{};
      int r = code, s = 0;
      for (int i = NV - 1; i >= 0; --i){ a[std::size_t(i)] = r % (d + 1); r /= d + 1; s += a[std::size_t(i)]; }
      if (s == d) e[std::size_t(c++)] = a;
    }
  }
  constexpr int index(const std::array<int, NV>& a) const {
    for (int k = 0; k < n; ++k){
      bool same = true;
      for (int i = 0; i < NV; ++i) same = same && e[std::size_t(k)][std::size_t(i)] == a[std::size_t(i)];
      if (same) return k;
    }
    return -1;
  }
};

// The homogeneous forms of a wall cell's field in its barycentrics: P2 of
// degree 2, the rest of degree K = 4. Each monomial of degree 2 or 3 as the
// product of one of the degree below and a variable (lo, var); per variable,
// a cubic's (linear's) monomial times it as a quartic's (quadratic's), up3
// (up1)
template<int NV>
struct Forms {
  static constexpr Degree<NV, 1> d1{};
  static constexpr Degree<NV, 2> d2{};
  static constexpr Degree<NV, 3> d3{};
  static constexpr Degree<NV, 4> d4{};
  static constexpr int n2 = d2.n, n3 = d3.n, n4 = d4.n;
  std::array<int, n2> lo2{}, var2{};
  std::array<int, n3> lo3{}, var3{};
  std::array<std::array<int, n3>, NV> up3{};
  std::array<std::array<int, NV>, NV> up1{};
  template<typename Lo, typename Hi, std::size_t N>
  static constexpr void split(const Lo& lo, const Hi& hi, std::array<int, N>& l, std::array<int, N>& v){
    for (std::size_t k = 0; k < N; ++k){
      std::array<int, NV> a = hi.e[k];
      int i = NV - 1;
      while (a[std::size_t(i)] == 0) --i;
      --a[std::size_t(i)];
      l[k] = lo.index(a);
      v[k] = i;
    }
  }
  constexpr Forms(){
    split(d1, d2, lo2, var2);
    split(d2, d3, lo3, var3);
    for (int i = 0; i < NV; ++i){
      for (int k = 0; k < n3; ++k){
        std::array<int, NV> a = d3.e[std::size_t(k)];
        ++a[std::size_t(i)];
        up3[std::size_t(i)][std::size_t(k)] = d4.index(a);
      }
      for (int k = 0; k < NV; ++k){
        std::array<int, NV> a{};
        ++a[std::size_t(k)];
        ++a[std::size_t(i)];
        up1[std::size_t(i)][std::size_t(k)] = d2.index(a);
      }
    }
  }
};

// The face slot of vertex x on the face lacking k (2D: the cell's)
template<int NV>
constexpr int face_slot(const int k, const int x){
  if constexpr (NV == 3) return Shape<NV>::face0 + x;
  else return Shape<NV>::face0 + 3*k + (x < k ? x : x - 1);
}

// The macro vertices of each codimension-2 entity: an edge's two, a vertex twice
template<int NV>
constexpr std::array<std::array<int, 2>, Shape<NV>::n_ent> entity_verts(){
  std::array<std::array<int, 2>, Shape<NV>::n_ent> t{};
  if constexpr (NV == 4){
    int e = 0;
    for (int a = 0; a < 4; ++a)
      for (int b = a + 1; b < 4; ++b) t[std::size_t(e++)] = {a, b};
  }
  else
    for (int v = 0; v < 3; ++v) t[std::size_t(v)] = {v, v};
  return t;
}

// The P2 slots of each rest bit: a facet's vertices and edge midpoints, an
// edge's ends and midpoint, a vertex
template<typename Cell>
constexpr std::array<unsigned, Shape<Cell::n_verts>::n_rest> entity_slots(){
  constexpr int nv = Cell::n_verts;
  constexpr auto ev = entity_verts<nv>();
  constexpr auto edges = split_eval::lex_edges<nv>();
  std::array<unsigned, Shape<nv>::n_rest> t{};
  for (int k = 0; k < nv; ++k){
    for (int j = 0; j < nv; ++j)
      if (j != k) t[std::size_t(k)] |= 1u << j;
    for (std::size_t e = 0; e < edges.size(); ++e)
      if (edges[e][0] != k && edges[e][1] != k) t[std::size_t(k)] |= 1u << Cell::mid_[e];
  }
  for (std::size_t q = 0; q < ev.size(); ++q){
    const int a = ev[q][0], b = ev[q][1];
    unsigned m = (1u << a) | (1u << b);
    if (a != b) m |= 1u << Cell::mid_[std::size_t(split_eval::lex_index<nv>(a, b))];
    t[std::size_t(nv) + q] = m;
  }
  for (int v = Shape<nv>::n_bits; v < Shape<nv>::n_rest; ++v)
    t[std::size_t(v)] = 1u << (v - Shape<nv>::n_bits);
  return t;
}

// The entities of candidates (rest bits) at rest in a P2 block: all their
// P2 nodes within tol, D values a node
template<typename Cell>
inline unsigned rest_entities(const double* u, const unsigned candidates, const double tol)
{
  constexpr int D = Cell::n_verts - 1;
  constexpr std::size_t n2 = Cell::n_dofs_max;
  static constexpr auto slots = entity_slots<Cell>();
  unsigned nodes = 0;
  for (std::size_t q = 0; q < n2; ++q){
    bool r = true;
    for (int d = 0; d < D; ++d) r = r && std::abs(u[std::size_t(d)*n2 + q]) <= tol;
    if (r) nodes |= 1u << q;
  }
  unsigned rest = 0;
  for (unsigned m = candidates; m; m &= m - 1){
    const int k = __builtin_ctz(m);
    if ((nodes & slots[std::size_t(k)]) == slots[std::size_t(k)]) rest |= 1u << k;
  }
  return rest;
}

// A wall cell: its wall terms' normals in the pool (D doubles each; a bubble
// facet's normal, then its column of the facets' pseudo-inverse), its wall
// terms by bit, its pattern's table (none: the wall terms alone) and its wall
// entities at rest (the candidates until they are fixed)
struct WallCell {
  std::uint32_t normals = 0;
  std::uint32_t terms = 0;
  std::uint16_t table = 0;
  std::uint16_t rest = 0;
};
constexpr std::uint16_t no_table = 0xffff;

// A pattern's table factored, of rank r padded to a block, per sub-cell j in
// a basis of its own: z = Vd d + sum_s Vs_s h_s (Vd NV x r, Vs_s D x r, rows
// of r) for each j, then for each j the coefficients of Q_j (a row a cubic
// monomial and component) by rows of r: of those z gives alone, idx the rows
// (the pivots), those of the others, in fours, their rows idx too
struct Table {
  const double* p = nullptr;
  int r = 0;
  const std::int16_t* idx = nullptr;
};
template<int NV>
constexpr std::size_t table_inputs(const int r){
  return std::size_t(r)*std::size_t(NV + Shape<NV>::n_slots*(NV - 1));
}
template<int NV>
constexpr std::size_t table_size(const int r){
  return std::size_t(NV)*(table_inputs<NV>(r) + std::size_t(r)*std::size_t(Forms<NV>::n3*(NV - 1)));
}
// A sub-cell's indices: the rank, the other rows' count padded to four, the
// pivots' rows, the others' rows (the last repeated to pad)
template<int NV>
constexpr int idx_stride = Forms<NV>::n3*(NV - 1) + 6;

// The wall cells of a mesh, their normals and the tables of their patterns
template<int NV>
struct Tables {
  std::vector<WallCell> cell;           // by wall cell
  std::vector<double> normals;
  std::vector<std::uint16_t> pattern;   // by table
  std::vector<double> W;
  std::vector<std::size_t> at;          // by table: its place in W
  std::vector<int> rank;                // and its rank, padded
  std::vector<std::int16_t> idx;        // NV idx_stride a table
  Table of(const WallCell& c) const {
    return {W.data() + at[c.table], rank[c.table], idx.data() + std::size_t(c.table)*NV*idx_stride<NV>};
  }
  const double* normals_of(const WallCell& c) const { return normals.data() + std::size_t(c.normals)*(NV - 1); }
  // The tables of the nonzero patterns at rest, ascending, each cell's set
  void build(const bool verbose);
};

// The problem on the regular reference simplex: the stiffness of the split's
// nodes (the interior ones first, a scalar per pair), the condition rows over
// every node (D columns a node), their targets per vertex divergence and
// pattern bits, the nodes' barycentrics and the edge matrix
struct Reference {
  Eigen::MatrixXd S, B, T, Eh;
  std::vector<int> bit;
  std::vector<std::array<double, 4>> lam;
  int n_int = 0;
};
template<int NV>
Reference reference();

// The tables of the patterns given, in their order: per interior node in
// edge coordinates the vertex divergences' (D x NV) and the slots' (D x D a
// slot, in edge coordinates too)
template<int NV>
std::vector<double> pattern_tables(const std::vector<std::uint16_t>& patterns);

// A table of pattern_tables factored, its rank padded to r, and its indices
template<int NV>
std::vector<double> factor(const double* w, int& r, std::vector<std::int16_t>& idx);

// The cell's barycentric gradients and its edge matrix E (column i: vertex
// i+1 minus vertex 0), the inverse of the gradients 1..D as rows
template<typename Cell>
struct Frame {
  static constexpr int nv = Cell::n_verts, D = nv - 1;
  Vector3d G[nv];
  double E[D][D];
  Frame() = default;
  explicit Frame(const Cell& cell){
    for (int k = 0; k < nv; ++k) G[k] = cell.bary_grad(k);
    if constexpr (D == 2){
      const double a = G[1][0], b = G[1][1], c = G[2][0], d = G[2][1];
      const double r = 1./(a*d - b*c);
      E[0][0] = d*r;  E[0][1] = -b*r;
      E[1][0] = -c*r; E[1][1] = a*r;
    }
    else {
      const Vector3d& p = G[1];
      const Vector3d& q = G[2];
      const Vector3d& s = G[3];
      const Vector3d c0 = q.cross(s), c1 = s.cross(p), c2 = p.cross(q);
      const double r = 1./p.dot(c0);
      for (int x = 0; x < 3; ++x){
        E[x][0] = c0[x]*r;
        E[x][1] = c1[x]*r;
        E[x][2] = c2[x]*r;
      }
    }
  }
};

// div of a P2 block at the vertices, D values a node:
// d_i = 4 (G_i.u_i + sum_{edges ib} G_b.u_ib) - sum_n G_n.u_n
template<typename Cell>
inline void vertex_divergence(const Frame<Cell>& f, const double* u, double* d)
{
  constexpr int nv = Cell::n_verts, D = nv - 1;
  constexpr std::size_t n2 = Cell::n_dofs_max;
  static constexpr auto edges = split_eval::lex_edges<nv>();
  double s = 0.;
  for (int i = 0; i < nv; ++i){
    double a = 0.;
    for (int x = 0; x < D; ++x) a += f.G[i][x]*u[std::size_t(x)*n2 + std::size_t(i)];
    d[i] = 4.*a;
    s += a;
  }
  for (std::size_t e = 0; e < edges.size(); ++e){
    const int a = edges[e][0], b = edges[e][1];
    const std::size_t m = std::size_t(Cell::mid_[e]);
    double ga = 0., gb = 0.;
    for (int x = 0; x < D; ++x){
      const double v = u[std::size_t(x)*n2 + m];
      ga += f.G[a][x]*v;
      gb += f.G[b][x]*v;
    }
    d[a] += 4.*gb;
    d[b] += 4.*ga;
  }
  for (int i = 0; i < nv; ++i) d[i] -= s;
}

// Slot s's monomial at the barycentrics l, and its derivatives in them
template<int NV>
inline double slot_weight(const int s, const double* l, double* dl)
{
  using S = Shape<NV>;
  static constexpr auto edges = split_eval::lex_edges<NV>();
  for (int m = 0; m < NV; ++m) dl[m] = 0.;
  if (s < S::n_ends){
    const int w = edges[std::size_t(s/2)][std::size_t(s%2)], v = edges[std::size_t(s/2)][std::size_t(1 - s%2)];
    const double lw = l[w], lv = l[v];
    dl[w] = 3.*lw*lw*lv - 2.*lw*lv*lv;
    dl[v] = lw*lw*lw - 2.*lw*lw*lv;
    return lw*lw*lv*(lw - lv);
  }
  if (s < S::bubble){
    // the face lacking k, times its vertex x
    int k, x;
    if constexpr (NV == 3){ k = 3; x = s - S::face0; }
    else {
      k = (s - S::face0)/3;
      const int i = (s - S::face0) % 3;
      x = i < k ? i : i + 1;
    }
    int f[3], q = 0;
    for (int m = 0; m < NV; ++m) if (m != k) f[q++] = m;
    const double a = l[f[0]], b = l[f[1]], c = l[f[2]], lx = l[x];
    dl[f[0]] = b*c*lx;
    dl[f[1]] = a*c*lx;
    dl[f[2]] = a*b*lx;
    dl[x] += a*b*c;
    return a*b*c*lx;
  }
  double p = 1.;
  for (int m = 0; m < NV; ++m){
    double o = 1.;
    for (int i = 0; i < NV; ++i) if (i != m) o *= l[i];
    dl[m] = o;
    p *= l[m];
  }
  return p;
}

// Every slot's monomial at the barycentrics l
template<int NV>
inline void slot_weights(const double* l, double* w)
{
  using S = Shape<NV>;
  static constexpr auto edges = split_eval::lex_edges<NV>();
  for (int e = 0; e < S::n_edges; ++e){
    const double a = l[edges[std::size_t(e)][0]], b = l[edges[std::size_t(e)][1]], ab = a*b;
    w[2*e] = a*ab*(a - b);
    w[2*e + 1] = b*ab*(b - a);
  }
  if constexpr (NV == 3){
    const double p = l[0]*l[1]*l[2];
    for (int x = 0; x < 3; ++x) w[S::face0 + x] = p*l[x];
  }
  else {
    const double l01 = l[0]*l[1], l23 = l[2]*l[3];
    const double f[4] = {l[1]*l23, l[0]*l23, l01*l[3], l01*l[2]};
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < 3; ++i) w[S::face0 + 3*k + i] = f[k]*l[i < k ? i : i + 1];
    w[S::bubble] = l01*l23;
  }
}

// The wall terms' coefficients of two fields in a wall cell, D a slot, a
// field a lane
template<int NV>
struct Terms {
  static constexpr int D = NV - 1;
  static constexpr int n_pairs = Shape<NV>::n_ends + 2*Shape<NV>::n_triples + (NV == 4 ? NV : 0);
  std::uint32_t slots = 0;
  pair c[Shape<NV>::n_slots][D];
  // Each term's part of a slot: the slot, its direction and the two fields' factors
  int n = 0;
  int s[n_pairs];
  const double* dir[n_pairs];
  pair q[n_pairs];
};

// The P2 slot of the edge between two vertices
template<typename Cell>
constexpr std::array<std::array<int, 4>, 4> mid_slots(){
  constexpr int nv = Cell::n_verts;
  std::array<std::array<int, 4>, 4> t{};
  for (int a = 0; a < nv; ++a)
    for (int c = a + 1; c < nv; ++c)
      t[std::size_t(a)][std::size_t(c)] = t[std::size_t(c)][std::size_t(a)]
        = Cell::mid_[std::size_t(split_eval::lex_index<nv>(a, c))];
  return t;
}

// From the P2 blocks u and v, the cell's terms and their normals n
template<typename Cell>
inline void term_coefficients(const double* n, const std::uint32_t terms, const double* u, const double* v,
                             Terms<Cell::n_verts>& o)
{
  constexpr int nv = Cell::n_verts, D = nv - 1;
  constexpr std::size_t n2 = Cell::n_dofs_max;
  using S = Shape<nv>;
  static constexpr auto edges = split_eval::lex_edges<nv>();
  static constexpr auto apex = edge_apices<nv>();
  static constexpr auto mid = mid_slots<Cell>();
  // P2's slope from a to c, 4 u_m - u_c, along nn
  const auto slope = [&](const double* nn, const int a, const int c){
    const std::size_t m = std::size_t(mid[std::size_t(a)][std::size_t(c)]);
    pair p = splat(0.);
    for (int d = 0; d < D; ++d){
      const std::size_t i = std::size_t(d)*n2;
      p += nn[d]*(4.*pair{u[i + m], v[i + m]} - pair{u[i + std::size_t(c)], v[i + std::size_t(c)]});
    }
    return p;
  };
  std::uint32_t slots = terms & ((1u << S::n_ends) - 1);
  int np = 0;
  for (std::uint32_t m = slots; m; m &= m - 1, n += D, ++np){
    const int t = __builtin_ctz(m);
    const int w = edges[std::size_t(t/2)][std::size_t(t%2)], c = edges[std::size_t(t/2)][std::size_t(1 - t%2)];
    const pair p = slope(n, w, c);
    for (int d = 0; d < D; ++d) o.c[t][d] = -p*n[d];
    o.q[np] = -p;
    o.s[np] = t;
    o.dir[np] = n;
  }
  // Side faces: two slots each, which another side face of the face may share
  const std::uint32_t triples = terms >> S::triple0 & ((1u << S::n_triples) - 1);
  for (std::uint32_t m = triples; m; m &= m - 1){
    const int t = __builtin_ctz(m), e = t/(nv - 2);
    const int a = edges[std::size_t(e)][0], c = edges[std::size_t(e)][1];
    const int k = nv == 3 ? 3 : 6 - a - c - apex[std::size_t(e)][std::size_t(t % (nv - 2))];
    for (const int s : {face_slot<nv>(k, a), face_slot<nv>(k, c)}){
      if (slots >> s & 1) continue;
      slots |= 1u << s;
      for (int d = 0; d < D; ++d) o.c[s][d] = splat(0.);
    }
  }
  for (std::uint32_t m = triples; m; m &= m - 1, n += D){
    const int t = __builtin_ctz(m), e = t/(nv - 2);
    const int a = edges[std::size_t(e)][0], c = edges[std::size_t(e)][1];
    const int x = apex[std::size_t(e)][std::size_t(t % (nv - 2))];
    const int k = nv == 3 ? 3 : 6 - a - c - x;
    const int sa = face_slot<nv>(k, a), sc = face_slot<nv>(k, c);
    const pair pa = slope(n, a, x), pc = slope(n, c, x);
    const pair qa = -(2.*pa + pc), qc = -(pa + 2.*pc);
    for (int d = 0; d < D; ++d){
      o.c[sa][d] += qa*n[d];
      o.c[sc][d] += qc*n[d];
    }
    o.s[np] = sa;
    o.s[np + 1] = sc;
    o.dir[np] = o.dir[np + 1] = n;
    o.q[np] = qa;
    o.q[np + 1] = qc;
    np += 2;
  }
  // The bubble: per facet at rest, its normal and pseudo-inverse column
  if constexpr (nv == 4){
    const std::uint32_t facets = terms >> S::bubble0;
    if (facets){
      slots |= 1u << S::bubble;
      for (int d = 0; d < D; ++d) o.c[S::bubble][d] = splat(0.);
    }
    for (std::uint32_t m = facets; m; m &= m - 1, n += 2*D, ++np){
      const int k = __builtin_ctz(m);
      pair t = splat(0.);
      for (int i = 0; i < nv; ++i)
        if (i != k) t += slope(n, i, k);
      const pair q = -2.*t;
      for (int d = 0; d < D; ++d) o.c[S::bubble][d] += q*n[D + d];
      o.s[np] = S::bubble;
      o.dir[np] = n + D;
      o.q[np] = q;
    }
  }
  o.slots = slots;
  o.n = np;
}

// Each slot's monomial as quartic monomials: added, and subtracted (-1: none)
template<int NV>
constexpr std::array<std::array<int, 2>, Shape<NV>::n_slots> slot_monomials(){
  using S = Shape<NV>;
  constexpr Degree<NV, 4> d4{};
  constexpr auto edges = split_eval::lex_edges<NV>();
  std::array<std::array<int, 2>, S::n_slots> t{};
  for (int s = 0; s < S::n_slots; ++s){
    std::array<int, NV> a{}, b{};
    if (s < S::n_ends){
      const int w = edges[std::size_t(s/2)][std::size_t(s%2)], v = edges[std::size_t(s/2)][std::size_t(1 - s%2)];
      a[std::size_t(w)] = 3;
      a[std::size_t(v)] = 1;
      b[std::size_t(w)] = 2;
      b[std::size_t(v)] = 2;
      t[std::size_t(s)][1] = d4.index(b);
    }
    else {
      // the face lacking k times its vertex x, or the bubble
      const int f = s - S::face0, k = s == S::bubble ? NV : NV == 3 ? 3 : f/3, i = NV == 3 ? f : f % 3;
      for (int m = 0; m < NV; ++m) a[std::size_t(m)] = m != k ? 1 : 0;
      if (s < S::bubble) ++a[std::size_t(i < k ? i : i + 1)];
      t[std::size_t(s)][1] = -1;
    }
    t[std::size_t(s)][0] = d4.index(a);
  }
  return t;
}

// The wall terms of two fields as a homogeneous quartic, D pairs a monomial
template<int NV>
inline void term_form(const Terms<NV>& o, pair* w)
{
  constexpr int D = NV - 1;
  static constexpr auto sm = slot_monomials<NV>();
  for (int k = 0; k < Forms<NV>::n4*D; ++k) w[k] = splat(0.);
  for (std::uint32_t m = o.slots; m; m &= m - 1){
    const int s = __builtin_ctz(m);
    const int a = sm[std::size_t(s)][0], b = sm[std::size_t(s)][1];
    for (int x = 0; x < D; ++x) w[a*D + x] += o.c[s][x];
    if (b >= 0)
      for (int x = 0; x < D; ++x) w[b*D + x] -= o.c[s][x];
  }
}

// P2 of two blocks u and v as a homogeneous quadratic, D pairs a monomial:
// l_i^2 takes u_i, l_i l_k 4 u_ik - u_i - u_k
template<typename Cell>
inline void quadratic(const double* u, const double* v, pair* p)
{
  constexpr int nv = Cell::n_verts, D = nv - 1;
  constexpr std::size_t n2 = Cell::n_dofs_max;
  static constexpr Forms<nv> fo{};
  static constexpr auto mid = mid_slots<Cell>();
  const auto at = [&](const int x, const int q){
    const std::size_t i = std::size_t(x)*n2 + std::size_t(q);
    return pair{u[i], v[i]};
  };
  for (int k = 0; k < Forms<nv>::n2; ++k){
    const int a = fo.lo2[std::size_t(k)], b = fo.var2[std::size_t(k)];
    for (int x = 0; x < D; ++x)
      p[k*D + x] = a == b ? at(x, a) : 4.*at(x, mid[std::size_t(a)][std::size_t(b)]) - at(x, a) - at(x, b);
  }
}

// The table's reduced inputs z of two fields: from their vertex divergences
// du, dv and the wall terms' parts, each part's direction in edge
// coordinates
template<typename Cell>
__attribute__((noinline)) void reduce(const Frame<Cell>& f, const Table t, const double* du, const double* dv,
                   const Terms<Cell::n_verts>& o, double* __restrict zu, double* __restrict zv)
{
  constexpr int nv = Cell::n_verts, D = nv - 1;
  const int r = t.r;
  pair g[Terms<nv>::n_pairs][D], qu[Terms<nv>::n_pairs], qv[Terms<nv>::n_pairs];
  const double* V[Terms<nv>::n_pairs];
  for (int k = 0; k < o.n; ++k){
    if (k > 0 && o.dir[k] == o.dir[k - 1])
      for (int y = 0; y < D; ++y) g[k][y] = g[k - 1][y];
    else
      for (int y = 0; y < D; ++y){
        double p = 0.;
        for (int x = 0; x < D; ++x) p += f.G[y + 1][x]*o.dir[k][x];
        g[k][y] = splat(p);
      }
    qu[k] = splat(o.q[k][0]);
    qv[k] = splat(o.q[k][1]);
    V[k] = t.p + std::size_t((nv + o.s[k]*D)*r);
  }
  for (int c = 0; c < r; c += block){
    pair p0 = splat(0.), p1 = splat(0.), q0 = splat(0.), q1 = splat(0.);
    for (int i = 0; i < nv; ++i){
      const pair a = load(t.p + i*r + c), b = load(t.p + i*r + c + 2);
      p0 += a*du[i];
      p1 += b*du[i];
      q0 += a*dv[i];
      q1 += b*dv[i];
    }
    for (int k = 0; k < o.n; ++k){
      pair a = load(V[k] + c)*g[k][0], b = load(V[k] + c + 2)*g[k][0];
      for (int y = 1; y < D; ++y){
        a += load(V[k] + y*r + c)*g[k][y];
        b += load(V[k] + y*r + c + 2)*g[k][y];
      }
      p0 += qu[k]*a;
      p1 += qu[k]*b;
      q0 += qv[k]*a;
      q1 += qv[k]*b;
    }
    store(zu + c, p0);
    store(zu + c + 2, p1);
    store(zv + c, q0);
    store(zv + c + 2, q1);
  }
}

// Sub-cell j's Q_j of two fields from their reduced inputs in its basis, D
// pairs a cubic monomial, a field a lane
template<int NV>
inline void sub_coefficients(const Table t, const int j, const double* zu, const double* zv, pair* y)
{
  constexpr int D = NV - 1;
  const int r = t.r;
  const std::int16_t* ix = t.idx + j*idx_stride<NV>;
  const int rank = ix[0], n = ix[1];
  const std::int16_t* piv = ix + 2;
  const std::int16_t* rows = piv + rank;
  for (int c = 0; c < rank; ++c) y[piv[c]] = pair{zu[c], zv[c]};
  const double* U = t.p + std::size_t(NV)*table_inputs<NV>(r) + std::size_t(j*Forms<NV>::n3*D*r);
  for (int g = 0; g < n; g += 4, U += 4*r){
    pair p[4], q[4];
    for (int i = 0; i < 4; ++i) p[i] = q[i] = splat(0.);
    for (int c = 0; c < r; c += 2){
      const pair bu = load(zu + c), bv = load(zv + c);
      for (int i = 0; i < 4; ++i){
        const pair w = load(U + i*r + c);
        p[i] += w*bu;
        q[i] += w*bv;
      }
    }
    for (int i = 0; i < 4; ++i) y[rows[g + i]] = sums(p[i], q[i]);
  }
}

// Sub-cell j's inputs table
template<int NV>
inline Table sub_inputs(const Table t, const int j){
  return {t.p + std::size_t(j)*table_inputs<NV>(t.r), t.r, t.idx};
}

// The wall terms o as a homogeneous quartic c plus, in sub-cell j, lambda_j E
// Q_j from its coefficients y (none: no correction)
template<typename Cell>
inline void quartic(const Frame<Cell>& f, const Terms<Cell::n_verts>& o, const int j, const pair* y, pair* c)
{
  constexpr int nv = Cell::n_verts, D = nv - 1;
  static constexpr Forms<nv> fo{};
  term_form<nv>(o, c);
  if (!y) return;
#pragma GCC unroll 20
  for (int k = 0; k < Forms<nv>::n3; ++k){
    pair* q = c + fo.up3[std::size_t(j)][std::size_t(k)]*D;
    for (int x = 0; x < D; ++x){
      pair p = f.E[x][0]*y[k*D];
      for (int z = 1; z < D; ++z) p += f.E[x][z]*y[k*D + z];
      q[x] += p;
    }
  }
}

// The quadratic p plus the quartic c of two fields at the barycentrics l,
// added to their values u, v and gradients gu, gv: the values by Euler's sum
// l_i d_i of each form, the gradients sum_{i>0} (d_i - d_0) G_i
template<typename Cell>
inline void add_form(const Frame<Cell>& f, const std::array<double, 4>& l, const pair* p, const pair* c,
                     Vector3d& u, Vector3d& v, Matrix3d& gu, Matrix3d& gv)
{
  constexpr int nv = Cell::n_verts, D = nv - 1;
  using F = Forms<nv>;
  static constexpr F fo{};
  pair m1[nv], m2[F::n2], m3[F::n3];
  for (int i = 0; i < nv; ++i) m1[i] = splat(l[std::size_t(i)]);
#pragma GCC unroll 10
  for (int k = 0; k < F::n2; ++k) m2[k] = m1[fo.lo2[std::size_t(k)]]*m1[fo.var2[std::size_t(k)]];
#pragma GCC unroll 20
  for (int k = 0; k < F::n3; ++k) m3[k] = m2[fo.lo3[std::size_t(k)]]*m1[fo.var3[std::size_t(k)]];
  pair dq[nv][D], dp[nv][D];
#pragma GCC unroll 4
  for (int i = 0; i < nv; ++i){
    for (int x = 0; x < D; ++x) dq[i][x] = dp[i][x] = splat(0.);
#pragma GCC unroll 20
    for (int k = 0; k < F::n3; ++k){
      const pair m = double(F::d3.e[std::size_t(k)][std::size_t(i)] + 1)*m3[k];
      const pair* q = c + fo.up3[std::size_t(i)][std::size_t(k)]*D;
      for (int x = 0; x < D; ++x) dq[i][x] += m*q[x];
    }
    for (int k = 0; k < nv; ++k){
      const pair m = (k == i ? 2. : 1.)*m1[k];
      const pair* q = p + fo.up1[std::size_t(i)][std::size_t(k)]*D;
      for (int x = 0; x < D; ++x) dp[i][x] += m*q[x];
    }
  }
  for (int x = 0; x < D; ++x){
    pair a = m1[0]*(0.25*dq[0][x] + 0.5*dp[0][x]);
    for (int i = 1; i < nv; ++i) a += m1[i]*(0.25*dq[i][x] + 0.5*dp[i][x]);
    u[x] += a[0];
    v[x] += a[1];
  }
  for (int x = 0; x < D; ++x){
    const pair d0 = dq[0][x] + dp[0][x];
    pair e[D];
    for (int i = 1; i < nv; ++i) e[i - 1] = dq[i][x] + dp[i][x] - d0;
    for (int y = 0; y < D; ++y){
      pair s = e[0]*f.G[1][y];
      for (int i = 2; i < nv; ++i) s += e[i - 1]*f.G[i][y];
      gu(x, y) += s[0];
      gv(x, y) += s[1];
    }
  }
}

// The values alone of two fields at the barycentrics l, added to u, v: the
// quadratic p, the wall terms o and, in sub-cell j, lambda_j E Q_j from its
// coefficients y (none: no correction)
template<typename Cell>
inline void add_values(const Frame<Cell>& f, const std::array<double, 4>& l, const int j, const pair* p,
                       const Terms<Cell::n_verts>& o, const pair* y, Vector3d& u, Vector3d& v)
{
  constexpr int nv = Cell::n_verts, D = nv - 1;
  using F = Forms<nv>;
  static constexpr F fo{};
  double m2[F::n2];
#pragma GCC unroll 10
  for (int k = 0; k < F::n2; ++k) m2[k] = l[std::size_t(fo.lo2[std::size_t(k)])]*l[std::size_t(fo.var2[std::size_t(k)])];
  pair a[D];
  for (int x = 0; x < D; ++x) a[x] = m2[0]*p[x];
#pragma GCC unroll 10
  for (int k = 1; k < F::n2; ++k)
    for (int x = 0; x < D; ++x) a[x] += m2[k]*p[k*D + x];
  double w[Shape<nv>::n_slots];
  slot_weights<nv>(l.data(), w);
  for (std::uint32_t m = o.slots; m; m &= m - 1){
    const int s = __builtin_ctz(m);
    for (int x = 0; x < D; ++x) a[x] += w[s]*o.c[s][x];
  }
  if (y){
    pair q[D];
    for (int x = 0; x < D; ++x) q[x] = splat(0.);
#pragma GCC unroll 20
    for (int k = 0; k < F::n3; ++k){
      const double m = m2[fo.lo3[std::size_t(k)]]*l[std::size_t(fo.var3[std::size_t(k)])];
      for (int x = 0; x < D; ++x) q[x] += m*y[k*D + x];
    }
    for (int x = 0; x < D; ++x){
      pair s = f.E[x][0]*q[0];
      for (int z = 1; z < D; ++z) s += f.E[x][z]*q[z];
      a[x] += l[std::size_t(j)]*s;
    }
  }
  for (int x = 0; x < D; ++x){
    u[x] += a[x][0];
    v[x] += a[x][1];
  }
}

// The correction alone of two fields in sub-cell j at the barycentrics l,
// added to their values u, v and, with Grad, their gradients gu, gv; t the
// cell's table, du and dv the fields' vertex divergences, terms their wall terms
template<typename Cell, bool Grad>
inline void add_in(const Frame<Cell>& f, const Table t, const int j, const std::array<double, 4>& l,
                   const double* du, const double* dv, const Terms<Cell::n_verts>& terms, Vector3d& u, Vector3d& v,
                   Matrix3d& gu, Matrix3d& gv)
{
  constexpr int nv = Cell::n_verts, D = nv - 1;
  using F = Forms<nv>;
  pair y[F::n3*D], c[F::n4*D], p[F::n2*D];
  for (pair& x : p) x = splat(0.);
  Terms<nv> o;
  o.slots = 0;
  double zu[Shape<nv>::r_max], zv[Shape<nv>::r_max];
  reduce<Cell>(f, sub_inputs<nv>(t, j), du, dv, terms, zu, zv);
  sub_coefficients<nv>(t, j, zu, zv, y);
  if constexpr (Grad){
    quartic<Cell>(f, o, j, y, c);
    add_form<Cell>(f, l, p, c, u, v, gu, gv);
  }
  else add_values<Cell>(f, l, j, p, o, y, u, v);
}

// The same in the sub-cell of the least barycentric
template<typename Cell, bool Grad>
inline void add(const Frame<Cell>& f, const Table t, const std::array<double, 4>& l, const double* du,
                const double* dv, const Terms<Cell::n_verts>& terms, Vector3d& u, Vector3d& v, Matrix3d& gu,
                Matrix3d& gv)
{
  double mu[Cell::n_verts];
  const int j = split_eval::sub_cell<Cell::n_verts>(l, mu);
  add_in<Cell, Grad>(f, t, j, l, du, dv, terms, u, v, gu, gv);
}

// A wall cell made ready for evaluations in it, of two fields: its frame, wall
// terms, vertex divergences, P2 as a quadratic and, as they are asked,
// its sub-cells' Q_j and, for gradients, its sub-cells' quartics of the rest;
// those known by bit: Q_j j, the quartic nv + j
template<typename Cell>
struct Ready {
  static constexpr int nv = Cell::n_verts, D = nv - 1;
  using F = Forms<nv>;
  Frame<Cell> frame;
  Terms<nv> terms;
  double du[nv], dv[nv];
  std::array<pair, F::n2*D> p;
  mutable std::array<pair, nv*F::n3*D> y;
  mutable std::array<pair, nv*F::n4*D> c;
  mutable unsigned known;
};

// A thread's last wall cell made ready, valid while key, what its blocks were read at, holds
template<typename Cell, typename Key>
struct Memo {
  Key key{};
  int id = -1;
  Ready<Cell> ready;
};

// Wall cell c from its P2 blocks u and v
template<typename Cell>
inline void ready(const Cell& cell, const Tables<Cell::n_verts>& t, const WallCell& c, const double* u,
                  const double* v, Ready<Cell>& r)
{
  r.frame = Frame<Cell>(cell);
  term_coefficients<Cell>(t.normals_of(c), c.terms, u, v, r.terms);
  quadratic<Cell>(u, v, r.p.data());
  r.known = 0;
  if (c.table == no_table) return;
  vertex_divergence<Cell>(r.frame, u, r.du);
  vertex_divergence<Cell>(r.frame, v, r.dv);
}

// The whole field of two stamps or fields in ready wall cell c at the
// barycentrics l, in sub-cell j, added to the values u, v and, with Grad,
// gradients gu, gv
template<typename Cell, bool Grad>
inline void add_ready_in(const Tables<Cell::n_verts>& t, const WallCell& c, const Ready<Cell>& r, const int j,
                         const std::array<double, 4>& l, Vector3d& u, Vector3d& v, Matrix3d& gu, Matrix3d& gv)
{
  constexpr int nv = Cell::n_verts, D = nv - 1;
  using F = Forms<nv>;
  pair* y = nullptr;
  if (c.table != no_table){
    y = r.y.data() + j*F::n3*D;
    if (!(r.known >> j & 1)){
      const Table w = t.of(c);
      double zu[Shape<nv>::r_max], zv[Shape<nv>::r_max];
      reduce<Cell>(r.frame, sub_inputs<nv>(w, j), r.du, r.dv, r.terms, zu, zv);
      sub_coefficients<nv>(w, j, zu, zv, y);
      r.known |= 1u << j;
    }
  }
  if constexpr (Grad){
    pair* q = r.c.data() + j*F::n4*D;
    if (!(r.known >> (nv + j) & 1)){
      quartic<Cell>(r.frame, r.terms, j, y, q);
      r.known |= 1u << (nv + j);
    }
    add_form<Cell>(r.frame, l, r.p.data(), q, u, v, gu, gv);
  }
  else add_values<Cell>(r.frame, l, j, r.p.data(), r.terms, y, u, v);
}

// The same in the sub-cell of the least barycentric, any without a table
template<typename Cell, bool Grad>
inline void add_ready(const Tables<Cell::n_verts>& t, const WallCell& c, const Ready<Cell>& r,
                      const std::array<double, 4>& l, Vector3d& u, Vector3d& v, Matrix3d& gu, Matrix3d& gv)
{
  double mu[Cell::n_verts];
  const int j = c.table != no_table ? split_eval::sub_cell<Cell::n_verts>(l, mu) : 0;
  add_ready_in<Cell, Grad>(t, c, r, j, l, u, v, gu, gv);
}

}  // namespace wall_split

#endif
