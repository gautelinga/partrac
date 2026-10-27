#ifndef __NEAR_WALL_HPP
#define __NEAR_WALL_HPP

// The near-wall rules, which raise the velocity's order in the cells with a
// wall vertex where the wall is at rest. A P1 velocity becomes P2 by the edge
// rule: the midpoint of every edge rising from a wall vertex is a fixed linear
// map of the velocity at the edge's fluid end, so u.n grows as delta^2 above
// the facet. The map differs by dimension -- a matrix in 2D,
// u_v/2 + (u_v.n) q in 3D -- so the per-end record and the table are
// specialised. A P2 velocity gains quartic terms that vanish at every P2 node,
// their coefficients linear in the cell's P2 values (UpgradeRule). Nothing here
// needs dolfin: the wall facets are the facet table's walls and their normals
// come from the coordinates.

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include "Tet.hpp"
#include "Triangle.hpp"
#include "typedefs.hpp"

class CellDofs;

namespace near_wall {

// Whether wall_upgrade reads the velocity as P2, else P1
template<typename Cell>
bool velocity_p2(const Uint ncoeffs){ return ncoeffs == Uint(Cell::n_dofs_max); }

// The rest tolerance of a velocity field: round-off of the largest value a cell can read
double rest_tolerance(const std::vector<double>& u);
// The same for a field on the vertex table, where an image vertex's own value
// is never gathered, its master's is
double rest_tolerance_by_vertex(const std::vector<double>& u, const std::vector<std::uint32_t>& vclass,
                                const std::size_t nverts, const Uint dim);

// What an edge end stores, and how many ends a cell has
template<typename Cell> struct WallRule;

template<> struct WallRule<Triangle> {
  // Midpoint of a wall-vertex edge: M u_v
  struct End { double mxx, mxy, myx, myy; };
  // Ends by edge (01, 02, 12), then (first, second)
  static constexpr int n_ends = 6;
  static constexpr std::array<std::array<int, 2>, 3> edge_ends = {{{0, 1}, {0, 2}, {1, 2}}};
};

template<> struct WallRule<Tet> {
  // Midpoint of a wall-vertex edge: u_v/2 + (u_v.n) q
  struct End { double qx, qy, qz, nx, ny, nz; };
  // Ends by edge (01, 02, 03, 12, 13, 23), then (first, second)
  static constexpr int n_ends = 12;
  static constexpr std::array<std::array<int, 2>, 6> edge_ends =
    {{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
};

template<typename Cell>
struct WallEdges {
  std::array<typename WallRule<Cell>::End, WallRule<Cell>::n_ends> ends;
  std::uint8_t wall;   // bit k: vertex k lies on a wall
  WallEdges() {}       // unset: build_wall_edges fills the table in parallel
};

// A block's wall vertices at rest, by bit: those of wall whose values are all
// within tol, n values a direction
template<typename Cell>
inline unsigned rest_vertices(const double* u, const Uint n, const std::uint8_t wall, const double tol)
{
  unsigned rest = 0;
  for (int k = 0; k < Cell::n_verts; ++k){
    bool r = wall >> k & 1;
    for (int d = 0; d < Cell::n_verts - 1; ++d) r = r && std::abs(u[d*n + k]) <= tol;
    if (r) rest |= 1u << k;
  }
  return rest;
}

// The near-wall rule, per dimension: the P2 block of a P1 one whose wall
// vertices are at rest, false where none is
template<typename Cell>
bool wall_block(const double* u, double* u2, const WallEdges<Cell>& w, const double tol);

template<>
inline bool wall_block<Triangle>(const double* u, double* u2, const WallEdges<Triangle>& w,
                          const double tol)
{
  constexpr int n1 = 3;
  // Walls: listed vertices at rest
  unsigned rest = 0;
  for (int k = 0; k < 3; ++k)
    if ((w.wall >> k & 1) && std::abs(u[k]) <= tol && std::abs(u[n1 + k]) <= tol)
      rest |= 1u << k;
  if (rest == 0) return false;

  constexpr auto edge_ends = WallRule<Triangle>::edge_ends;
  for (int k = 0; k < 3; ++k){
    const bool r = rest >> k & 1;
    u2[k] = r ? 0. : u[k];
    u2[6 + k] = r ? 0. : u[n1 + k];
  }
  for (int e = 0; e < 3; ++e){
    const int a = edge_ends[e][0], c = edge_ends[e][1];
    const int m = Triangle::mid_[e];
    const bool wa = rest >> a & 1, wc = rest >> c & 1;
    const WallRule<Triangle>::End& we = w.ends[2*e + (wa ? 0 : 1)];
    if (wa != wc){
      const int v = wa ? c : a;
      const double ux = u[v], uy = u[n1 + v];
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
inline bool wall_block<Tet>(const double* u, double* u2, const WallEdges<Tet>& w,
                     const double tol)
{
  constexpr int n1 = 4;
  // Walls: listed vertices at rest
  unsigned rest = 0;
  for (int k = 0; k < 4; ++k)
    if ((w.wall >> k & 1) && std::abs(u[k]) <= tol && std::abs(u[n1 + k]) <= tol
        && std::abs(u[2*n1 + k]) <= tol)
      rest |= 1u << k;
  if (rest == 0) return false;

  constexpr auto edge_ends = WallRule<Tet>::edge_ends;
  for (int k = 0; k < 4; ++k){
    const bool r = rest >> k & 1;
    for (int c = 0; c < 3; ++c)
      u2[10*c + k] = r ? 0. : u[c*n1 + k];
  }
  for (int e = 0; e < 6; ++e){
    const int a = edge_ends[e][0], b = edge_ends[e][1];
    const int m = Tet::mid_[e];
    const bool wa = rest >> a & 1, wb = rest >> b & 1;
    const WallRule<Tet>::End& we = w.ends[2*e + (wa ? 0 : 1)];
    if (wa != wb){
      const int v = wa ? b : a;
      const double ux = u[v], uy = u[n1 + v], uz = u[2*n1 + v];
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

// wall_block with the wall vertices at rest given, by bit
template<typename Cell>
void raise_block(const double* u, double* u2, const WallEdges<Cell>& w, const unsigned rest);

template<>
inline void raise_block<Triangle>(const double* u, double* u2, const WallEdges<Triangle>& w,
                                  const unsigned rest)
{
  constexpr int n1 = 3;
  constexpr auto edge_ends = WallRule<Triangle>::edge_ends;
  for (int k = 0; k < 3; ++k){
    const bool r = rest >> k & 1;
    u2[k] = r ? 0. : u[k];
    u2[6 + k] = r ? 0. : u[n1 + k];
  }
  for (int e = 0; e < 3; ++e){
    const int a = edge_ends[e][0], c = edge_ends[e][1];
    const int m = Triangle::mid_[e];
    const bool wa = rest >> a & 1, wc = rest >> c & 1;
    const WallRule<Triangle>::End& we = w.ends[2*e + (wa ? 0 : 1)];
    if (wa != wc){
      const int v = wa ? c : a;
      const double ux = u[v], uy = u[n1 + v];
      u2[m] = we.mxx*ux + we.mxy*uy;
      u2[6 + m] = we.myx*ux + we.myy*uy;
    }
    else {
      u2[m] = 0.5*(u2[a] + u2[c]);
      u2[6 + m] = 0.5*(u2[6 + a] + u2[6 + c]);
    }
  }
}

template<>
inline void raise_block<Tet>(const double* u, double* u2, const WallEdges<Tet>& w, const unsigned rest)
{
  constexpr int n1 = 4;
  constexpr auto edge_ends = WallRule<Tet>::edge_ends;
  for (int k = 0; k < 4; ++k){
    const bool r = rest >> k & 1;
    for (int c = 0; c < 3; ++c)
      u2[10*c + k] = r ? 0. : u[c*n1 + k];
  }
  for (int e = 0; e < 6; ++e){
    const int a = edge_ends[e][0], b = edge_ends[e][1];
    const int m = Tet::mid_[e];
    const bool wa = rest >> a & 1, wb = rest >> b & 1;
    const WallRule<Tet>::End& we = w.ends[2*e + (wa ? 0 : 1)];
    if (wa != wb){
      const int v = wa ? b : a;
      const double ux = u[v], uy = u[n1 + v], uz = u[2*n1 + v];
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
}

// Entry i of a table
template<typename A>
constexpr auto& at(A& a, const int i){ return a[std::size_t(i)]; }

// An edge's vertex pairs, in WallRule's order: edge_of[a][b], -1 for a == b
template<typename Cell>
constexpr auto edge_table(){
  constexpr int nv = Cell::n_verts;
  std::array<std::array<int, nv>, nv> t{};
  for (int a = 0; a < nv; ++a)
    for (int b = 0; b < nv; ++b) t[std::size_t(a)][std::size_t(b)] = -1;
  for (std::size_t e = 0; e < WallRule<Cell>::edge_ends.size(); ++e){
    const int a = WallRule<Cell>::edge_ends[e][0], b = WallRule<Cell>::edge_ends[e][1];
    t[std::size_t(a)][std::size_t(b)] = t[std::size_t(b)][std::size_t(a)] = int(e);
  }
  return t;
}

// The vertices off each edge, ascending
template<typename Cell>
constexpr auto apex_table(){
  constexpr int nv = Cell::n_verts;
  std::array<std::array<int, nv - 2>, WallRule<Cell>::edge_ends.size()> t{};
  for (std::size_t e = 0; e < t.size(); ++e){
    std::size_t j = 0;
    for (int v = 0; v < nv; ++v)
      if (v != WallRule<Cell>::edge_ends[e][0] && v != WallRule<Cell>::edge_ends[e][1]) t[e][j++] = v;
  }
  return t;
}

// The face terms' face (the vertex it lacks; 2D: nv, the cell) and vertex, by face then vertex
template<typename Cell>
constexpr auto face_table(){
  constexpr int nv = Cell::n_verts, nf = nv == 3 ? 1 : nv;
  std::array<std::array<int, 2>, std::size_t(nf*3)> t{};
  std::size_t j = 0;
  for (int f = 0; f < nf; ++f){
    const int k = nv == 3 ? nv : f;
    for (int x = 0; x < nv; ++x)
      if (x != k) t[j++] = {k, x};
  }
  return t;
}

// The wall upgrade of a P2 velocity: in a cell with a wall vertex at rest, P2
// plus quartic terms that vanish at every P2 node -- per edge end at a wall
// vertex w, lambda_w^2 lambda_v (lambda_w - lambda_v); per face over a wall
// edge, its three barycentrics' product times one of them; in 3D the cell's
// bubble. The terms' coefficients are linear in the cell's P2 values.
template<typename Cell> struct UpgradeRule {
  static constexpr int nv = Cell::n_verts;
  static constexpr int n_edges = nv*(nv - 1)/2, n_apex = nv - 2;
  // Face terms: three a face, the face lacking vertex k first by k; 2D the cell's
  static constexpr int n_faces = nv == 3 ? 3 : 3*nv, n_bubble = nv == 3 ? 0 : 1;
  static constexpr auto edge_of = edge_table<Cell>();
  static constexpr auto apex = apex_table<Cell>();
  static constexpr auto face_term = face_table<Cell>();
  // The slot of face k's term at vertex x
  static constexpr int face_slot(const int k, const int x){ return 3*(k < nv ? k : 0) + (x < k ? x : x - 1); }
};

// The face of a, b and v by the vertex it lacks; 2D: nv, the cell
template<typename Cell>
constexpr int off_vertex(const int a, const int b, const int v){
  return Cell::n_verts == 3 ? 3 : 6 - a - b - v;
}

// Terms a cell: edge ends, faces, bubble
template<typename Cell>
constexpr int n_terms = 2*UpgradeRule<Cell>::n_edges + UpgradeRule<Cell>::n_faces + UpgradeRule<Cell>::n_bubble;

// The edge end at vertex a of the edge to b
template<typename Cell>
constexpr int end_of(const int a, const int b){
  const int e = UpgradeRule<Cell>::edge_of[std::size_t(a)][std::size_t(b)];
  return 2*e + (WallRule<Cell>::edge_ends[std::size_t(e)][0] == a ? 0 : 1);
}

// A P2 wall cell's terms: where its coefficients start in a stamp's table, D
// a term, and which it has by bit (ends, faces, bubble)
struct TermSpan {
  std::uint32_t offset, mask;
};

// A P2 wall cell's terms at the two stamps in play: their coefficients and which terms
struct TermView {
  const double* prev;
  const double* next;
  std::uint32_t mask;
};

// A P2 wall cell. Its normals are consecutive in the table's pool, in the order
// of the bits that name them: the edge ends' by term, the faces' by triple,
// then each bubble facet's normal and its column of the facet rows' pseudo-inverse
struct UpgradeCell {
  std::uint32_t cell;      // the cell's id
  std::uint32_t normals;   // its first normal in the pool
  std::uint16_t triples;   // bit: an edge and apex with a face normal
  std::uint8_t facets;     // bit k: facet k a bubble row
  std::uint8_t wall;       // bit k: vertex k lies on a wall
};

// The P2 wall cells, their terms and their normals
template<typename Cell>
struct UpgradeTable {
  using Normal = std::array<double, Cell::n_verts - 1>;
  std::vector<UpgradeCell> cells;
  std::vector<TermSpan> spans;   // by wall cell
  std::vector<Normal> normals;
  std::uint32_t n_coef = 0;      // a stamp's coefficients
};

// The coefficients of the cell's terms, D a term at its index, from its P2
// block, its wall vertices at rest by bit
template<typename Cell>
void upgrade_coefficients(const double* u, double* c, const UpgradeCell& w, const std::uint32_t terms,
                          const typename UpgradeTable<Cell>::Normal* n, const unsigned rest);

// A P2 velocity field's terms in every wall cell, from its values by node (D
// a node): the cell's terms only, in their order, zero in a cell with no wall
// vertex at rest
template<typename Cell>
void field_terms(const UpgradeTable<Cell>& upgrade, const CellDofs& u_dofs, const double* u,
                 const double tol, std::vector<double>& coef);

// The barycentric products the terms are made of: each face's, and the cell's
template<typename Cell>
struct TermProducts {
  std::array<double, 4> face;   // the face lacking vertex k; 2D: the cell's at 3
  double cell;
  TermProducts(const std::array<double, 4>& l){
    if constexpr (Cell::n_verts == 3){
      face = {0., 0., 0., l[0]*l[1]*l[2]};
      cell = 0.;
    }
    else {
      const double l01 = l[0]*l[1], l23 = l[2]*l[3];
      face[0] = l[1]*l23;
      face[1] = l[0]*l23;
      face[2] = l01*l[3];
      face[3] = l01*l[2];
      cell = l01*l23;
    }
  }
};

// Every term at the barycentrics l
template<typename Cell>
inline void term_weights(const std::array<double, 4>& l, const TermProducts<Cell>& P, double* W)
{
  using R = UpgradeRule<Cell>;
  constexpr auto edge_ends = WallRule<Cell>::edge_ends;
  constexpr int f0 = 2*R::n_edges;
  for (int e = 0; e < R::n_edges; ++e){
    const double la = at(l, at(edge_ends, e)[0]), lb = at(l, at(edge_ends, e)[1]);
    const double ab = la*lb;
    W[2*e] = la*ab*(la - lb);
    W[2*e + 1] = lb*ab*(lb - la);
  }
  for (int t = 0; t < R::n_faces; ++t)
    W[f0 + t] = at(P.face, at(R::face_term, t)[0])*at(l, at(R::face_term, t)[1]);
  if constexpr (R::n_bubble > 0) W[f0 + R::n_faces] = P.cell;
}

// Term t's gradient, from the barycentrics' gradients g
template<typename Cell>
inline Vector3d term_gradient(const int t, const std::array<double, 4>& l, const TermProducts<Cell>& P,
                              const std::array<Vector3d, Cell::n_verts>& g)
{
  constexpr int nv = Cell::n_verts;
  using R = UpgradeRule<Cell>;
  constexpr auto edge_ends = WallRule<Cell>::edge_ends;
  constexpr int f0 = 2*R::n_edges, b0 = f0 + R::n_faces;
  if (t < f0){
    const int w = at(at(edge_ends, t/2), t%2), v = at(at(edge_ends, t/2), 1 - t%2);
    const double lw = at(l, w), lv = at(l, v);
    return (3.*lw*lw*lv - 2.*lw*lv*lv)*at(g, w) + (lw*lw*lw - 2.*lw*lw*lv)*at(g, v);
  }
  if (t < b0){
    // B l_x: dB/dl_m l_x for m on the face, and B at x
    const int k = at(R::face_term, t - f0)[0], x = at(R::face_term, t - f0)[1];
    int f[3], q = 0;
    for (int j = 0; j < nv; ++j) if (j != k) f[q++] = j;
    const double lx = at(l, x);
    const double l0 = at(l, f[0]), l1 = at(l, f[1]), l2 = at(l, f[2]);
    return (l1*l2*lx)*at(g, f[0]) + (l0*l2*lx)*at(g, f[1]) + (l0*l1*lx)*at(g, f[2]) + at(P.face, k)*at(g, x);
  }
  Vector3d out = Vector3d::Zero();
  if constexpr (R::n_bubble > 0)
    for (int m = 0; m < nv; ++m) out += at(P.face, m)*at(g, m);
  return out;
}

// Two fields' terms at a point, added to their values u and v: D coefficients
// a term at its index, W every term's weight
template<typename Cell>
inline void add_terms(const double* W, const std::uint32_t mask, const double* cu, const double* cv,
                      Vector3d& u, Vector3d& v)
{
  constexpr int D = Cell::n_verts - 1;
  for (std::uint32_t m = mask; m; m &= m - 1){
    const int t = __builtin_ctz(m);
    for (int d = 0; d < D; ++d){
      u[d] += W[t]*cu[D*t + d];
      v[d] += W[t]*cv[D*t + d];
    }
  }
}

// The same to their gradients, g the barycentrics' gradients
template<typename Cell>
inline void add_term_gradients(const std::array<double, 4>& l, const TermProducts<Cell>& P,
                               const std::array<Vector3d, Cell::n_verts>& g, const std::uint32_t mask,
                               const double* cu, const double* cv, Matrix3d& gu, Matrix3d& gv)
{
  constexpr int D = Cell::n_verts - 1;
  for (std::uint32_t m = mask; m; m &= m - 1){
    const int t = __builtin_ctz(m);
    const Vector3d gk = term_gradient<Cell>(t, l, P, g);
    for (int d = 0; d < D; ++d)
      for (int x = 0; x < D; ++x){
        gu(d, x) += cu[D*t + d]*gk[x];
        gv(d, x) += cv[D*t + d]*gk[x];
      }
  }
}

// The upgrade's wall cells and the normals of their terms, by the rules of
// build_wall_edges; a facet at a corner vertex carries no bubble
template<typename Cell>
void build_wall_terms(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                      const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                      const std::vector<std::int32_t>& facet_neigh,
                      const std::vector<std::uint32_t>& vclass,
                      std::vector<std::int32_t>& wall_index,
                      UpgradeTable<Cell>& table, const bool verbose);

// The wall cells and their edge ends, from the mesh arrays and the facet table.
// A wall facet is a facet_wall entry, and its normal is the outward one the
// vertex it faces fixes the sign of; the facets are visited in the order of
// their sorted vertex ids, which is the order the entity numbering gives them,
// so no sum depends on how the cells are numbered. vclass sends a vertex to
// the master of its periodic images; verbose prints the wall cell count.
template<typename Cell>
void build_wall_edges(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                      const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                      const std::vector<std::int32_t>& facet_neigh,
                      const std::vector<std::uint32_t>& vclass,
                      std::vector<std::int32_t>& wall_index,
                      std::vector<WallEdges<Cell>>& wall_cells, const bool verbose);

// The wall tables: vclass where the format left it empty, each cell's wall
// index (-1: no wall vertex) and, for a P2 velocity, the terms' table, else
// the P1 edge ends
template<typename Cell>
void build_tables(const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                  const std::size_t ncells, const std::size_t nverts, const Uint gdim,
                  const std::vector<bool>& periodic, const Vector3d& x_min, const Vector3d& x_max,
                  const double periodic_tol, const std::vector<std::int32_t>& facet_neigh,
                  const bool p2, std::vector<std::uint32_t>& vclass,
                  std::vector<std::int32_t>& wall_index, std::vector<WallEdges<Cell>>& wall_cells,
                  UpgradeTable<Cell>& upgrade, const bool verbose);

// 0 bulk, 1 on a wall, 2 next to a cell on a wall; from the facet table
void label_cell_type(std::vector<int>& cell_type, const std::vector<std::int32_t>& across,
                     const Uint nv);

}  // namespace near_wall

#endif
