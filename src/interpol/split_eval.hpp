#ifndef __SPLIT_EVAL_HPP
#define __SPLIT_EVAL_HPP

// The spatial kernel of the divergence-free split field: on each cell the
// unique continuous P2 field on its barycentric (Alfeld) split that carries the
// cell's P2 boundary data and has div u = 0 pointwise. It is one fixed
// reference matrix through a Piola map, so this holds the matrix, a cell's
// interior values through it, the net flux a cell's data must have, and the
// evaluation in the sub-cell of the smallest barycentric coordinate. Free
// functions over one array of node values, so a loader of stamps and one of
// frequency components share them. The tables come from the cell class's own
// P2 basis, so they cannot drift from quadbasis.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <eigen3/Eigen/Dense>

#include "typedefs.hpp"

namespace split_eval {

// Vertex pairs of a simplex's edges, lexicographic: (0,1), (0,2), ...
template<int NV>
constexpr std::array<std::array<int, 2>, NV*(NV-1)/2> lex_edges(){
  std::array<std::array<int, 2>, NV*(NV-1)/2> e{};
  int k = 0;
  for (int a = 0; a < NV; ++a)
    for (int b = a + 1; b < NV; ++b)
      e[k++] = {a, b};
  return e;
}

// Lexicographic index of the edge (a, b), a < b
template<int NV>
constexpr int lex_index(const int a, const int b){
  int k = 0;
  for (int i = 0; i < NV; ++i)
    for (int j = i + 1; j < NV; ++j, ++k)
      if (i == a && j == b) return k;
  return -1;
}

// The endpoints of the edge each P2 slot holds; a vertex slot gives itself twice
template<typename Cell>
constexpr std::array<std::array<int, 2>, Cell::n_dofs_max> slot_edges(){
  constexpr int nv = Cell::n_verts;
  std::array<std::array<int, 2>, Cell::n_dofs_max> s{};
  for (int k = 0; k < nv; ++k) s[k] = {k, k};
  const auto e = lex_edges<nv>();
  for (std::size_t p = 0; p < e.size(); ++p) s[std::size_t(Cell::mid_[p])] = e[p];
  return s;
}

// The interior nodes of the split: the barycenter, then the midpoint of the
// segment from each macro vertex to it
template<typename Cell>
struct Split {
  static constexpr int D = Cell::n_verts - 1;
  static constexpr int nv = Cell::n_verts;
  static constexpr std::size_t n_bnd = Cell::n_dofs_max;     // the macro P2 nodes
  static constexpr std::size_t n_int = std::size_t(nv) + 1;
  static constexpr std::size_t n_sub = Cell::n_dofs_max;     // a sub-cell's P2 nodes
  static constexpr std::size_t n_ref = n_int*std::size_t(D)*n_bnd*std::size_t(D);
  using RefMatrix = std::array<double, n_ref>;               // row-major
};

// A sub-cell's P2 nodes: the macro slot each one is, or -(1 + interior node).
// Sub-cell i leaves out macro vertex i, so its own vertices are the others in
// order and then the barycenter.
template<typename Cell>
constexpr std::array<std::array<int, Cell::n_dofs_max>, Cell::n_verts> sub_nodes(){
  constexpr int nv = Cell::n_verts;
  constexpr int d = nv - 1;
  const auto se = slot_edges<Cell>();
  std::array<std::array<int, Cell::n_dofs_max>, nv> t{};
  for (int i = 0; i < nv; ++i){
    int f[d] = {};
    int m = 0;
    for (int j = 0; j < nv; ++j) if (j != i) f[m++] = j;
    for (std::size_t k = 0; k < Cell::n_dofs_max; ++k){
      const int a = se[k][0], b = se[k][1];
      if (a == b)             t[i][k] = a < d ? f[a] : -1;
      else if (b == d)        t[i][k] = -(2 + f[a]);
      else                    t[i][k] = Cell::mid_[std::size_t(lex_index<nv>(f[a], f[b]))];
    }
  }
  return t;
}

// The three P2 slots whose values give a facet's mean velocity: in 2D the
// facet's two ends and its midpoint, in 3D the three edge midpoints
template<typename Cell>
constexpr std::array<std::array<int, 3>, Cell::n_verts> facet_nodes(){
  constexpr int nv = Cell::n_verts;
  std::array<std::array<int, 3>, nv> t{};
  for (int o = 0; o < nv; ++o){
    int f[nv - 1] = {};
    int m = 0;
    for (int j = 0; j < nv; ++j) if (j != o) f[m++] = j;
    if constexpr (nv == 3)
      t[o] = {f[0], Cell::mid_[std::size_t(lex_index<nv>(f[0], f[1]))], f[1]};
    else
      t[o] = {Cell::mid_[std::size_t(lex_index<nv>(f[0], f[1]))],
              Cell::mid_[std::size_t(lex_index<nv>(f[0], f[2]))],
              Cell::mid_[std::size_t(lex_index<nv>(f[1], f[2]))]};
  }
  return t;
}

// Their weights: Simpson's rule along an edge, the mean of the midpoints on a face
template<typename Cell>
constexpr std::array<double, 3> facet_weights(){
  if constexpr (Cell::n_verts == 3) return {1./6., 4./6., 1./6.};
  else                              return {1./3., 1./3., 1./3.};
}

// u_int = R g on the reference simplex, from the cell class's own P2 basis: the
// unique P2 field on the Alfeld split with boundary data g and div u = 0 at the
// d+1 vertices of every sub-cell, where a P2 field's divergence is P1. rank is
// the number of interior dofs the constraints pin, n_int*D when it is unique.
template<typename Cell>
typename Split<Cell>::RefMatrix reference_matrix(int& rank)
{
  using S = Split<Cell>;
  constexpr int D = S::D;
  constexpr int nv = S::nv;
  constexpr int nrows = nv*nv;
  constexpr int na = int(S::n_int)*D;
  constexpr int nb = int(S::n_bnd)*D;
  static constexpr auto sn = sub_nodes<Cell>();

  // The reference simplex, and the barycenter its split adds
  double V[nv][3] = {};
  for (int k = 1; k < nv; ++k) V[k][k - 1] = 1.;
  double z[3] = {};
  for (int k = 0; k < nv; ++k)
    for (int c = 0; c < D; ++c) z[c] += V[k][c]/double(nv);

  Eigen::Matrix<double, nrows, na> A = Eigen::Matrix<double, nrows, na>::Zero();
  Eigen::Matrix<double, nrows, nb> B = Eigen::Matrix<double, nrows, nb>::Zero();
  for (int o = 0; o < nv; ++o){
    const double* s[nv] = {};
    int m = 0;
    for (int j = 0; j < nv; ++j) if (j != o) s[m++] = V[j];
    s[nv - 1] = z;
    const Cell sub = [&]{
      if constexpr (nv == 3) return Cell(s[0], s[1], s[2]);
      else                   return Cell(s[0], s[1], s[2], s[3]);
    }();
    // div u is P1 on the sub-cell, so it vanishes exactly when it does at the
    // sub-cell's own vertices
    for (int q = 0; q < nv; ++q){
      double bary[4] = {0., 0., 0., 0.};
      bary[q] = 1.;
      std::array<double, Cell::n_dofs_max> Nx, Ny, Nz;
      if constexpr (nv == 3) sub.quadderiv(bary[0], bary[1], bary[2], Nx.data(), Ny.data());
      else                   sub.quadderiv(bary[0], bary[1], bary[2], bary[3],
                                           Nx.data(), Ny.data(), Nz.data());
      const double* dN[3] = {Nx.data(), Ny.data(), Nz.data()};
      const int row = o*nv + q;
      for (std::size_t k = 0; k < S::n_sub; ++k){
        const int t = sn[std::size_t(o)][k];
        for (int c = 0; c < D; ++c){
          if (t >= 0) B(row, t*D + c) += dN[c][k];
          else        A(row, (-1 - t)*D + c) += dN[c][k];
        }
      }
    }
  }

  const Eigen::ColPivHouseholderQR<Eigen::Matrix<double, nrows, na>> qr(A);
  rank = int(qr.rank());
  const Eigen::Matrix<double, na, nb> X = qr.solve(-B);
  typename S::RefMatrix R{};
  for (int i = 0; i < na; ++i)
    for (int j = 0; j < nb; ++j) R[std::size_t(i)*std::size_t(nb) + std::size_t(j)] = X(i, j);
  return R;
}

// One cell's interior values: g the macro P2 node values, D doubles a node, J
// the cell's affine Jacobian and Jinv its inverse, both row-major D x D. The
// map u = J uhat leaves the divergence in the reference cell, so the reference
// matrix serves every cell.
template<typename Cell>
inline void interior_values(const typename Split<Cell>::RefMatrix& R,
                            const double* J, const double* Jinv,
                            const double* g, double* u_int)
{
  using S = Split<Cell>;
  constexpr int D = S::D;
  constexpr std::size_t nb = S::n_bnd*std::size_t(D);
  constexpr std::size_t ni = S::n_int*std::size_t(D);
  double gh[nb], uh[ni];
  for (std::size_t n = 0; n < S::n_bnd; ++n)
    for (int c = 0; c < D; ++c){
      double v = 0.;
      for (int k = 0; k < D; ++k) v += Jinv[c*D + k]*g[n*std::size_t(D) + std::size_t(k)];
      gh[n*std::size_t(D) + std::size_t(c)] = v;
    }
  for (std::size_t r = 0; r < ni; ++r){
    double v = 0.;
    const double* Rr = R.data() + r*nb;
    for (std::size_t q = 0; q < nb; ++q) v += Rr[q]*gh[q];
    uh[r] = v;
  }
  for (std::size_t n = 0; n < S::n_int; ++n)
    for (int c = 0; c < D; ++c){
      double v = 0.;
      for (int k = 0; k < D; ++k) v += J[c*D + k]*uh[n*std::size_t(D) + std::size_t(k)];
      u_int[n*std::size_t(D) + std::size_t(c)] = v;
    }
}

// A cell's macro P2 node values, D doubles a node, in the basis's slot order
template<typename Cell, int Dim>
inline void gather_bnd(const std::uint32_t* row, const double* u, double* g){
  for (std::size_t k = 0; k < Split<Cell>::n_bnd; ++k){
    const double* v = u + std::size_t(row[k])*std::size_t(Dim);
    for (int c = 0; c < Dim; ++c) g[k*std::size_t(Dim) + std::size_t(c)] = v[c];
  }
}

// A cell's net flux and its largest facet flux, times measure: 2D
// |e|(u_a + 4 u_m + u_b)/6 . n, 3D |f|/3 (the three edge midpoints) . n
template<typename Cell>
inline void flux_parts(const Cell& cell, const double* g, const double measure,
                       double& net, double& big)
{
  constexpr int D = Split<Cell>::D;
  static constexpr auto fn = facet_nodes<Cell>();
  static constexpr auto fw = facet_weights<Cell>();
  net = 0.;
  big = 0.;
  for (int o = 0; o < Cell::n_verts; ++o){
    const Vector3d n = cell.bary_grad(o);
    double q = 0.;
    for (int c = 0; c < D; ++c){
      double v = 0.;
      for (int j = 0; j < 3; ++j)
        v += fw[std::size_t(j)]*g[std::size_t(fn[std::size_t(o)][std::size_t(j)])*D
                                  + std::size_t(c)];
      q -= n[c]*v;
    }
    q *= measure;
    net += q;
    big = std::max(big, std::abs(q));
  }
}

// The macro vertex each sub-local vertex of sub-cell i is; the last is the
// barycenter, which takes i's own slot here
template<int NV>
constexpr std::array<std::array<int, NV>, NV> sub_verts(){
  std::array<std::array<int, NV>, NV> t{};
  for (int i = 0; i < NV; ++i){
    int m = 0;
    for (int j = 0; j < NV; ++j) if (j != i) t[i][m++] = j;
    t[i][NV - 1] = i;
  }
  return t;
}

// The sub-cell of the smallest barycentric coordinate, and the barycentrics mu
// of the point in it: mu_z = NV lambda_i, mu_j = lambda_j - lambda_i
template<int NV>
inline int sub_cell(const std::array<double, 4>& bary, double* mu){
  static constexpr auto sv = sub_verts<NV>();
  // A running minimum in a register, the index by a mask
  double lo = bary[0];
  int i = 0;
  for (int k = 1; k < NV; ++k){
    const double b = bary[k];
    const int m = -int(b < lo);
    i = (k & m) | (i & ~m);
    lo = b < lo ? b : lo;
  }
  for (int m = 0; m < NV - 1; ++m) mu[m] = bary[sv[i][m]] - bary[i];
  mu[NV - 1] = double(NV)*bary[i];
  return i;
}

// The sub-cell's P2 basis in mu, in the cell class's own slot order
template<typename Cell>
inline void sub_basis(const double* mu, double* N){
  static constexpr auto se = slot_edges<Cell>();
  for (int k = 0; k < Cell::n_verts; ++k) N[k] = mu[k]*(2.*mu[k] - 1.);
  for (std::size_t k = Cell::n_verts; k < Cell::n_dofs_max; ++k)
    N[k] = 4.*mu[se[k][0]]*mu[se[k][1]];
}

// Its gradients, by the chain rule on the cell's gradients of lambda
template<typename Cell>
inline void sub_deriv(const Cell& cell, const int i, const double* mu,
                      double* Nx, double* Ny, double* Nz)   // Nz unused in 2D
{
  constexpr int NV = Cell::n_verts;
  static constexpr auto se = slot_edges<Cell>();
  static constexpr auto sv = sub_verts<NV>();
  Vector3d gl[NV];
  for (int k = 0; k < NV; ++k) gl[k] = cell.bary_grad(k);
  Vector3d g[NV];
  const Vector3d gi = gl[i];
  for (int m = 0; m < NV - 1; ++m) g[m] = gl[sv[i][m]] - gi;
  g[NV - 1] = double(NV)*gi;
  for (int k = 0; k < NV; ++k){
    const double a = 4.*mu[k] - 1.;
    Nx[k] = a*g[k][0];
    Ny[k] = a*g[k][1];
    if constexpr (NV == 4) Nz[k] = a*g[k][2];
  }
  for (std::size_t k = NV; k < Cell::n_dofs_max; ++k){
    const int a = se[k][0], b = se[k][1];
    Nx[k] = 4.*(mu[b]*g[a][0] + mu[a]*g[b][0]);
    Ny[k] = 4.*(mu[b]*g[a][1] + mu[a]*g[b][1]);
    if constexpr (NV == 4) Nz[k] = 4.*(mu[b]*g[a][2] + mu[a]*g[b][2]);
  }
  if constexpr (NV == 3) (void)Nz;
}

// One table as three, so both addresses are safe to form: the macro slot (0
// where the node is interior), the interior node (0 where it is a macro one),
// and which of the two it is
template<typename Cell>
constexpr std::array<std::array<std::uint8_t, Cell::n_dofs_max*3>, Cell::n_verts> sub_gather(){
  constexpr auto t = sub_nodes<Cell>();
  constexpr std::size_t N = Split<Cell>::n_sub;
  std::array<std::array<std::uint8_t, Cell::n_dofs_max*3>, Cell::n_verts> g{};
  for (int i = 0; i < Cell::n_verts; ++i)
    for (std::size_t k = 0; k < N; ++k){
      const int s = t[std::size_t(i)][k];
      g[i][k] = std::uint8_t(s >= 0 ? s : 0);
      g[i][N + k] = std::uint8_t(s >= 0 ? 0 : -1 - s);
      g[i][2*N + k] = std::uint8_t(s >= 0 ? 1 : 0);
    }
  return g;
}

// One sub-cell's node values: the macro P2 nodes of the facet opposite i out of
// the field, the rest out of this cell's interior block; a component per
// contiguous run, as block_value reads it
template<typename Cell, int Dim>
inline void gather_sub(const int i, const std::uint32_t* row, const double* u,
                       const double* u_int, double* block){
  static constexpr auto t = sub_gather<Cell>();
  constexpr std::size_t N = Split<Cell>::n_sub;
  const std::uint8_t* g = t[std::size_t(i)].data();
  for (std::size_t k = 0; k < N; ++k){
    const double* v = g[2*N + k] ? u + std::size_t(row[g[k]])*std::size_t(Dim)
                                 : u_int + std::size_t(g[N + k])*std::size_t(Dim);
    for (int c = 0; c < Dim; ++c) block[std::size_t(c)*N + k] = v[c];
  }
}

// The same for the two stamps a step blends between
template<typename Cell, int Dim>
inline void gather_sub_stamps(const int i, const std::uint32_t* row,
                              const double* u_prev, const double* u_next,
                              const double* int_prev, const double* int_next,
                              double* block_prev, double* block_next){
  static constexpr auto t = sub_gather<Cell>();
  constexpr std::size_t N = Split<Cell>::n_sub;
  const std::uint8_t* g = t[std::size_t(i)].data();
  for (std::size_t k = 0; k < N; ++k){
    const bool macro = g[2*N + k] != 0;
    const std::size_t o = macro ? std::size_t(row[g[k]])*std::size_t(Dim)
                                : std::size_t(g[N + k])*std::size_t(Dim);
    const double* p = (macro ? u_prev : int_prev) + o;
    const double* n = (macro ? u_next : int_next) + o;
    for (int c = 0; c < Dim; ++c){
      block_prev[std::size_t(c)*N + k] = p[c];
      block_next[std::size_t(c)*N + k] = n[c];
    }
  }
}


}  // namespace split_eval

#endif
