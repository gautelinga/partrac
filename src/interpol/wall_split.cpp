#include "wall_split.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <eigen3/Eigen/Dense>

namespace wall_split {

namespace {

using Mat = Eigen::MatrixXd;

// Gauss-Legendre on [0, 1], weights summing to 1
void gauss_legendre(const int n, std::vector<double>& x, std::vector<double>& w){
  x.resize(std::size_t(n));
  w.resize(std::size_t(n));
  for (int i = 0; i < n; ++i){
    double z = std::cos(M_PI*(double(i) + 0.75)/(double(n) + 0.5)), dp = 0.;
    for (int it = 0; it < 100; ++it){
      double p0 = 1., p1 = z;
      for (int k = 2; k <= n; ++k){
        const double p2 = (double(2*k - 1)*z*p1 - double(k - 1)*p0)/double(k);
        p0 = p1;
        p1 = p2;
      }
      dp = double(n)*(z*p1 - p0)/(z*z - 1.);
      const double dz = p1/dp;
      z -= dz;
      if (std::abs(dz) < 1e-16) break;
    }
    x[std::size_t(i)] = 0.5*(1. - z);
    w[std::size_t(i)] = 1./((1. - z*z)*dp*dp);
  }
}

// The collapsed Gauss rule on the unit simplex, as barycentrics and weights summing to 1
template<int NV>
void gauss_simplex(const int n, std::vector<std::array<double, NV>>& P, std::vector<double>& W){
  std::vector<double> x, w;
  gauss_legendre(n, x, w);
  P.clear();
  W.clear();
  double sum = 0.;
  for (int a = 0; a < n; ++a)
    for (int b = 0; b < n; ++b)
      for (int c = 0; c < (NV == 4 ? n : 1); ++c){
        const double xa = x[std::size_t(a)], xb = x[std::size_t(b)], xc = x[std::size_t(c)];
        std::array<double, NV> p{};
        double wt;
        if constexpr (NV == 4){
          const double p1 = xa, p2 = xb*(1. - xa), p3 = xc*(1. - xa)*(1. - xb);
          p = {1. - p1 - p2 - p3, p1, p2, p3};
          wt = w[std::size_t(a)]*w[std::size_t(b)]*w[std::size_t(c)]*(1. - xa)*(1. - xa)*(1. - xb);
        }
        else {
          const double p1 = xa, p2 = xb*(1. - xa);
          p = {1. - p1 - p2, p1, p2};
          wt = w[std::size_t(a)]*w[std::size_t(b)]*(1. - xa);
        }
        P.push_back(p);
        W.push_back(wt);
        sum += wt;
      }
  for (double& v : W) v /= sum;
}

// The degree-K Lagrange basis of a sub-cell's lattice at mu: values,
// derivatives and second derivatives in mu
template<int NV>
struct Basis {
  static constexpr int n = Lattice<NV>::n_loc;
  std::array<double, n> phi;
  std::array<std::array<double, NV>, n> d;
  std::array<std::array<std::array<double, NV>, NV>, n> h;
  explicit Basis(const std::array<double, NV>& mu){
    static constexpr Lattice<NV> lat{};
    for (int l = 0; l < n; ++l){
      double f[NV], df[NV], ddf[NV];
      for (int i = 0; i < NV; ++i){
        double v = 1., dv = 0., ddv = 0.;
        for (int m = 0; m < lat.all[std::size_t(l)][std::size_t(i)]; ++m){
          const double t = (double(K)*mu[std::size_t(i)] - double(m))/double(m + 1);
          const double tp = double(K)/double(m + 1);
          ddv = ddv*t + 2.*dv*tp;
          dv = dv*t + v*tp;
          v = v*t;
        }
        f[i] = v;
        df[i] = dv;
        ddf[i] = ddv;
      }
      double p = 1.;
      for (int i = 0; i < NV; ++i) p *= f[i];
      phi[std::size_t(l)] = p;
      for (int i = 0; i < NV; ++i){
        for (int k = 0; k < NV; ++k){
          double q = i == k ? ddf[i] : df[i]*df[k];
          for (int o = 0; o < NV; ++o) if (o != i && o != k) q *= f[o];
          h[std::size_t(l)][std::size_t(i)][std::size_t(k)] = q;
        }
        double q = df[i];
        for (int o = 0; o < NV; ++o) if (o != i) q *= f[o];
        d[std::size_t(l)][std::size_t(i)] = q;
      }
    }
  }
};

// The regular reference simplex, edge 1
template<int NV>
Mat reference_vertices(){
  Mat V(NV, NV - 1);
  if constexpr (NV == 4){
    const double s = 1./std::sqrt(8.);
    V << s, s, s,  s, -s, -s,  -s, s, -s,  -s, -s, s;
  }
  else
    V << 0., 0.,  1., 0.,  0.5, 0.5*std::sqrt(3.);
  return V;
}

// A condition row: sub-cell j, the macro barycentrics, a value (div) or a
// normal (n.grad div) row, and the pattern bit it belongs to
template<int NV>
struct Row {
  int j;
  std::array<double, NV> lam;
  bool normal;
  int bit;
};

// The degree-p lattice of the facet facing k, as macro barycentrics; its centroid for p = 0
template<int NV>
std::vector<std::array<double, NV>> facet_lattice(const int k, const int p){
  int fl[NV - 1], m = 0;
  for (int i = 0; i < NV; ++i) if (i != k) fl[m++] = i;
  std::vector<std::array<double, NV>> out;
  if (p == 0){
    std::array<double, NV> l{};
    for (int i = 0; i < NV - 1; ++i) l[std::size_t(fl[i])] = 1./double(NV - 1);
    out.push_back(l);
    return out;
  }
  int pw = 1;
  for (int i = 0; i < NV - 1; ++i) pw *= p + 1;
  for (int code = 0; code < pw; ++code){
    int b[NV - 1], r = code, s = 0;
    for (int i = NV - 2; i >= 0; --i){ b[i] = r % (p + 1); r /= p + 1; s += b[i]; }
    if (s != p) continue;
    std::array<double, NV> l{};
    for (int i = 0; i < NV - 1; ++i) l[std::size_t(fl[i])] += double(b[i])/double(p);
    out.push_back(l);
  }
  return out;
}

template<int NV>
std::vector<Row<NV>> condition_rows(){
  constexpr auto ev = entity_verts<NV>();
  std::vector<Row<NV>> rows;
  // A codimension-2 entity: div at its degree K-1 lattice in every sub-cell holding it
  for (std::size_t q = 0; q < ev.size(); ++q){
    const int a = ev[q][0], b = ev[q][1];
    for (int j = 0; j < NV; ++j){
      if (j == a || j == b) continue;
      const int np = a == b ? 1 : K;
      for (int s = 0; s < np; ++s){
        const double t = a == b ? 0. : double(s)/double(K - 1);
        std::array<double, NV> l{};
        l[std::size_t(a)] += 1. - t;
        l[std::size_t(b)] += t;
        rows.push_back({j, l, false, NV + int(q)});
      }
    }
  }
  // A facet: div at its degree K-1 lattice, n.grad div at its degree K-2 one, in its sub-cell
  for (int k = 0; k < NV; ++k){
    for (const auto& l : facet_lattice<NV>(k, K - 1)) rows.push_back({k, l, false, k});
    for (const auto& l : facet_lattice<NV>(k, K - 2)) rows.push_back({k, l, true, k});
  }
  return rows;
}

// Sub-cell j's active nodes' Q: lambda_j Q is the node's basis function, Q a
// homogeneous cubic in the cell's barycentrics, its monomials the rows: K NV
// times the factors (K mu_i - m)/(m + 1) but the centroid's first, mu_i =
// l_sv - l_j, mu_c = NV l_j, m as m sum l
template<int NV>
Mat cubic_form(const int j){
  static constexpr Lattice<NV> lat{};
  static constexpr auto sv = split_eval::sub_verts<NV>();
  constexpr Degree<NV, 3> d3{};
  Mat out = Mat::Zero(d3.n, Lattice<NV>::n_act);
  for (int l = 0; l < Lattice<NV>::n_act; ++l){
    const auto& a = lat.act[std::size_t(l)];
    std::map<std::array<int, NV>, double> p{{std::array<int, NV>{}, double(K*NV)}};
    for (int i = 0; i < NV; ++i)
      for (int m = i == NV - 1 ? 1 : 0; m < a[std::size_t(i)]; ++m){
        const double s = 1./double(m + 1);
        double c[NV];
        for (int v = 0; v < NV; ++v) c[v] = -double(m)*s;
        if (i < NV - 1){
          c[sv[std::size_t(j)][std::size_t(i)]] += double(K)*s;
          c[j] -= double(K)*s;
        }
        else c[j] += double(K*NV)*s;
        std::map<std::array<int, NV>, double> q;
        for (const auto& [e, w] : p)
          for (int v = 0; v < NV; ++v){
            std::array<int, NV> f = e;
            ++f[std::size_t(v)];
            q[f] += c[v]*w;
          }
        p = q;
      }
    for (const auto& [e, w] : p) out(d3.index(e), l) = w;
  }
  return out;
}

}  // namespace

template<int NV>
Reference reference()
{
  constexpr int D = NV - 1;
  using S = Shape<NV>;
  constexpr int nl = Lattice<NV>::n_loc;
  static constexpr auto sv = split_eval::sub_verts<NV>();
  const int ni = S::n_int;
  const Mat V = reference_vertices<NV>();
  // Barycentric gradients of the reference simplex, rows
  Mat Ah = Mat::Ones(NV, NV);
  Ah.topRows(D) = V.transpose();
  const Mat G = Ah.inverse().leftCols(D);
  Mat Eh(D, D);
  for (int i = 0; i < D; ++i) Eh.col(i) = (V.row(i + 1) - V.row(0)).transpose();
  double vol = std::abs(Eh.determinant());
  for (int i = 2; i <= D; ++i) vol /= double(i);
  // grad mu of each sub-cell
  std::array<Mat, NV> Gm;
  for (int j = 0; j < NV; ++j){
    Gm[std::size_t(j)] = Mat(NV, D);
    for (int m = 0; m < NV - 1; ++m) Gm[std::size_t(j)].row(m) = G.row(sv[std::size_t(j)][std::size_t(m)]) - G.row(j);
    Gm[std::size_t(j)].row(NV - 1) = double(NV)*G.row(j);
  }
  // Every node of the split, the interior ones first in their order, then
  // those on the macro faces: sub-cell j's local node l is node[j][l]
  std::array<std::array<int, nl>, NV> node;
  std::vector<std::array<double, NV>> lam(static_cast<std::size_t>(ni));
  {
    std::map<std::array<int, NV>, int> face;
    for (int j = 0; j < NV; ++j)
      for (int l = 0, q = 0; l < nl; ++l){
        const auto& a = S::lat.all[std::size_t(l)];
        std::array<int, NV> key{};
        for (int m = 0, s = 0; m < NV; ++m)
          key[std::size_t(m)] = m == j ? a[NV - 1] : NV*a[std::size_t(s++)] + a[NV - 1];
        int id;
        if (a[NV - 1] >= 1) id = S::in.idx[std::size_t(j)][std::size_t(q++)];
        else {
          const auto it = face.emplace(key, int(lam.size())).first;
          id = it->second;
          if (id == int(lam.size())) lam.emplace_back();
        }
        for (int m = 0; m < NV; ++m) lam[std::size_t(id)][std::size_t(m)] = double(key[std::size_t(m)])/double(NV*K);
        node[std::size_t(j)][std::size_t(l)] = id;
      }
  }
  const int nn = int(lam.size());

  // The stiffness of the nodes, a scalar per pair
  Mat St = Mat::Zero(nn, nn);
  {
    std::vector<std::array<double, NV>> qp;
    std::vector<double> qw;
    gauss_simplex<NV>(K + 1, qp, qw);
    for (int j = 0; j < NV; ++j)
      for (std::size_t q = 0; q < qp.size(); ++q){
        const Basis<NV> b(qp[q]);
        Mat g(nl, D);
        for (int l = 0; l < nl; ++l){
          Eigen::RowVectorXd r = Eigen::RowVectorXd::Zero(D);
          for (int m = 0; m < NV; ++m) r += b.d[std::size_t(l)][std::size_t(m)]*Gm[std::size_t(j)].row(m);
          g.row(l) = r;
        }
        for (int l = 0; l < nl; ++l)
          for (int k = 0; k < nl; ++k)
            St(node[std::size_t(j)][std::size_t(l)], node[std::size_t(j)][std::size_t(k)])
              += qw[q]*g.row(l).dot(g.row(k))*vol/double(NV);
      }
  }
  // The rows on the reference simplex over every node, D a node, and their
  // targets per vertex divergence
  const std::vector<Row<NV>> rows = condition_rows<NV>();
  const int nr = int(rows.size());
  Mat B = Mat::Zero(nr, nn*D), T = Mat::Zero(nr, NV);
  for (int r = 0; r < nr; ++r){
    const Row<NV>& row = rows[std::size_t(r)];
    const int j = row.j;
    std::array<double, NV> mu{};
    for (int m = 0; m < NV - 1; ++m)
      mu[std::size_t(m)] = row.lam[std::size_t(sv[std::size_t(j)][std::size_t(m)])] - row.lam[std::size_t(j)];
    mu[NV - 1] = double(NV)*row.lam[std::size_t(j)];
    const Basis<NV> b(mu);
    const Mat& gm = Gm[std::size_t(j)];
    Eigen::RowVectorXd n;
    if (row.normal) n = G.row(row.bit)/G.row(row.bit).norm();
    for (int l = 0; l < nl; ++l){
      Eigen::RowVectorXd c = Eigen::RowVectorXd::Zero(D);
      if (!row.normal)
        for (int m = 0; m < NV; ++m) c += b.d[std::size_t(l)][std::size_t(m)]*gm.row(m);
      else
        for (int m = 0; m < NV; ++m)
          for (int o = 0; o < NV; ++o)
            c += b.h[std::size_t(l)][std::size_t(m)][std::size_t(o)]*gm.row(o).dot(n)*gm.row(m);
      B.block(r, node[std::size_t(j)][std::size_t(l)]*D, 1, D) += c;
    }
    for (int i = 0; i < NV; ++i)
      T(r, i) = row.normal ? -G.row(i).dot(n) : -row.lam[std::size_t(i)];
  }

  Reference out;
  out.S = St;
  out.B = B;
  out.T = T;
  out.Eh = Eh;
  out.n_int = ni;
  for (const Row<NV>& r : rows) out.bit.push_back(r.bit);
  for (const auto& l : lam){
    std::array<double, 4> q{};
    for (int m = 0; m < NV; ++m) q[std::size_t(m)] = l[std::size_t(m)];
    out.lam.push_back(q);
  }
  return out;
}

template<int NV>
std::vector<double> pattern_tables(const std::vector<std::uint16_t>& patterns)
{
  constexpr int D = NV - 1;
  using S = Shape<NV>;
  constexpr int ns = S::n_slots;
  const Reference ref = reference<NV>();
  const Mat& St = ref.S;
  const Mat& B = ref.B;
  const Mat& T = ref.T;
  const Mat& Eh = ref.Eh;
  const Mat Ehi = Eh.inverse();
  const int ni = ref.n_int, nn = int(ref.lam.size()), nr = int(ref.bit.size());
  // S = L L^T on the interior nodes: with y = L^T x a component, the least |y| under the rows
  const Mat Lt = St.topLeftCorner(ni, ni).llt().matrixU();

  // The slots' monomials at the nodes, and the rows of each, per component
  // of its coefficient
  Mat Phi(nn, ns);
  for (int a = 0; a < nn; ++a)
    for (int s = 0; s < ns; ++s){
      double dl[NV];
      Phi(a, s) = slot_weight<NV>(s, ref.lam[std::size_t(a)].data(), dl);
    }
  Mat R = Mat::Zero(nr, ns*D);
  for (int a = 0; a < nn; ++a)
    for (int s = 0; s < ns; ++s)
      for (int y = 0; y < D; ++y) R.col(s*D + y) += Phi(a, s)*B.col(a*D + y);

  // Least reference H1 norm under the active rows: the pseudo-inverse, its
  // singular values cut at sqrt(1e-11) of the largest
  std::vector<double> out(patterns.size()*S::stride);
#pragma omp parallel for schedule(dynamic, 1)
  for (std::size_t p = 0; p < patterns.size(); ++p){
    std::vector<int> act;
    for (int r = 0; r < nr; ++r)
      if (patterns[p] >> ref.bit[std::size_t(r)] & 1) act.push_back(r);
    const int na = int(act.size());
    Mat C(na, ni*D), Ta(na, NV), Ra(na, ns*D);
    for (int r = 0; r < na; ++r){
      Ta.row(r) = T.row(act[std::size_t(r)]);
      Ra.row(r) = R.row(act[std::size_t(r)]);
    }
    for (int x = 0; x < D; ++x){
      Mat Bx(na, ni);
      for (int r = 0; r < na; ++r)
        for (int a = 0; a < ni; ++a) Bx(r, a) = B(act[std::size_t(r)], a*D + x);
      C.middleCols(x*ni, ni) = Lt.transpose().triangularView<Eigen::Lower>().solve(Bx.transpose()).transpose();
    }
    const Eigen::JacobiSVD<Mat> svd(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const Eigen::VectorXd sv = svd.singularValues();
    const double cut = std::sqrt(1e-11)*sv[0];
    Eigen::VectorXd inv(sv.size());
    for (int i = 0; i < sv.size(); ++i) inv[i] = sv[i] > cut ? 1./sv[i] : 0.;
    const Mat Y = svd.matrixV()*(inv.asDiagonal()*svd.matrixU().transpose());
    // The interior values per row target, D a node
    Mat P(ni*D, na);
    for (int x = 0; x < D; ++x){
      const Mat Xx = Lt.triangularView<Eigen::Upper>().solve(Y.middleRows(x*ni, ni));
      for (int a = 0; a < ni; ++a) P.row(a*D + x) = Xx.row(a);
    }
    const Mat Wd = P*Ta;
    const Mat Wc = -P*Ra;
    // In edge coordinates, the slots' coefficients too
    double* o = out.data() + p*S::stride;
    for (int a = 0; a < ni; ++a){
      double* q = o + std::size_t(a)*S::node;
      const Mat w = Ehi*Wd.middleRows(a*D, D);
      for (int x = 0; x < D; ++x)
        for (int i = 0; i < NV; ++i) q[x*NV + i] = w(x, i);
      for (int s = 0; s < ns; ++s){
        const Mat c = Ehi*Wc.block(a*D, s*D, D, D)*Eh;
        for (int x = 0; x < D; ++x)
          for (int y = 0; y < D; ++y) q[S::n_d + std::size_t(s*D*D + x*D + y)] = c(x, y);
      }
    }
  }
  return out;
}

template<int NV>
std::vector<double> factor(const double* w, int& r, std::vector<std::int16_t>& idx)
{
  constexpr int D = NV - 1;
  using S = Shape<NV>;
  constexpr int ni = S::n_int, nc = NV + S::n_slots*D;
  // Per interior node its D rows: the vertex divergences', then each slot's
  Mat A(ni*D, nc);
  for (int a = 0; a < ni; ++a)
    for (int x = 0; x < D; ++x){
      const double* q = w + std::size_t(a)*S::node;
      for (int i = 0; i < NV; ++i) A(a*D + x, i) = q[x*NV + i];
      for (int s = 0; s < S::n_slots; ++s)
        for (int y = 0; y < D; ++y) A(a*D + x, NV + s*D + y) = q[S::n_d + std::size_t(s*D*D + x*D + y)];
    }
  // Its singular values cut at 1e-12 of the largest
  const Eigen::JacobiSVD<Mat> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  const Eigen::VectorXd sv = svd.singularValues();
  int rank = 0;
  while (rank < sv.size() && sv[rank] > 1e-12*sv[0]) ++rank;
  r = (rank + block - 1)/block*block;
  const Mat U = svd.matrixU().leftCols(rank)*sv.head(rank).asDiagonal();
  const Mat V = svd.matrixV().leftCols(rank);
  std::vector<double> out(table_size<NV>(r), 0.);
  idx.assign(std::size_t(NV*idx_stride<NV>), 0);
  constexpr int nq = Forms<NV>::n3*D;
  for (int j = 0; j < NV; ++j){
    // Q_j from its active nodes, a row a cubic monomial and component
    const Mat Mf = cubic_form<NV>(j);
    Mat Q = Mat::Zero(nq, rank);
    for (int k = 0; k < Forms<NV>::n3; ++k)
      for (int x = 0; x < D; ++x)
        for (int l = 0; l < S::n_act; ++l)
          Q.row(k*D + x) += Mf(k, l)*U.row(S::in.idx[std::size_t(j)][std::size_t(l)]*D + x);
    // Its pivot rows T, by column pivoting: in the basis T z, T's rows are z
    // and the others Q T^-1
    std::vector<int> piv(std::size_t(rank), 0), rest;
    std::vector<bool> is(std::size_t(nq), false);
    if (rank > 0){
      const Eigen::ColPivHouseholderQR<Mat> qr(Q.transpose());
      for (int c = 0; c < rank; ++c){
        piv[std::size_t(c)] = qr.colsPermutation().indices()(c);
        is[std::size_t(piv[std::size_t(c)])] = true;
      }
    }
    for (int q = 0; q < nq; ++q) if (!is[std::size_t(q)]) rest.push_back(q);
    Mat T(rank, rank);
    for (int c = 0; c < rank; ++c) T.row(c) = Q.row(piv[std::size_t(c)]);
    const Mat Qt = T.transpose().partialPivLu().solve(Q.transpose()).transpose();
    const Mat Vj = V*T.transpose();
    double* v = out.data() + std::size_t(j)*table_inputs<NV>(r);
    for (int k = 0; k < nc; ++k)
      for (int c = 0; c < rank; ++c) v[std::size_t(k*r + c)] = Vj(k, c);
    const int n = int(rest.size()), np = (n + 3)/4*4;
    double* u = out.data() + std::size_t(NV)*table_inputs<NV>(r) + std::size_t(j*nq*r);
    std::int16_t* ix = idx.data() + j*idx_stride<NV>;
    ix[0] = std::int16_t(rank);
    ix[1] = std::int16_t(np);
    for (int c = 0; c < rank; ++c) ix[2 + c] = std::int16_t(piv[std::size_t(c)]);
    for (int i = 0; i < np; ++i){
      const int q = rest[std::size_t(std::min(i, n - 1))];
      ix[2 + rank + i] = std::int16_t(q);
      for (int c = 0; c < rank; ++c) u[std::size_t(i*r + c)] = Qt(q, c);
    }
  }
  return out;
}

template<int NV>
void Tables<NV>::build(const bool verbose)
{
  std::vector<std::uint8_t> seen(std::size_t(1) << Shape<NV>::n_bits, 0);
  for (const WallCell& c : cell) seen[c.rest & Shape<NV>::pattern_mask] = 1;
  pattern.clear();
  for (std::size_t p = 1; p < seen.size(); ++p)
    if (seen[p]) pattern.push_back(std::uint16_t(p));
  std::vector<std::uint16_t> index(seen.size(), no_table);
  for (std::size_t t = 0; t < pattern.size(); ++t) index[pattern[t]] = std::uint16_t(t);
  const std::vector<double> dense = pattern_tables<NV>(pattern);
  std::vector<std::vector<double>> f(pattern.size());
  std::vector<std::vector<std::int16_t>> fi(pattern.size());
  rank.assign(pattern.size(), 0);
#pragma omp parallel for schedule(dynamic, 1)
  for (std::size_t t = 0; t < pattern.size(); ++t)
    f[t] = factor<NV>(dense.data() + t*Shape<NV>::stride, rank[t], fi[t]);
  idx.clear();
  for (const auto& i : fi) idx.insert(idx.end(), i.begin(), i.end());
  W.clear();
  at.assign(pattern.size(), 0);
  for (std::size_t t = 0; t < pattern.size(); ++t){
    at[t] = W.size();
    W.insert(W.end(), f[t].begin(), f[t].end());
  }
#pragma omp parallel for schedule(static)
  for (std::size_t w = 0; w < cell.size(); ++w) cell[w].table = index[cell[w].rest & Shape<NV>::pattern_mask];
  if (verbose)
    std::cout << "Wall patterns: " << pattern.size() << std::endl;
}

template Reference reference<3>();
template Reference reference<4>();
template std::vector<double> pattern_tables<3>(const std::vector<std::uint16_t>&);
template std::vector<double> pattern_tables<4>(const std::vector<std::uint16_t>&);
template std::vector<double> factor<3>(const double*, int&, std::vector<std::int16_t>&);
template std::vector<double> factor<4>(const double*, int&, std::vector<std::int16_t>&);
template struct Tables<3>;
template struct Tables<4>;

}  // namespace wall_split
