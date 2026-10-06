#include "openfoam_nodes.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>

#include <omp.h>

#include "Error.hpp"

namespace openfoam_nodes {

using openfoam_load::CaseData;
using openfoam_load::FieldData;
using openfoam_split::SplitData;

namespace {

using V3 = std::array<double, 3>;

V3 unit(const V3& a){
  const double n = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  const double s = std::max(n, 1e-300);
  return {a[0]/s, a[1]/s, a[2]/s};
}

// CSR of the (row, value) pairs visit(emit) emits, twice called: each row's
// values ascending and distinct
template<typename Visit>
void csr(const std::size_t nrows, Visit&& visit, std::vector<std::size_t>& start, std::vector<std::size_t>& val){
  std::vector<std::size_t> at(nrows + 1, 0);
  visit([&](const std::size_t r, std::size_t){ ++at[r + 1]; });
  for (std::size_t i = 0; i < nrows; ++i) at[i + 1] += at[i];
  std::vector<std::size_t> v(at[nrows]);
  std::vector<std::size_t> fill(at.begin(), at.end() - 1);
  visit([&](const std::size_t r, const std::size_t x){ v[fill[r]++] = x; });
  start.assign(nrows + 1, 0);
  val.clear();
  val.reserve(v.size());
  for (std::size_t i = 0; i < nrows; ++i){
    const auto b = v.begin() + std::ptrdiff_t(at[i]), e = v.begin() + std::ptrdiff_t(at[i + 1]);
    std::sort(b, e);
    val.insert(val.end(), b, std::unique(b, e));
    start[i + 1] = val.size();
  }
}

// pointConstraint from the normals of the symmetry patches at a point: one
// normal I - nn, two the line along their cross product, three nothing
std::array<double, 9> constraint(const std::vector<V3>& normals){
  int first = 0;
  V3 second{0., 0., 0.};
  for (const V3& n : normals){
    if (first == 0){
      first = 1;
      second = n;
    }
    else if (first == 1){
      const V3 pn{n[1]*second[2] - n[2]*second[1], n[2]*second[0] - n[0]*second[2],
                  n[0]*second[1] - n[1]*second[0]};
      const double m = std::sqrt(pn[0]*pn[0] + pn[1]*pn[1] + pn[2]*pn[2]);
      if (m > 1e-3){
        first = 2;
        second = {pn[0]/m, pn[1]/m, pn[2]/m};
      }
    }
    else if (first == 2){
      if (std::abs(n[0]*second[0] + n[1]*second[1] + n[2]*second[2]) > 1e-3) first = 3;
    }
  }
  std::array<double, 9> T{};
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j){
      const double nn = second[i]*second[j];
      T[3*i + j] = first == 1 ? (i == j ? 1. : 0.) - nn : first == 2 ? nn : 0.;
    }
  return T;
}

// The cyclic master of node n's point, n not a centre
std::size_t master_point(const Geometry& g, const SplitData& s, const std::size_t n){
  return std::size_t(g.master[std::size_t(s.node_point[n])]);
}

// The masters of the nodes' points, ascending; row_of a master's row among
// them, -1 for a point no node needs
std::vector<std::size_t> node_masters(const Geometry& g, const SplitData& s, std::vector<std::int64_t>& row_of){
  row_of.assign(g.c->npoints(), -1);
  for (std::size_t n = 0; n < s.nnodes(); ++n)
    if (s.node_point[n] >= 0) row_of[master_point(g, s, n)] = 0;
  std::vector<std::size_t> masters;
  for (std::size_t p = 0; p < row_of.size(); ++p)
    if (row_of[p] == 0){
      row_of[p] = std::int64_t(masters.size());
      masters.push_back(p);
    }
  return masters;
}

}  // namespace

InverseDistance::InverseDistance(const Geometry& g, const SplitData& s){
  const CaseData& c = *g.c;
  const std::size_t np = c.npoints();
  const std::size_t ni = std::size_t(c.n_internal);
  const std::size_t nb = c.nfaces() - ni;
  ncells_ = std::size_t(c.ncells);
  n_internal_ = ni;
  const std::vector<std::int32_t>& master = g.master;

  // Patch kinds
  std::vector<char> blend(c.patches.size(), 0);
  symmetry_.assign(c.patches.size(), 0);
  for (std::size_t i = 0; i < c.patches.size(); ++i){
    const std::string& t = c.patches[i].type;
    blend[i] = t != "cyclic" && t != "empty";
    symmetry_[i] = t == "symmetry" || t == "symmetryPlane";
    patch_start_.push_back(std::size_t(c.patches[i].start));
  }
  bface_patch_ = g.bface_patch;
  bface_owner_.assign(c.owner.begin() + std::ptrdiff_t(ni), c.owner.end());
  bface_normal_.resize(nb);
  for (std::size_t b = 0; b < nb; ++b)
    bface_normal_[b] = {g.bface_normal[3*b], g.bface_normal[3*b + 1], g.bface_normal[3*b + 2]};
  const auto blends = [&](const std::size_t f){ return blend[std::size_t(g.bface_patch[f - ni])]; };

  // A row a master point: its boundary faces on patches that blend if any
  // image has one, else its cells
  std::vector<std::int64_t> row_of(np, -1);
  start_.push_back(0);
  for (std::size_t n = 0; n < s.nnodes(); ++n){
    if (s.node_point[n] < 0) continue;
    const std::size_t m = master_point(g, s, n);
    if (row_of[m] >= 0) continue;
    row_of[m] = std::int64_t(den_.size());
    bool boundary = false;
    for (std::size_t i = g.img_start[m]; i < g.img_start[m + 1]; ++i){
      const std::size_t q = g.img[i];
      for (std::size_t k = g.pb_start[q]; k < g.pb_start[q + 1]; ++k) boundary = boundary || blends(g.pb[k]);
    }
    double den = 0.;
    for (std::size_t i = g.img_start[m]; i < g.img_start[m + 1]; ++i){
      const std::size_t q = g.img[i];
      const auto& st = boundary ? g.pb_start : g.pc_start;
      const auto& val = boundary ? g.pb : g.pc;
      for (std::size_t k = st[q]; k < st[q + 1]; ++k){
        const std::size_t e = val[k];
        if (boundary && !blends(e)) continue;
        double d2 = 0.;
        for (std::size_t d = 0; d < 3; ++d){
          const double y = boundary ? c.bface_centres[3*(e - ni) + d] : c.cell_centres[3*e + d];
          const double dd = c.points[3*q + d] - y;
          d2 += dd*dd;
        }
        const double w = 1./std::sqrt(d2);
        src_.push_back(boundary ? ncells_ + (e - ni) : e);
        w_.push_back(w);
        den += w;
      }
    }
    den_.push_back(den);
    start_.push_back(src_.size());
  }

  // The symmetry patches' normals at each master, patch by patch
  std::map<std::int32_t, std::vector<V3>> normals;
  for (std::size_t i = 0; i < c.patches.size(); ++i){
    const auto& p = c.patches[i];
    if (!symmetry_[i] || p.size == 0) continue;
    std::map<std::int32_t, V3> acc;
    V3 sum{0., 0., 0.};
    for (std::size_t f = std::size_t(p.start); f < std::size_t(p.start + p.size); ++f){
      const std::size_t b = f - ni;
      for (std::size_t d = 0; d < 3; ++d) sum[d] += c.bface_areas[3*b + d];
      for (auto q = c.face_start[f]; q < c.face_start[f + 1]; ++q){
        V3& a = acc.emplace(master[std::size_t(c.face_points[std::size_t(q)])], V3{0., 0., 0.}).first->second;
        for (std::size_t d = 0; d < 3; ++d) a[d] += bface_normal_[b][d];
      }
    }
    for (const auto& kv : acc)
      normals[kv.first].push_back(p.type == "symmetryPlane" ? unit(sum) : unit(kv.second));
  }
  constraint_.assign(den_.size(), -1);
  for (const auto& kv : normals){
    const std::int64_t r = row_of[std::size_t(kv.first)];
    if (r < 0) continue;
    constraint_[std::size_t(r)] = std::int32_t(t_.size());
    t_.push_back(constraint(kv.second));
  }

  kind_.assign(s.node_kind.begin(), s.node_kind.end());
  ref_.resize(s.nnodes());
  for (std::size_t n = 0; n < s.nnodes(); ++n)
    ref_[n] = s.node_point[n] >= 0 ? std::size_t(row_of[master_point(g, s, n)]) : std::size_t(s.node_cell[n]);
}

void InverseDistance::apply(const FieldData& f, const std::vector<int>& comps,
                            std::vector<double>& out) const {
  const std::size_t k = std::size_t(f.ncomp);
  if (f.internal.size() != k*ncells_)
    partrac::fail(f.name, " at ", f.time, ": ", f.internal.size()/k, " cell values for ", ncells_, " cells");
  // Each boundary face's value as its patch field holds it
  const std::size_t nb = bface_patch_.size();
  std::vector<double> bv(nb*k, 0.);
  for (std::size_t b = 0; b < nb; ++b){
    const std::size_t pi = std::size_t(bface_patch_[b]);
    const auto& fp = f.patches[pi];
    const double* own = f.internal.data() + k*std::size_t(bface_owner_[b]);
    double* v = bv.data() + k*b;
    if (fp.has_value && !symmetry_[pi]){
      const std::size_t i = b + n_internal_ - patch_start_[pi];
      for (std::size_t d = 0; d < k; ++d) v[d] = fp.values[k*i + d];
      continue;
    }
    for (std::size_t d = 0; d < k; ++d) v[d] = own[d];
    if (symmetry_[pi] && k == 3){
      const auto& n = bface_normal_[b];
      const double un = own[0]*n[0] + own[1]*n[1] + own[2]*n[2];
      for (std::size_t d = 0; d < 3; ++d) v[d] = own[d] - un*n[d];
    }
  }
  // The rows
  const std::size_t nrows = den_.size();
  std::vector<double> rv(nrows*k, 0.);
  for (std::size_t r = 0; r < nrows; ++r){
    double num[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
    for (std::size_t q = start_[r]; q < start_[r + 1]; ++q){
      const std::size_t e = src_[q];
      const double* v = e < ncells_ ? f.internal.data() + k*e : bv.data() + k*(e - ncells_);
      for (std::size_t d = 0; d < k; ++d) num[d] += w_[q]*v[d];
    }
    double* o = rv.data() + k*r;
    for (std::size_t d = 0; d < k; ++d) o[d] = num[d]/(den_[r] > 0. ? den_[r] : 1.);
    if (k == 3 && constraint_[r] >= 0){
      const auto& T = t_[std::size_t(constraint_[r])];
      const double a[3] = {o[0], o[1], o[2]};
      for (std::size_t i = 0; i < 3; ++i) o[i] = T[3*i]*a[0] + T[3*i + 1]*a[1] + T[3*i + 2]*a[2];
    }
  }
  const std::size_t nc = comps.size();
  out.assign(kind_.size()*nc, 0.);
  for (std::size_t n = 0; n < kind_.size(); ++n){
    const double* v = (kind_[n] ? f.internal.data() : rv.data()) + k*ref_[n];
    for (std::size_t j = 0; j < nc; ++j) out[n*nc + j] = v[comps[j]];
  }
}

// The least-squares operator

namespace {

// Conditions whose values are never data, whatever fixesValue says
const char* const not_data[] = {"zeroGradient", "inletOutlet", "outletInlet", "calculated",
                                "extrapolatedCalculated", "symmetry", "symmetryPlane", "empty",
                                "cyclic", "slip", "fixedFluxPressure"};

// The last ring tried, and the eigenvalue below which the normal matrix is singular
constexpr int last_ring = 2;
constexpr double singular = 1e-12;

// A cell or a face, and the shift that carries it to its master's side
struct Entry {
  std::size_t id;
  V3 s;
};

// The eigenvalues (ascending) and eigenvectors (columns) of a symmetric n x n matrix, by Jacobi
void eigen(double a[4][4], const int n, double lam[4], double v[4][4]){
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) v[i][j] = i == j;
  double diag = 0.;
  for (int i = 0; i < n; ++i) diag += a[i][i]*a[i][i];
  for (int sweep = 0; sweep < 60; ++sweep){
    double off = 0.;
    for (int p = 0; p < n; ++p)
      for (int q = p + 1; q < n; ++q) off += a[p][q]*a[p][q];
    if (off <= 1e-32*diag || off == 0.) break;
    for (int p = 0; p < n; ++p)
      for (int q = p + 1; q < n; ++q){
        if (a[p][q] == 0.) continue;
        const double th = (a[q][q] - a[p][p])/(2*a[p][q]);
        const double t = (th >= 0 ? 1. : -1.)/(std::abs(th) + std::sqrt(th*th + 1.));
        const double c = 1./std::sqrt(t*t + 1.), s = t*c;
        for (int k = 0; k < n; ++k){
          const double kp = a[k][p], kq = a[k][q];
          a[k][p] = c*kp - s*kq;
          a[k][q] = s*kp + c*kq;
        }
        for (int k = 0; k < n; ++k){
          const double pk = a[p][k], qk = a[q][k];
          a[p][k] = c*pk - s*qk;
          a[q][k] = s*pk + c*qk;
        }
        for (int k = 0; k < n; ++k){
          const double kp = v[k][p], kq = v[k][q];
          v[k][p] = c*kp - s*kq;
          v[k][q] = s*kp + c*kq;
        }
      }
  }
  int order[4] = {0, 1, 2, 3};
  for (int i = 1; i < n; ++i)
    for (int j = i; j > 0 && a[order[j]][order[j]] < a[order[j - 1]][order[j - 1]]; --j) std::swap(order[j], order[j - 1]);
  double w[4][4];
  for (int i = 0; i < n; ++i){
    lam[i] = a[order[i]][order[i]];
    for (int k = 0; k < n; ++k) w[k][i] = v[k][order[i]];
  }
  for (int i = 0; i < n; ++i)
    for (int k = 0; k < n; ++k) v[k][i] = w[k][i];
}

struct Fit {
  bool fail = false;
  double cond = 1.;
};

// Weighted linear least squares u = a + G.q over the offsets q (d a point):
// each point's coefficient of a. Weights (h/|q|)^2, h the mean distance; the
// rank test on the normal matrix scaled to unit diagonal; a singular one by
// its pseudo-inverse; with grad, each point's coefficients of G too (d a point)
Fit wls(const double* q, const std::size_t n, const std::size_t d, const double tol, double* coef,
        double* grad = nullptr){
  const std::size_t m = d + 1;
  double h = 0.;
  for (std::size_t e = 0; e < n; ++e){
    double r2 = 0.;
    for (std::size_t i = 0; i < d; ++i) r2 += q[d*e + i]*q[d*e + i];
    coef[e] = std::sqrt(r2);
    h += coef[e];
  }
  h = n ? h/double(n) : 0.;
  if (!(h > 0.)) h = 1.;
  double N[4][4] = {};
  for (std::size_t e = 0; e < n; ++e){
    const double r = std::max(coef[e], 1e-300);
    const double w = (h/r)*(h/r);
    double A[4] = {1., 0., 0., 0.};
    for (std::size_t i = 0; i < d; ++i) A[i + 1] = q[d*e + i]/h;
    for (std::size_t i = 0; i < m; ++i)
      for (std::size_t j = i; j < m; ++j) N[i][j] += w*A[i]*A[j];
    coef[e] = w;
  }
  double D[4], Ns[4][4];
  for (std::size_t i = 0; i < m; ++i) D[i] = 1./std::sqrt(std::max(N[i][i], 1e-300));
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t j = i; j < m; ++j) Ns[i][j] = Ns[j][i] = N[i][j]*D[i]*D[j];
  double lam[4], V[4][4];
  eigen(Ns, int(m), lam, V);
  const double lmax = std::max(lam[m - 1], 1e-300);
  Fit fit;
  fit.fail = lam[0] < tol*lmax;
  fit.cond = lmax/std::max(lam[0], 1e-300);
  double amax = 0.;
  for (std::size_t i = 0; i < m; ++i) amax = std::max(amax, std::abs(lam[i]));
  // y = D Ns^+ D e_j
  const auto solve = [&](const std::size_t j, double* y){
    for (std::size_t i = 0; i < m; ++i){
      if (lam[0] < singular*lmax && !(std::abs(lam[i]) > singular*amax)) continue;
      const double c = V[j][i]*D[j]/lam[i];
      for (std::size_t k = 0; k < m; ++k) y[k] += c*V[k][i];
    }
    for (std::size_t k = 0; k < m; ++k) y[k] *= D[k];
  };
  double x[4] = {0., 0., 0., 0.};
  solve(0, x);
  // G's component a
  for (std::size_t a = 0; grad && a < d; ++a){
    double y[4] = {0., 0., 0., 0.};
    solve(a + 1, y);
    for (std::size_t e = 0; e < n; ++e){
      double g = y[0];
      for (std::size_t i = 0; i < d; ++i) g += y[i + 1]*q[d*e + i]/h;
      grad[d*e + a] = coef[e]*g/h;
    }
  }
  for (std::size_t e = 0; e < n; ++e){
    double a = x[0];
    for (std::size_t i = 0; i < d; ++i) a += x[i + 1]*q[d*e + i]/h;
    coef[e] *= a;
  }
  return fit;
}

// The row of every point that needs one, from one thread
struct Rows {
  std::vector<std::size_t> count;             // entries a row
  std::vector<std::uint32_t> col;
  std::vector<double> coef;
  std::vector<std::uint32_t> mirror;          // into table, 0 the identity
  std::vector<std::array<double, 9>> table{std::array<double, 9>{1., 0., 0., 0., 1., 0., 0., 0., 1.}};
  std::vector<std::int32_t> zero_wins;
  std::size_t fixed = 0, second_ring = 0, deficient = 0, mirrored = 0;
  bool gradient = false;                      // the fit's G too, into dcoef
  std::vector<double> dcoef;                  // d an entry
};

// What building a row needs, shared by the threads
struct Builder {
  const Geometry& g;
  const CaseData& c;
  const FieldData& f;
  const std::vector<char>& data;
  std::vector<char> symmetry, uniform;
  double tol = 0., rank_tol = 0., quantum = 0.;
  std::size_t k = 1;
  LeastSquaresReport* report = nullptr;

  // Scratch, per thread, kept across rows
  struct Work {
    std::vector<Entry> ring, next, faces, group, patch_ring;
    std::vector<double> q, coef, gcoef, gval, dcoef;
    std::vector<std::size_t> gcol, gid;
    std::vector<char> zero, pick;
    std::vector<std::pair<std::int32_t, V3>> normals;
    std::vector<V3> planes;
    std::vector<std::array<double, 9>> mirrors;
  };

  std::size_t ni() const { return std::size_t(c.n_internal); }
  V3 x(const std::size_t p) const { return {c.points[3*p], c.points[3*p + 1], c.points[3*p + 2]}; }
  V3 centre(const std::size_t cell) const {
    return {c.cell_centres[3*cell], c.cell_centres[3*cell + 1], c.cell_centres[3*cell + 2]};
  }
  V3 fcentre(const std::size_t f) const {
    const std::size_t b = f - ni();
    return {c.bface_centres[3*b], c.bface_centres[3*b + 1], c.bface_centres[3*b + 2]};
  }
  V3 normal(const std::size_t f) const {
    const std::size_t b = f - ni();
    return {g.bface_normal[3*b], g.bface_normal[3*b + 1], g.bface_normal[3*b + 2]};
  }
  V3 shift(const std::size_t p) const { return {g.shift[3*p], g.shift[3*p + 1], g.shift[3*p + 2]}; }
  std::int32_t patch(const std::size_t f) const { return g.bface_patch[f - ni()]; }
  const double* value(const std::size_t face) const {
    const std::size_t p = std::size_t(patch(face));
    return f.patches[p].values.data() + k*(face - std::size_t(c.patches[p].start));
  }

  // Entries once each, the first kept: an id at one shift
  void unique(std::vector<Entry>& v) const {
    std::size_t n = 0;
    for (std::size_t i = 0; i < v.size(); ++i){
      bool seen = false;
      for (std::size_t j = 0; j < n && !seen; ++j){
        if (v[j].id != v[i].id) continue;
        bool same = true;
        for (std::size_t d = 0; d < 3; ++d) same = same && std::llround(v[j].s[d]/quantum) == std::llround(v[i].s[d]/quantum);
        seen = same;
      }
      if (!seen) v[n++] = v[i];
    }
    v.resize(n);
  }

  // The ring and its cells' face neighbours
  void ring_out(const std::vector<Entry>& in, std::vector<Entry>& out) const {
    out = in;
    for (const Entry& e : in)
      for (std::size_t l = g.link_start[e.id]; l < g.link_start[e.id + 1]; ++l)
        out.push_back({g.link[l], {e.s[0] + g.link_shift[3*l], e.s[1] + g.link_shift[3*l + 1], e.s[2] + g.link_shift[3*l + 2]}});
    unique(out);
  }

  // The ring and the faces of the same patch sharing a point, through its images, with one of its faces
  void patch_ring_out(const std::vector<Entry>& in, std::vector<Entry>& out) const {
    out = in;
    for (const Entry& e : in){
      const std::int32_t pe = patch(e.id);
      for (auto j = c.face_start[e.id]; j < c.face_start[e.id + 1]; ++j){
        const std::size_t q = std::size_t(c.face_points[std::size_t(j)]);
        const std::size_t m = std::size_t(g.master[q]);
        for (std::size_t i = g.img_start[m]; i < g.img_start[m + 1]; ++i){
          const std::size_t im = g.img[i];
          for (std::size_t b = g.pb_start[im]; b < g.pb_start[im + 1]; ++b){
            const std::size_t f2 = g.pb[b];
            if (patch(f2) != pe) continue;
            Entry n{f2, e.s};
            for (std::size_t d = 0; d < 3; ++d) n.s[d] += g.shift[3*im + d] - g.shift[3*q + d];
            out.push_back(n);
          }
        }
      }
    }
    unique(out);
  }

  // A varying patch's value at x: a fit in its tangent plane over the faces
  // touching the node, their next ring where one-sided, their mean where
  // even that is; into w.gcol, w.gcoef
  void patch_fit(Work& w, const V3& at) const {
    V3 nrm{0., 0., 0.};
    for (const Entry& e : w.group){
      const V3 n = normal(e.id);
      for (std::size_t d = 0; d < 3; ++d) nrm[d] += n[d];
    }
    nrm = unit(nrm);
    std::array<V3, 2> T;
    std::size_t nt = 2;
    const auto cross = [](const V3& a, const V3& b){
      return V3{a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
    };
    if (c.empty_axis >= 0){
      V3 e{0., 0., 0.};
      e[std::size_t(c.empty_axis)] = 1.;
      T[0] = unit(cross(nrm, e));
      nt = 1;
    }
    else {
      std::size_t a = 0;
      for (std::size_t d = 1; d < 3; ++d) if (std::abs(nrm[d]) < std::abs(nrm[a])) a = d;
      V3 e{0., 0., 0.};
      e[a] = 1.;
      T[0] = unit(cross(nrm, e));
      T[1] = cross(nrm, T[0]);
    }
    const std::vector<Entry>* ring = &w.group;
    for (int level = 1; level <= last_ring; ++level){
      if (level > 1){
        patch_ring_out(w.group, w.patch_ring);
        ring = &w.patch_ring;
      }
      const std::size_t n = ring->size();
      w.q.resize(nt*n);
      w.coef.resize(n);
      for (std::size_t i = 0; i < n; ++i){
        const Entry& e = (*ring)[i];
        const V3 fc = fcentre(e.id);
        V3 d;
        for (std::size_t j = 0; j < 3; ++j) d[j] = fc[j] + e.s[j] - at[j];
        for (std::size_t t = 0; t < nt; ++t) w.q[nt*i + t] = d[0]*T[t][0] + d[1]*T[t][1] + d[2]*T[t][2];
      }
      const Fit fit = wls(w.q.data(), n, nt, rank_tol, w.coef.data());
      if (fit.fail && level < last_ring) continue;
      for (std::size_t i = 0; i < n; ++i){
        w.gcol.push_back((*ring)[i].id);
        w.gcoef.push_back(fit.fail ? 1./double(n) : w.coef[i]);
      }
      return;
    }
  }

  // The row of master m, on a face whose condition fixes the value: per
  // patch its value at the node, the patches then averaged, or where they
  // differ and one is zero, the zero ones
  void fixed_row(Work& w, const std::size_t m, Rows& out) const {
    std::stable_sort(w.faces.begin(), w.faces.end(), [&](const Entry& a, const Entry& b){ return patch(a.id) < patch(b.id); });
    const V3 at = x(m);
    w.gcol.clear();
    w.gcoef.clear();
    w.gid.clear();
    w.gval.clear();
    std::size_t ng = 0;
    for (std::size_t i = 0; i < w.faces.size();){
      std::size_t j = i;
      while (j < w.faces.size() && patch(w.faces[j].id) == patch(w.faces[i].id)) ++j;
      const std::size_t from = w.gcol.size();
      if (uniform[std::size_t(patch(w.faces[i].id))]){
        for (std::size_t e = i; e < j; ++e){
          w.gcol.push_back(w.faces[e].id);
          w.gcoef.push_back(1./double(j - i));
        }
      }
      else {
        w.group.assign(w.faces.begin() + std::ptrdiff_t(i), w.faces.begin() + std::ptrdiff_t(j));
        patch_fit(w, at);
      }
      w.gval.resize(k*(ng + 1), 0.);
      for (std::size_t e = from; e < w.gcol.size(); ++e){
        w.gid.push_back(ng);
        const double* v = value(w.gcol[e]);
        for (std::size_t d = 0; d < k; ++d) w.gval[k*ng + d] += w.gcoef[e]*v[d];
      }
      ++ng;
      i = j;
    }
    double mdev = 0.;
    std::size_t nzero = 0;
    std::vector<char>& zero = w.zero;
    zero.assign(ng, 0);
    for (std::size_t a = 0; a < ng; ++a){
      double vmax = 0.;
      for (std::size_t d = 0; d < k; ++d){
        const double v = w.gval[k*a + d];
        mdev = std::max(mdev, std::abs(v - w.gval[d]));
        vmax = std::max(vmax, std::abs(v));
      }
      zero[a] = vmax <= tol;
      nzero += std::size_t(zero[a]);
    }
    const bool won = mdev > tol && nzero > 0;
    if (won) out.zero_wins.push_back(std::int32_t(m));
    std::size_t n = 0;
    for (std::size_t e = 0; e < w.gcol.size(); ++e){
      const std::size_t a = w.gid[e];
      const double wt = won ? (zero[a] ? 1./double(nzero) : 0.) : 1./double(ng);
      const double v = w.gcoef[e]*wt;
      if (v == 0.) continue;
      out.col.push_back(std::uint32_t(std::size_t(c.ncells) + w.gcol[e] - ni()));
      out.coef.push_back(v);
      out.mirror.push_back(0);
      ++n;
    }
    out.count.push_back(n);
    ++out.fixed;
    if (report){
      report->ring[m] = 0;
      report->cond[m] = 0.;
    }
  }

  // The row of master m by a fit over its cells, mirrored in the symmetry
  // planes through it, the next ring where the cells are too nearly coplanar
  void fitted_row(Work& w, const std::size_t m, Rows& out) const {
    const V3 at = x(m);
    // The symmetry planes through the node, one normal a patch, parallel ones once
    std::vector<V3>& planes = w.planes;
    planes.clear();
    {
      auto& acc = w.normals;
      acc.clear();
      for (std::size_t i = g.img_start[m]; i < g.img_start[m + 1]; ++i){
        const std::size_t im = g.img[i];
        for (std::size_t b = g.pb_start[im]; b < g.pb_start[im + 1]; ++b){
          const std::size_t f2 = g.pb[b];
          const std::int32_t p = patch(f2);
          if (!symmetry[std::size_t(p)]) continue;
          auto it = std::find_if(acc.begin(), acc.end(), [&](const auto& a){ return a.first == p; });
          if (it == acc.end()){
            acc.push_back({p, V3{0., 0., 0.}});
            it = acc.end() - 1;
          }
          const V3 n = normal(f2);
          for (std::size_t d = 0; d < 3; ++d) it->second[d] += n[d];
        }
      }
      std::sort(acc.begin(), acc.end(), [](const auto& a, const auto& b){ return a.first < b.first; });
      for (const auto& a : acc){
        const V3 n = unit(a.second);
        bool fresh = true;
        for (const V3& o : planes) fresh = fresh && std::abs(n[0]*o[0] + n[1]*o[1] + n[2]*o[2]) < 1 - 1e-8;
        if (fresh) planes.push_back(n);
      }
    }
    // Their reflections, and the products of every subset of them
    w.mirrors.clear();
    const std::size_t np = planes.size();
    for (std::size_t r = 1; r <= np; ++r){
      std::vector<char>& pick = w.pick;
      pick.assign(np, 0);
      std::fill(pick.begin(), pick.begin() + std::ptrdiff_t(r), 1);
      do {
        std::array<double, 9> M{1., 0., 0., 0., 1., 0., 0., 0., 1.};
        for (std::size_t i = 0; i < np; ++i){
          if (!pick[i]) continue;
          const V3& n = planes[i];
          std::array<double, 9> R, P{};
          for (std::size_t a = 0; a < 3; ++a)
            for (std::size_t b = 0; b < 3; ++b) R[3*a + b] = (a == b) - 2*n[a]*n[b];
          for (std::size_t a = 0; a < 3; ++a)
            for (std::size_t b = 0; b < 3; ++b)
              for (std::size_t e = 0; e < 3; ++e) P[3*a + b] += R[3*a + e]*M[3*e + b];
          M = P;
        }
        w.mirrors.push_back(M);
      } while (std::prev_permutation(pick.begin(), pick.end()));
    }
    // The first ring: the cells of every image
    w.ring.clear();
    for (std::size_t i = g.img_start[m]; i < g.img_start[m + 1]; ++i){
      const std::size_t im = g.img[i];
      for (std::size_t j = g.pc_start[im]; j < g.pc_start[im + 1]; ++j) w.ring.push_back({g.pc[j], shift(im)});
    }
    const std::size_t d = g.axes.size();
    const std::size_t nm = w.mirrors.size() + 1;
    for (int level = 1; level <= last_ring; ++level){
      if (level > 1){
        ring_out(w.ring, w.next);
        w.ring.swap(w.next);
      }
      const std::size_t n = w.ring.size();
      w.q.resize(d*n*nm);
      w.coef.resize(n*nm);
      if (out.gradient) w.dcoef.resize(d*n*nm);
      for (std::size_t mi = 0; mi < nm; ++mi)
        for (std::size_t i = 0; i < n; ++i){
          const Entry& e = w.ring[i];
          const V3 cc = centre(e.id);
          V3 y;
          for (std::size_t j = 0; j < 3; ++j) y[j] = cc[j] + e.s[j] - at[j];
          if (mi){
            const auto& M = w.mirrors[mi - 1];
            const V3 z = y;
            for (std::size_t a = 0; a < 3; ++a) y[a] = M[3*a]*z[0] + M[3*a + 1]*z[1] + M[3*a + 2]*z[2];
          }
          for (std::size_t a = 0; a < d; ++a) w.q[d*(mi*n + i) + a] = y[std::size_t(g.axes[a])];
        }
      const Fit fit = wls(w.q.data(), n*nm, d, rank_tol, w.coef.data(), out.gradient ? w.dcoef.data() : nullptr);
      if (fit.fail && level < last_ring) continue;
      if (out.gradient) out.dcoef.insert(out.dcoef.end(), w.dcoef.begin(), w.dcoef.end());
      std::uint32_t base = 0;
      if (nm > 1){
        base = std::uint32_t(out.table.size());
        out.table.insert(out.table.end(), w.mirrors.begin(), w.mirrors.end());
        ++out.mirrored;
      }
      for (std::size_t mi = 0; mi < nm; ++mi)
        for (std::size_t i = 0; i < n; ++i){
          out.col.push_back(std::uint32_t(w.ring[i].id));
          out.coef.push_back(w.coef[mi*n + i]);
          out.mirror.push_back(mi ? base + std::uint32_t(mi - 1) : 0);
        }
      out.count.push_back(n*nm);
      if (level > 1) ++out.second_ring;
      if (fit.fail) ++out.deficient;
      if (report){
        report->ring[m] = std::int8_t(fit.fail ? -1 : level);
        report->cond[m] = fit.cond;
      }
      return;
    }
  }

  void row(Work& w, const std::size_t m, Rows& out) const {
    w.faces.clear();
    for (std::size_t i = g.img_start[m]; i < g.img_start[m + 1]; ++i){
      const std::size_t im = g.img[i];
      for (std::size_t b = g.pb_start[im]; b < g.pb_start[im + 1]; ++b)
        if (data[std::size_t(patch(g.pb[b]))]) w.faces.push_back({g.pb[b], shift(im)});
    }
    if (w.faces.empty()) fitted_row(w, m, out);
    else fixed_row(w, m, out);
  }
};

}  // namespace

bool is_data(const openfoam_load::FieldPatch& p){
  if (!p.fixes_value) return false;
  for (const char* n : not_data) if (p.condition == n) return false;
  return true;
}

std::vector<PatchClass> classify_patches(const CaseData& c, const FieldData& u){
  if (u.patches.size() != c.patches.size() || u.ncomp != 3)
    partrac::fail(u.name, " at ", u.time, ": not a vector field of the case's ", c.patches.size(), " patches");
  double umax = 0.;
  for (const double v : u.internal) umax = std::max(umax, std::abs(v));
  for (const auto& p : u.patches) for (const double v : p.values) umax = std::max(umax, std::abs(v));
  const double tol = 1e-12*umax;
  std::vector<PatchClass> k(c.patches.size(), PatchClass::other);
  for (std::size_t i = 0; i < c.patches.size(); ++i){
    const auto& patch = c.patches[i];
    const auto& p = u.patches[i];
    if (patch.type == "cyclic"){
      k[i] = PatchClass::cyclic;
      continue;
    }
    const std::size_t size = std::size_t(patch.size), b0 = std::size_t(patch.start - c.n_internal);
    if (!is_data(p) || p.values.size() != 3*size) continue;
    bool zero = true, along = true;
    for (std::size_t f = 0; f < size; ++f){
      const double* v = &p.values[3*f];
      const double* a = &c.bface_areas[3*(b0 + f)];
      const V3 n = unit({a[0], a[1], a[2]});
      for (std::size_t d = 0; d < 3; ++d) zero = zero && std::abs(v[d]) <= tol;
      along = along && std::abs(v[0]*n[0] + v[1]*n[1] + v[2]*n[2]) <= tol;
    }
    if (zero) k[i] = PatchClass::wall;
    else if (along) k[i] = PatchClass::moving_wall;
  }
  return k;
}

const char* class_name(const PatchClass k){
  switch (k){
  case PatchClass::wall: return "no-slip wall";
  case PatchClass::moving_wall: return "moving wall";
  case PatchClass::cyclic: return "cyclic";
  default: return "other";
  }
}

std::size_t walled_cells(const CaseData& c, const std::vector<PatchClass>& k){
  std::vector<char> on_wall(c.npoints(), 0), free(std::size_t(c.ncells), 0);
  const auto points = [&](const std::size_t f){
    return std::make_pair(c.face_points.begin() + c.face_start[f], c.face_points.begin() + c.face_start[f + 1]);
  };
  for (std::size_t i = 0; i < c.patches.size(); ++i){
    if (k[i] != PatchClass::wall) continue;
    for (std::size_t f = std::size_t(c.patches[i].start); f < std::size_t(c.patches[i].start + c.patches[i].size); ++f)
      for (auto p = points(f).first; p != points(f).second; ++p) on_wall[std::size_t(*p)] = 1;
  }
  for (std::size_t f = 0; f < c.nfaces(); ++f){
    const auto pts = points(f);
    if (std::all_of(pts.first, pts.second, [&](const std::int32_t p){ return on_wall[std::size_t(p)]; })) continue;
    free[std::size_t(c.owner[f])] = 1;
    if (f < std::size_t(c.n_internal)) free[std::size_t(c.neighbour[f])] = 1;
  }
  return std::size_t(std::count(free.begin(), free.end(), 0));
}

void check_walled(const CaseData& c, const std::vector<PatchClass>& k, const int tets_per_hex, std::ostream& log){
  const std::size_t n = walled_cells(c, k);
  if (!n) return;
  if (n != std::size_t(c.ncells)){
    log << "Warning: " << n << " of " << c.ncells << " cells (" << 100.*double(n)/double(c.ncells)
        << "%) have every point on a no-slip wall, as in a gap one cell wide; "
        << (tets_per_hex == 6 ? "under split=6 their velocity is zero throughout"
                              : "under split=12 only their centre values are free, falling to zero at their points")
        << std::endl;
    return;
  }
  const char* why = "every mesh point lies on a no-slip wall, as in a case one cell thick between two walls";
  if (tets_per_hex == 6)
    partrac::fail(c.dir, ": ", why, "; under split=6 every node is fixed at zero, and so is the velocity; "
                  "split=12 keeps the cell values");
  log << "Warning: " << why << "; under split=12 only the cell centres carry the flow, each cell's velocity "
      << "falling from its centre value to zero at its points" << std::endl;
}

Geometry::Geometry(const CaseData& cd) : c(&cd){
  const std::size_t np = cd.npoints();
  const std::size_t nc = std::size_t(cd.ncells);
  const std::size_t ni = std::size_t(cd.n_internal);
  const std::size_t nb = cd.nfaces() - ni;
  master = openfoam_split::cyclic_masters(cd);
  shift.resize(3*np);
  for (std::size_t p = 0; p < np; ++p)
    for (std::size_t d = 0; d < 3; ++d) shift[3*p + d] = cd.points[3*std::size_t(master[p]) + d] - cd.points[3*p + d];
  bface_patch.assign(nb, -1);
  for (std::size_t i = 0; i < cd.patches.size(); ++i)
    std::fill(bface_patch.begin() + (cd.patches[i].start - cd.n_internal),
              bface_patch.begin() + (cd.patches[i].start + cd.patches[i].size - cd.n_internal), std::int32_t(i));
  bface_normal.resize(3*nb);
  for (std::size_t b = 0; b < nb; ++b){
    const V3 n = unit({cd.bface_areas[3*b], cd.bface_areas[3*b + 1], cd.bface_areas[3*b + 2]});
    for (std::size_t d = 0; d < 3; ++d) bface_normal[3*b + d] = n[d];
  }
  csr(np, [&](auto&& emit){
    for (std::size_t f = 0; f < cd.nfaces(); ++f)
      for (auto q = cd.face_start[f]; q < cd.face_start[f + 1]; ++q){
        const std::size_t p = std::size_t(cd.face_points[std::size_t(q)]);
        emit(p, std::size_t(cd.owner[f]));
        if (f < ni) emit(p, std::size_t(cd.neighbour[f]));
      }
  }, pc_start, pc);
  csr(np, [&](auto&& emit){
    for (std::size_t f = ni; f < cd.nfaces(); ++f)
      for (auto q = cd.face_start[f]; q < cd.face_start[f + 1]; ++q) emit(std::size_t(cd.face_points[std::size_t(q)]), f);
  }, pb_start, pb);
  csr(np, [&](auto&& emit){
    for (std::size_t p = 0; p < np; ++p) emit(std::size_t(master[p]), p);
  }, img_start, img);

  // Face neighbours: internal faces both ways, then each cyclic pair's
  link_start.assign(nc + 1, 0);
  for (std::size_t f = 0; f < ni; ++f){
    ++link_start[std::size_t(cd.owner[f]) + 1];
    ++link_start[std::size_t(cd.neighbour[f]) + 1];
  }
  for (const auto& a : cd.patches){
    if (a.type != "cyclic" || !a.owner) continue;
    const auto& b = cd.patches[std::size_t(a.neighbour)];
    for (std::int64_t i = 0; i < a.size; ++i){
      ++link_start[std::size_t(cd.owner[std::size_t(a.start + i)]) + 1];
      ++link_start[std::size_t(cd.owner[std::size_t(b.start + i)]) + 1];
    }
  }
  for (std::size_t i = 0; i < nc; ++i) link_start[i + 1] += link_start[i];
  link.resize(link_start[nc]);
  link_shift.assign(3*link.size(), 0.);
  std::vector<std::size_t> at(link_start.begin(), link_start.end() - 1);
  const auto add = [&](const std::int32_t from, const std::int32_t to, const double* s, const double sign){
    const std::size_t l = at[std::size_t(from)]++;
    link[l] = std::size_t(to);
    if (s) for (std::size_t d = 0; d < 3; ++d) link_shift[3*l + d] = sign*s[d];
  };
  for (std::size_t f = 0; f < ni; ++f) add(cd.owner[f], cd.neighbour[f], nullptr, 1.);
  for (std::size_t f = 0; f < ni; ++f) add(cd.neighbour[f], cd.owner[f], nullptr, 1.);
  for (const auto& a : cd.patches){
    if (a.type != "cyclic" || !a.owner) continue;
    const auto& b = cd.patches[std::size_t(a.neighbour)];
    const std::size_t fa = std::size_t(a.start), fb = std::size_t(b.start), n = std::size_t(a.size);
    std::vector<double> sep(3*n);
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t d = 0; d < 3; ++d)
        sep[3*i + d] = cd.bface_centres[3*(fa + i - ni) + d] - cd.bface_centres[3*(fb + i - ni) + d];
    for (std::size_t i = 0; i < n; ++i) add(cd.owner[fa + i], cd.owner[fb + i], &sep[3*i], 1.);
    for (std::size_t i = 0; i < n; ++i) add(cd.owner[fb + i], cd.owner[fa + i], &sep[3*i], -1.);
  }

  for (int d = 0; d < 3; ++d) if (d != cd.empty_axis) axes.push_back(d);
  size = openfoam_load::box_diagonal(cd);
}

LeastSquares::LeastSquares(const Geometry& g, const SplitData& s, const FieldData& f,
                           LeastSquaresReport* report, const double rank_tol){
  const CaseData& c = *g.c;
  const std::size_t np = c.npoints();
  field_ = f.name;
  ncomp_ = f.ncomp;
  ncells_ = std::size_t(c.ncells);
  n_internal_ = std::size_t(c.n_internal);
  if (f.patches.size() != c.patches.size())
    partrac::fail(f.name, " at ", f.time, ": ", f.patches.size(), " patches for the mesh's ", c.patches.size());
  if (f.ncomp != 1 && f.ncomp != 3)
    partrac::fail(f.name, " at ", f.time, " has ", f.ncomp, " components: a scalar or a vector is reconstructed");
  const std::size_t k = std::size_t(f.ncomp);

  // The patches: data or not, and a data patch uniform or varying at load
  const std::size_t npatch = c.patches.size();
  data_.assign(npatch, 0);
  Builder b{g, c, f, data_, std::vector<char>(npatch, 0), std::vector<char>(npatch, 1)};
  double vmax = 0.;
  for (std::size_t i = 0; i < npatch; ++i){
    const auto& p = f.patches[i];
    name_.push_back(c.patches[i].name);
    condition_.push_back(p.condition);
    patch_start_.push_back(std::size_t(c.patches[i].start));
    patch_size_.push_back(std::size_t(c.patches[i].size));
    data_[i] = is_data(p);
    b.symmetry[i] = c.patches[i].type == "symmetry" || c.patches[i].type == "symmetryPlane";
    if (!data_[i]) continue;
    if (!p.has_value || p.values.size() != k*patch_size_[i])
      partrac::fail(f.name, " at ", f.time, ": patch ", c.patches[i].name, " (", p.condition,
                    ") fixes the value but holds ", p.values.size(), " values for ", c.patches[i].size, " faces");
    for (const double v : p.values) vmax = std::max(vmax, std::abs(v));
  }
  b.tol = 1e-12*vmax;
  for (std::size_t i = 0; i < npatch; ++i){
    if (!data_[i]) continue;
    const auto& v = f.patches[i].values;
    for (std::size_t j = 0; j < v.size(); ++j)
      if (std::abs(v[j] - v[j % k]) > b.tol) b.uniform[i] = 0;
  }
  b.rank_tol = rank_tol;
  b.quantum = 1e-9*g.size;
  b.k = k;
  b.report = report;
  if (report){
    report->ring.assign(np, -2);
    report->cond.assign(np, 0.);
    report->zero_wins.clear();
  }

  // The masters the nodes need, ascending
  std::vector<std::int64_t> row_of;
  const std::vector<std::size_t> masters = node_masters(g, s, row_of);

  // The rows, each thread a contiguous run of masters
  std::vector<Rows> part;
#pragma omp parallel
  {
#pragma omp single
    part.resize(std::size_t(omp_get_num_threads()));
    const std::size_t t = std::size_t(omp_get_thread_num()), nt = part.size();
    const std::size_t lo = masters.size()*t/nt, hi = masters.size()*(t + 1)/nt;
    Builder::Work w;
    Rows& out = part[t];
    out.count.reserve(hi - lo);
    for (std::size_t i = lo; i < hi; ++i) b.row(w, masters[i], out);
  }

  // Joined, the mirrors into one table
  std::size_t nnz = 0;
  bool any_mirror = false;
  for (const Rows& r : part){
    nnz += r.col.size();
    any_mirror = any_mirror || r.table.size() > 1;
  }
  start_.assign(1, 0);
  start_.reserve(masters.size() + 1);
  col_.reserve(nnz);
  coef_.reserve(nnz);
  if (any_mirror){
    mirror_.reserve(nnz);
    mirror_table_.push_back({1., 0., 0., 0., 1., 0., 0., 0., 1.});
  }
  std::map<std::array<double, 9>, std::uint32_t> known;
  for (const Rows& r : part){
    for (const std::size_t n : r.count) start_.push_back(start_.back() + n);
    col_.insert(col_.end(), r.col.begin(), r.col.end());
    coef_.insert(coef_.end(), r.coef.begin(), r.coef.end());
    if (any_mirror){
      std::vector<std::uint32_t> id(r.table.size(), 0);
      for (std::size_t i = 1; i < r.table.size(); ++i){
        const auto it = known.emplace(r.table[i], std::uint32_t(mirror_table_.size()));
        if (it.second) mirror_table_.push_back(r.table[i]);
        id[i] = it.first->second;
      }
      for (const std::uint32_t m : r.mirror) mirror_.push_back(id[m]);
    }
    fixed += r.fixed;
    second_ring += r.second_ring;
    deficient += r.deficient;
    mirrored += r.mirrored;
    zero_won += r.zero_wins.size();
    if (report) report->zero_wins.insert(report->zero_wins.end(), r.zero_wins.begin(), r.zero_wins.end());
  }

  // The nodes
  kind_.assign(s.node_kind.begin(), s.node_kind.end());
  ref_.resize(s.nnodes());
  for (std::size_t n = 0; n < s.nnodes(); ++n)
    ref_[n] = std::uint32_t(s.node_point[n] >= 0 ? row_of[master_point(g, s, n)] : s.node_cell[n]);
  if (report)
    for (std::size_t p = 0; p < np; ++p){
      const std::size_t m = std::size_t(g.master[p]);
      if (row_of[m] < 0) continue;
      report->ring[p] = report->ring[m];
      report->cond[p] = report->cond[m];
    }
}

void LeastSquares::apply(const FieldData& f, const std::vector<int>& comps, std::vector<double>& out) const {
  if (f.ncomp != ncomp_)
    partrac::fail(f.name, " at ", f.time, " has ", f.ncomp, " components, the first stamp's ", field_, " ", ncomp_);
  const std::size_t k = std::size_t(f.ncomp);
  if (f.internal.size() != k*ncells_)
    partrac::fail(f.name, " at ", f.time, ": ", f.internal.size()/k, " cell values for ", ncells_, " cells");
  if (f.patches.size() != data_.size())
    partrac::fail(f.name, " at ", f.time, ": ", f.patches.size(), " patches for the mesh's ", data_.size());
  // The data: the cells, then the faces of the patches whose values are data
  std::size_t nb = 0;
  for (std::size_t i = 0; i < data_.size(); ++i) nb = std::max(nb, patch_start_[i] + patch_size_[i] - n_internal_);
  std::vector<double> b(k*(ncells_ + nb), 0.);
  std::copy(f.internal.begin(), f.internal.end(), b.begin());
  for (std::size_t i = 0; i < data_.size(); ++i){
    const auto& p = f.patches[i];
    const bool d = is_data(p);
    if (d != bool(data_[i]))
      partrac::fail(f.name, " at ", f.time, ": patch ", name_[i], " is ", p.condition, ", which ",
                    d ? "fixes" : "does not fix", " the value; at the first stamp it was ", condition_[i],
                    ", which ", d ? "did not" : "did", ", and the node values keep the first stamp's conditions");
    if (!d) continue;
    if (!p.has_value || p.values.size() != k*patch_size_[i])
      partrac::fail(f.name, " at ", f.time, ": patch ", name_[i], " (", p.condition, ") holds ",
                    p.values.size()/k, " values for ", patch_size_[i], " faces");
    std::copy(p.values.begin(), p.values.end(), b.begin() + std::ptrdiff_t(k*(ncells_ + patch_start_[i] - n_internal_)));
  }
  const std::size_t nc = comps.size();
  const std::size_t nn = kind_.size();
  out.assign(nn*nc, 0.);
  const bool mirrors = k == 3 && !mirror_.empty();
#pragma omp parallel for schedule(static)
  for (std::size_t n = 0; n < nn; ++n){
    double v[3] = {0., 0., 0.};
    if (kind_[n]){
      for (std::size_t d = 0; d < k; ++d) v[d] = f.internal[k*ref_[n] + d];
    }
    else {
      const std::size_t r = ref_[n];
      for (std::size_t e = start_[r]; e < start_[r + 1]; ++e){
        const double* x = b.data() + k*col_[e];
        const double w = coef_[e];
        if (mirrors && mirror_[e]){
          const auto& M = mirror_table_[mirror_[e]];
          for (std::size_t i = 0; i < 3; ++i) v[i] += w*(M[3*i]*x[0] + M[3*i + 1]*x[1] + M[3*i + 2]*x[2]);
        }
        else {
          for (std::size_t d = 0; d < k; ++d) v[d] += w*x[d];
        }
      }
    }
    for (std::size_t j = 0; j < nc; ++j) out[n*nc + j] = v[comps[j]];
  }
}

Gradient::Gradient(const Geometry& g, const SplitData& s, const double rank_tol){
  const CaseData& c = *g.c;
  const std::size_t npatch = c.patches.size();
  d_ = g.axes.size();
  ncells_ = std::size_t(c.ncells);
  const FieldData none;
  const std::vector<char> no_data(npatch, 0);
  Builder b{g, c, none, no_data, std::vector<char>(npatch, 0), std::vector<char>(npatch, 1)};
  for (std::size_t i = 0; i < npatch; ++i)
    b.symmetry[i] = c.patches[i].type == "symmetry" || c.patches[i].type == "symmetryPlane";
  b.rank_tol = rank_tol;
  b.quantum = 1e-9*g.size;

  // The masters the nodes need, ascending
  std::vector<std::int64_t> row_of;
  const std::vector<std::size_t> masters = node_masters(g, s, row_of);

  // The rows, each thread a contiguous run of masters
  std::vector<Rows> part;
#pragma omp parallel
  {
#pragma omp single
    part.resize(std::size_t(omp_get_num_threads()));
    const std::size_t t = std::size_t(omp_get_thread_num()), nt = part.size();
    const std::size_t lo = masters.size()*t/nt, hi = masters.size()*(t + 1)/nt;
    Builder::Work w;
    Rows& out = part[t];
    out.gradient = true;
    out.count.reserve(hi - lo);
    for (std::size_t i = lo; i < hi; ++i) b.fitted_row(w, masters[i], out);
  }
  start_.assign(1, 0);
  for (const Rows& r : part){
    for (const std::size_t n : r.count) start_.push_back(start_.back() + n);
    col_.insert(col_.end(), r.col.begin(), r.col.end());
    coef_.insert(coef_.end(), r.dcoef.begin(), r.dcoef.end());
    second_ring += r.second_ring;
    deficient += r.deficient;
    mirrored += r.mirrored;
  }
  ref_.resize(s.nnodes());
  for (std::size_t n = 0; n < s.nnodes(); ++n)
    ref_[n] = std::int32_t(s.node_point[n] >= 0 ? row_of[master_point(g, s, n)] : -1);
}

void Gradient::apply(const FieldData& f, std::vector<double>& out) const {
  if (f.ncomp != 1)
    partrac::fail(f.name, " at ", f.time, " has ", f.ncomp, " components; a gradient is of a scalar");
  if (f.internal.size() != ncells_)
    partrac::fail(f.name, " at ", f.time, ": ", f.internal.size(), " cell values for ", ncells_, " cells");
  const std::size_t nn = ref_.size(), d = d_;
  out.assign(nn*d, 0.);
#pragma omp parallel for schedule(static)
  for (std::size_t n = 0; n < nn; ++n){
    if (ref_[n] < 0) continue;
    const std::size_t r = std::size_t(ref_[n]);
    double G[3] = {0., 0., 0.};
    for (std::size_t e = start_[r]; e < start_[r + 1]; ++e){
      const double v = f.internal[col_[e]];
      for (std::size_t a = 0; a < d; ++a) G[a] += coef_[d*e + a]*v;
    }
    for (std::size_t a = 0; a < d; ++a) out[n*d + a] = G[a];
  }
}

std::size_t Gradient::bytes() const {
  return 8*(start_.capacity() + coef_.capacity()) + 4*(ref_.capacity() + col_.capacity());
}

std::size_t LeastSquares::bytes() const {
  return kind_.capacity() + 4*(ref_.capacity() + col_.capacity() + mirror_.capacity())
       + 8*(start_.capacity() + coef_.capacity() + patch_start_.capacity() + patch_size_.capacity())
       + 72*mirror_table_.capacity() + data_.capacity();
}

std::size_t InverseDistance::bytes() const {
  return kind_.capacity() + 8*(ref_.capacity() + start_.capacity() + src_.capacity() + w_.capacity()
                               + den_.capacity() + patch_start_.capacity())
       + 4*(constraint_.capacity() + bface_patch_.capacity() + bface_owner_.capacity())
       + 24*bface_normal_.capacity() + 72*t_.capacity() + symmetry_.capacity();
}

}  // namespace openfoam_nodes
