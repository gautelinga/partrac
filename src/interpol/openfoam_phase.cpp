#include "openfoam_phase.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace openfoam_phase {

double above_half(const double* v, const int nv){
  double s[4];
  int pos[4], neg[4], np = 0, nn = 0;
  for (int i = 0; i < nv; ++i){
    s[i] = v[i] - 0.5;
    if (s[i] > 0.) pos[np++] = i;
    else neg[nn++] = i;
  }
  if (!np) return 0.;
  if (!nn) return 1.;
  // A node alone on its side: the product of its edges' cut fractions
  const auto corner = [&](const int a, const int* o, const int no){
    double f = 1.;
    for (int j = 0; j < no; ++j) f *= s[a]/(s[a] - s[o[j]]);
    return f;
  };
  if (np == 1) return corner(pos[0], neg, nn);
  if (nn == 1) return 1. - corner(neg[0], pos, np);
  // Two and two: a prism, as three tets
  const double a = s[pos[0]], b = s[pos[1]], c = s[neg[0]], d = s[neg[1]];
  return a*a/((a - c)*(a - d)) - a*b*c/((a - c)*(a - d)*(b - c)) - b*b*d/((a - d)*(b - c)*(b - d));
}

Cells::Cells(const std::size_t ncells, const std::vector<std::int32_t>& cell_of, const std::vector<std::int32_t>& centre_)
  : centre(centre_)
{
  start.assign(ncells + 1, 0);
  for (const std::int32_t c : cell_of) ++start[std::size_t(c) + 1];
  for (std::size_t c = 0; c < ncells; ++c) start[c + 1] += start[c];
  simplex.resize(cell_of.size());
  std::vector<std::uint32_t> at(start.begin(), start.end() - 1);
  for (std::size_t i = 0; i < cell_of.size(); ++i) simplex[at[std::size_t(cell_of[i])]++] = std::uint32_t(i);
}

namespace {

// A simplex's measure: volume, or area in 2D
double measure(const std::uint32_t* t, const std::vector<double>& x, const int nv){
  const int d = nv - 1;
  const double* p0 = &x[std::size_t(d)*t[0]];
  double e[3][3] = {};
  for (int i = 1; i < nv; ++i)
    for (int j = 0; j < d; ++j) e[i - 1][j] = x[std::size_t(d)*t[i] + std::size_t(j)] - p0[j];
  if (nv == 3) return std::abs(e[0][0]*e[1][1] - e[0][1]*e[1][0])/2;
  return std::abs(e[0][0]*(e[1][1]*e[2][2] - e[1][2]*e[2][1]) - e[0][1]*(e[1][0]*e[2][2] - e[1][2]*e[2][0])
                  + e[0][2]*(e[1][0]*e[2][1] - e[1][1]*e[2][0]))/6;
}

}  // namespace

Report conserve(const Cells& cells, const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                const int nv, const std::vector<double>& alpha, std::vector<double>& values){
  Report r;
  const std::size_t nc = cells.ncells(), k = std::size_t(nv);
#pragma omp parallel
  {
    Report t;
    std::vector<double> vol, val;
    std::vector<int> at;
#pragma omp for schedule(dynamic, 256) nowait
    for (std::size_t c = 0; c < nc; ++c){
      const std::size_t s0 = cells.start[c], n = cells.start[c + 1] - s0;
      const std::int64_t centre = cells.centre[c];
      vol.resize(n);
      val.resize(k*n);
      at.assign(n, -1);
      double vc = 0.;
      for (std::size_t i = 0; i < n; ++i){
        const std::uint32_t* q = &topo[k*cells.simplex[s0 + i]];
        vol[i] = measure(q, coords, nv);
        vc += vol[i];
        for (std::size_t j = 0; j < k; ++j){
          val[k*i + j] = values[q[j]];
          if (std::int64_t(q[j]) == centre) at[i] = int(j);
        }
      }
      // The cell's enclosed volume with its centre at x
      const auto enclosed = [&](const double x){
        double e = 0.;
        for (std::size_t i = 0; i < n; ++i){
          if (at[i] >= 0) val[k*i + std::size_t(at[i])] = x;
          e += vol[i]*above_half(&val[k*i], nv);
        }
        return e;
      };
      const double want = vc*std::clamp(alpha[c], 0., 1.);
      const double tol = 1e-14*vc;
      const double x0 = centre >= 0 ? values[std::size_t(centre)] : 0.;
      const double e0 = enclosed(x0) - want;
      double x = x0, e = e0;
      if (centre < 0){
        if (std::abs(e0) > tol) ++t.unfanned;
      }
      else if ((alpha[c] <= 0. || alpha[c] >= 1.) && std::abs(e0) > tol) ++t.unreachable;
      else if (alpha[c] < purity || alpha[c] > 1. - purity) ++t.pure;
      else if (std::abs(e0) <= tol) ++t.corrected;
      else {
        // Within [0, 1], else its end nearer the target
        double a = 0., fa = enclosed(0.) - want, b = 1., fb = enclosed(1.) - want;
        if (fa >= 0. || fb <= 0.){
          x = fa >= 0. ? 0. : 1.;
          e = fa >= 0. ? fa : fb;
          if (std::abs(e) > tol) ++t.clipped;
          else ++t.corrected;
        }
        else {
          // Illinois between a and b
          int side = 0;
          for (int it = 0; it < 200; ++it){
            double m = (a*fb - b*fa)/(fb - fa);
            if (!(m > a && m < b)) m = 0.5*(a + b);
            const double fm = enclosed(m) - want;
            x = m;
            e = fm;
            if (std::abs(fm) <= tol) break;
            if (fm > 0.){
              b = m;
              fb = fm;
              if (side == -1) fa *= 0.5;
              side = -1;
            }
            else {
              a = m;
              fa = fm;
              if (side == 1) fb *= 0.5;
              side = 1;
            }
            if (b - a <= 4*std::numeric_limits<double>::epsilon()) break;
          }
          ++t.corrected;   // or stopped on the bracket's width, the logged error saying by how much
        }
        values[std::size_t(centre)] = x;
      }
      t.target += want;
      t.before += e0 + want;
      t.after += e + want;
      t.worst_before = std::max(t.worst_before, std::abs(e0)/vc);
      if (std::abs(e)/vc > t.worst_after || t.worst_cell < 0){
        t.worst_after = std::abs(e)/vc;
        t.worst_cell = std::int64_t(c);
      }
    }
#pragma omp critical
    {
      r.target += t.target;
      r.before += t.before;
      r.after += t.after;
      r.worst_before = std::max(r.worst_before, t.worst_before);
      if (t.worst_cell >= 0 && (t.worst_after > r.worst_after || r.worst_cell < 0
                                || (t.worst_after == r.worst_after && t.worst_cell < r.worst_cell))){
        r.worst_after = t.worst_after;
        r.worst_cell = t.worst_cell;
      }
      r.pure += t.pure;
      r.corrected += t.corrected;
      r.clipped += t.clipped;
      r.unreachable += t.unreachable;
      r.unfanned += t.unfanned;
    }
  }
  return r;
}

void centre_means(const Cells& cells, const std::vector<std::uint32_t>& topo, const int nv, const int d,
                  std::vector<double>& grad){
  const std::size_t k = std::size_t(nv), dd = std::size_t(d);
#pragma omp parallel
  {
    std::vector<std::uint32_t> nodes;
#pragma omp for schedule(static)
    for (std::size_t c = 0; c < cells.ncells(); ++c){
      const std::int64_t centre = cells.centre[c];
      if (centre < 0) continue;
      nodes.clear();
      for (std::size_t i = cells.start[c]; i < cells.start[c + 1]; ++i)
        for (std::size_t j = 0; j < k; ++j){
          const std::uint32_t q = topo[k*cells.simplex[i] + j];
          if (std::int64_t(q) != centre) nodes.push_back(q);
        }
      std::sort(nodes.begin(), nodes.end());
      nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
      double* g = &grad[dd*std::size_t(centre)];
      for (std::size_t j = 0; j < dd; ++j){
        double sum = 0.;
        for (const std::uint32_t q : nodes) sum += grad[dd*q + j];
        g[j] = nodes.empty() ? 0. : sum/double(nodes.size());
      }
    }
  }
}

}  // namespace openfoam_phase
