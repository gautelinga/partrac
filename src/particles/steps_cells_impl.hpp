#ifndef __STEPS_CELLS_IMPL_HPP
#define __STEPS_CELLS_IMPL_HPP

// The loop of RK4cells; included only by its per-interpolator step sources

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>
#include "CellsIntegrator.hpp"
#include "held_eval.hpp"
#include "regions.hpp"
#include "steps_impl.hpp"

namespace cells {

constexpr double band = 1e-9;       // levels down to -band are inside
constexpr double t_eps = 1e-12;     // of dt: time left that ends a step, roots on the entered facet ignored
constexpr double land_max = 0.1;    // a landing moving more than this fraction of the step rejects it
constexpr int cap_free = 10;        // steps of a major step before the cap is checked
constexpr double proj_max = 1e-6;   // a projection onto walls moving a level more than this fails
constexpr int chunk = 16;           // particles a thread takes at once
constexpr double inf = std::numeric_limits<double>::infinity();

// A particle's major step: the element unnormalised, the frame advanced
template<TransportElement E, typename Interp>
struct State {
  typename Interp::region_type R;
  Vector3d x;
  Vector3d n;
  Matrix3d F;
  double t;              // time reached
  double r;              // time left
  double hc = inf;       // the step after a rejection
  double umax2 = 0.;     // the fastest first stage, squared
  int ent = -1;          // the facet entered
  int count = 0;         // steps, crossings and rejections so far
  bool nox = false;      // the next step ignores its crossing
  // The first stage at (x, t): values, levels and their rates
  Vector3d u1;
  Matrix3d J1;
  typename Interp::levels_type lev, d1, d2;
  // The cap, for the region and fastest first stage it was set from; set at the first check
  double cap;
  double cap_u2;
  int cap_id;
};

// The held stages as the element needs them: a point's without the gradients
template<TransportElement E, typename Held>
struct Stages {
  Held& ev;
  bool operator()(const Vector3d& x, const double t){
    if constexpr (E == TransportElement::Point) return ev.velocity(x, t);
    else return ev(x, t);
  }
  Vector3d u(){ return ev.u(); }
  Matrix3d J(){ return ev.J(); }
  bool end(const Vector3d& x, const double t){ return ev.end(x, t); }
};

// The end evaluated as well, for the landing
template<TransportElement E, typename Held>
struct ToFacet : Stages<E, Held> {
  bool end(const Vector3d& x, const double t){ return this->ev.at(x, t); }
};

// The time a level l + b tau + a tau^2 falls through zero, a and b not both rising; inf if never
inline double fall_time(const double l, const double b, const double a){
  if (std::abs(a) < 1e-14*(std::abs(b) + 1e-300))
    return b < 0. ? -l/b : inf;
  const double disc = b*b - 4.*a*l;
  if (disc < 0.)
    return inf;
  const double sq = std::sqrt(disc);
  const double q = -0.5*(b + (b >= 0. ? sq : -sq));
  double tau = inf;
  for (const double rr : {q/a, l/q})
    if (std::isfinite(rr) && rr > 0. && b + 2.*a*rr <= 0. && rr < tau) tau = rr;
  return tau;
}

// The first time in (0, horizon] a level falls through zero along l + d1 tau
// + d2 tau^2 / 2, l clamped at zero, and its facet k; 0 for a level within the
// band falling now, a wall's aside; inf if none. The entered facet's roots below eps ignored
template<int n_levels, typename Interp, typename Levels>
inline double first_crossing(const Interp& intp, const typename Interp::region_type& R, const Levels& lev,
                             const Levels& d1, const Levels& d2,
                             const double horizon, const int ent, const double eps, int& k){
  double best = inf;
  for (int m = 0; m < n_levels; ++m){
    const double l = std::max(lev[m], 0.), b = d1[m], a = 0.5*d2[m];
    if (b < 0. && l <= band){
      // Resting on a wall: the held step
      if (intp.is_wall(R, m)) continue;
      k = m;
      return 0.;
    }
    // Rising; above zero at the horizon and concave, still falling there or never reaching zero
    if (b >= 0. && a >= 0.) continue;
    if (l + horizon*(b + a*horizon) > 0. && (a <= 0. || b + 2.*a*horizon <= 0. || b*b < 4.*a*l)) continue;
    const double tau = fall_time(l, b, a);
    if (tau < best && !(m == ent && tau < eps)){
      best = tau;
      k = m;
    }
  }
  return best <= horizon ? best : inf;
}

// The first stage at the state's start, its levels and their rates; false if
// the start is below the band
template<TransportElement E, typename Interp>
inline bool first_stage(Interp& intp, HeldEval<Interp>& ev, PointValues& pv, State<E, Interp>& s){
  if (!ev.at(s.x, s.t))
    return false;
  s.u1 = ev.u();
  s.J1 = ev.J();
  s.lev = ev.levels();
  s.umax2 = std::max(s.umax2, s.u1.squaredNorm());
  intp.level_rates(s.R, s.u1, s.d1);
  intp.level_rates(s.R, Vector3d(s.J1*s.u1 + pv.get_a()), s.d2);
  return true;
}

// Rare: the start below the band, so the region the position is in
template<TransportElement E, typename Interp>
__attribute__((noinline))
bool relocate(Interp& intp, HeldEval<Interp>& ev, PointValues& pv, State<E, Interp>& s, const int id, CellsCounts& c){
  typename Interp::region_type R;
  typename Interp::levels_type lev;
  if (!intp.region_of(id, s.x, band, R, lev))
    return false;
  if (R.id != id && id >= 0)
    ++c.relocations;
  s.R = R;
  s.ent = -1;
  ev.bind(R);
  return first_stage(intp, ev, pv, s);
}

// Rare: past the cap on steps a major step, from the fastest first stage and the region's size
template<TransportElement E, typename Interp>
__attribute__((noinline))
bool over_cap(Interp& intp, State<E, Interp>& s, const double dt){
  if (s.count == cap_free || s.R.id != s.cap_id || s.umax2 != s.cap_u2){
    s.cap = 5.*dt*std::sqrt(s.umax2)/intp.region_size(s.R.id) + cap_free;
    s.cap_id = s.R.id;
    s.cap_u2 = s.umax2;
  }
  return s.count >= s.cap;
}

// Into the region across facet k; false at a wall
template<TransportElement E, typename Interp>
inline bool enter(Interp& intp, HeldEval<Interp>& ev, State<E, Interp>& s, const int k, const Vector3d& x,
                  CellsCounts& c){
  typename Interp::region_type next;
  const int entered = intp.across(s.R, k, x, next);
  if (entered < 0)
    return false;
  s.R = next;
  s.ent = entered;
  s.hc = inf;
  ev.bind(next);
  ++c.crossings;
  return true;
}

// Rare: a level within the band falling, across now
template<TransportElement E, typename Interp>
__attribute__((noinline))
bool cross_now(Interp& intp, HeldEval<Interp>& ev, PointValues& pv, State<E, Interp>& s, const int k, CellsCounts& c){
  ++c.immediate;
  if (!enter(intp, ev, s, k, s.x, c))
    return false;
  // The entered facet untested: from a wider cell its level is below the band
  ev.skip = s.ent;
  const bool ok = first_stage(intp, ev, pv, s);
  ev.skip = -1;
  return ok || relocate(intp, ev, pv, s, s.R.id, c);
}

// A rejection: half the step next, the first stage kept
template<TransportElement E, typename Interp>
inline void reject(State<E, Interp>& s, const double h, CellsCounts& c){
  ++c.rejects;
  s.hc = 0.5*h;
  s.nox = false;
}

// The held step to the predicted crossing of facet k at h, its landing on the
// facet by a Newton correction to second order from the end's values, then
// the region across; false at a wall or a failed region. r_end: time left that ends the major step
template<TransportElement E, typename Interp>
__attribute__((noinline))
bool step_to_facet(Interp& intp, HeldEval<Interp>& ev, PointValues& pv, State<E, Interp>& s, const int k, const double h,
                   const double r_end, CellsCounts& c){
  Vector3d dx;
  [[maybe_unused]] Vector3d el;
  [[maybe_unused]] Matrix3d dF;
  ToFacet<E, HeldEval<Interp>> st{{ev}};
  ev.skip = k;
  const bool ok = rk_stages_from<RK4Tableau, E, true>(st, s.u1, s.J1, s.x, s.n, s.F, s.t, h, dx, el, dF);
  ev.skip = -1;
  ++c.steps;
  if (!ok){
    reject(s, h, c);
    return true;
  }
  // Facet k's level to zero along the path from the end
  const Vector3d ue = ev.u();
  const Matrix3d Je = ev.J();
  const Vector3d w = Je*ue + pv.get_a();
  const double lk = ev.levels()[std::size_t(k)];
  typename Interp::levels_type g1, g2;
  intp.level_rates(s.R, ue, g1);
  intp.level_rates(s.R, w, g2);
  const double rate = g1[std::size_t(k)];
  if (!(rate < 0.)){
    reject(s, h, c);
    return true;
  }
  const double disc = rate*rate - 2.*g2[std::size_t(k)]*lk;
  const double delta = disc >= 0. ? -2.*lk/(rate - std::sqrt(disc)) : -lk/rate;
  if (std::abs(delta) > land_max*h){
    reject(s, h, c);
    return true;
  }
  // Past the major step's end: the facet comes after the meeting
  if (h + delta > s.r*(1. + 1e-12)){
    ++c.rejects;
    s.nox = true;
    return true;
  }
  const Vector3d x = s.x + dx + ue*delta + w*(0.5*delta*delta);
  if (!enter(intp, ev, s, k, x, c))
    return false;
  s.x = x;
  if constexpr (E == TransportElement::Vector){
    const Vector3d Jn = Je*el;
    s.n = el + Jn*delta + Je*Jn*(0.5*delta*delta);
  }
  if constexpr (E == TransportElement::Tensor){
    const Matrix3d Fe = s.F + dF;
    const Matrix3d JF = Je*Fe;
    s.F = Fe + JF*delta + Je*JF*(0.5*delta*delta);
  }
  s.t += h + delta;
  s.r -= h + delta;
  s.nox = false;
  // The end's levels in the region entered, else in the region holding it
  if (s.r <= r_end)
    return ev.end(s.x, s.t) || relocate(intp, ev, pv, s, s.R.id, c);
  return first_stage(intp, ev, pv, s) || relocate(intp, ev, pv, s, s.R.id, c);
}

// The major step from the first stage at the start; false: the fallback
template<TransportElement E, typename Interp>
inline bool major_step(Interp& intp, HeldEval<Interp>& ev, PointValues& pv, State<E, Interp>& s, const double dt,
                       CellsCounts& c){
  constexpr int n_levels = Interp::n_levels;
  while (s.r > t_eps*dt){
    if (s.count >= cap_free && over_cap(intp, s, dt))
      return false;
    ++s.count;
    const double horizon = std::min(s.r, s.hc);
    int k = -1;
    const double tau = first_crossing<n_levels>(intp, s.R, s.lev, s.d1, s.d2, horizon, s.ent, t_eps*dt, k);
    if (tau == 0.){
      if (!cross_now(intp, ev, pv, s, k, c))
        return false;
      continue;
    }
    if (tau < inf && !s.nox){
      if (!step_to_facet(intp, ev, pv, s, k, tau, t_eps*dt, c))
        return false;
      continue;
    }
    // No crossing: RK4 held in the cell
    Vector3d dx;
    [[maybe_unused]] Vector3d el;
    [[maybe_unused]] Matrix3d dF;
    Stages<E, HeldEval<Interp>> st{ev};
    const bool ok = rk_stages_from<RK4Tableau, E, true>(st, s.u1, s.J1, s.x, s.n, s.F, s.t, horizon, dx, el, dF);
    ++c.steps;
    if (!ok){
      reject(s, horizon, c);
      continue;
    }
    s.x += dx;
    if constexpr (E == TransportElement::Vector) s.n = el;
    if constexpr (E == TransportElement::Tensor) s.F += dF;
    s.t += horizon;
    s.r -= horizon;
    s.hc = inf;
    s.nox = false;
    if (s.r > t_eps*dt && !first_stage(intp, ev, pv, s) && !relocate(intp, ev, pv, s, s.R.id, c))
      return false;
  }
  return true;
}

// Outside: below zero, or on a face a half-open region leaves to the next
template<typename Interp>
inline bool below(const double l){ return l < 0. || (Interp::half_open && l == 0.); }

// Onto the walls of mask at once and inside; false if they do not meet, the
// point leaves the region or a level moves more than proj_max
template<TransportElement E, typename Interp>
bool onto_walls(Interp& intp, State<E, Interp>& s, const int mask, typename Interp::levels_type& lev){
  constexpr int n_levels = Interp::n_levels;
  using Mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, n_levels, n_levels>;
  using Vec = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, n_levels, 1>;
  std::array<int, n_levels> idx;
  std::array<Vector3d, n_levels> g;
  int w = 0;
  for (int m = 0; m < n_levels; ++m)
    if ((mask >> m) & 1){
      idx[w] = m;
      g[w] = intp.level_grad(s.R, m);
      ++w;
    }
  Mat M(w, w);
  for (int i = 0; i < w; ++i)
    for (int j = 0; j < w; ++j) M(i, j) = g[i].dot(g[j]);
  const Eigen::FullPivLU<Mat> lu(M);
  if (!lu.isInvertible())
    return false;
  const typename Interp::levels_type lev0 = lev;
  for (double eps = 0.;; eps = eps > 0. ? 4.*eps : 1e-16){
    Vec rhs(w);
    for (int i = 0; i < w; ++i) rhs[i] = eps - lev[idx[i]];
    const Vec lam = lu.solve(rhs);
    for (int i = 0; i < w; ++i) s.x += g[i]*lam[i];
    intp.region_point(s.R, s.x, lev);
    bool in = true;
    for (int i = 0; i < w; ++i) in = in && !below<Interp>(lev[idx[i]]);
    if (in) break;
    if (eps > band) return false;
  }
  for (int m = 0; m < n_levels; ++m)
    if (lev[m] < -band || std::abs(lev[m] - lev0[m]) > proj_max) return false;
  return true;
}

// Rare: an end below a wall, onto it and inside until contains takes it; then
// onto every wall met at once while one is below (walls at an acute angle).
// False if that fails: the fallback
template<TransportElement E, typename Interp>
__attribute__((noinline))
bool project_end(Interp& intp, State<E, Interp>& s){
  constexpr int n_levels = Interp::n_levels;
  typename Interp::levels_type lev;
  intp.region_point(s.R, s.x, lev);
  int met = 0;
  for (int m = 0; m < n_levels; ++m){
    if (!below<Interp>(lev[m]) || !intp.is_wall(s.R, m)) continue;
    met |= 1 << m;
    const Vector3d g = intp.level_grad(s.R, m);
    const double g2 = g.squaredNorm();
    s.x -= g*(lev[m]/g2);
    intp.region_point(s.R, s.x, lev);
    for (double eps = 1e-16; below<Interp>(lev[m]); eps *= 4.){
      s.x += g*((eps - lev[m])/g2);
      intp.region_point(s.R, s.x, lev);
    }
  }
  for (int pass = 0; pass < n_levels; ++pass){
    int now = 0;
    for (int m = 0; m < n_levels; ++m)
      if (below<Interp>(lev[m]) && intp.is_wall(s.R, m)) now |= 1 << m;
    if (!now)
      return true;
    met |= now;
    if (!onto_walls(intp, s, met, lev))
      return false;
  }
  return false;
}

}  // namespace cells

template<TransportElement E, typename Interp>
PARTRAC_HOT_LOOP
std::vector<Uint> RK4CellsIntegrator::step(Interp& intp, ParticleSet& ps, const double t, const double dt){
  static_assert(has_regions<Interp>::value, "RK4cells needs an interpolator with regions");
  std::vector<Uint> outside_nodes;
  CellsCounts counts;
  #pragma omp parallel
  {
    std::vector<Uint> outside_nodes_loc;
    Uint n_accepted_loc = 0;
    Uint n_declined_loc = 0;
    CellsCounts c;
    // Kept per thread: a frequency loader's held nodes own memory
    PointValues pv(intp.get_U0());
    HeldEval<Interp> ev(intp, pv);
    ev.band = cells::band;

    #pragma omp for schedule(dynamic, cells::chunk)
    for (Uint i = 0; i < ps.N(); ++i){
      cells::State<E, Interp> s;
      s.x = ps.x(i);
      if constexpr (E == TransportElement::Vector) s.n = ps.rhohat(i);
      if constexpr (E == TransportElement::Tensor) s.F = ps.frame(i);
      s.t = t;
      s.r = dt;
      const int id = ps.get_cell_id(i);
      bool ok = false;
      if (id >= 0){
        s.R = intp.region_in(id, s.x);
        ev.bind(s.R);
        ok = cells::first_stage(intp, ev, pv, s);
      }
      // A start below the band: located
      ok = ok || cells::relocate(intp, ev, pv, s, id, c);
      // An end below a wall: onto it
      if (ok && cells::major_step(intp, ev, pv, s, dt, c)
          && (!cells::below<Interp>(min_level<Interp::n_levels>(ev.levels())) || cells::project_end(intp, s))){
        ps.set_x(i, s.x);
        ps.set_t_loc(i, ps.t_loc(i) + dt);
        if constexpr (E == TransportElement::Vector){
          const double len = s.n.norm();
          ps.set_rhohat(i, s.n/len);
          ps.set_w(i, ps.w(i) + log(len));
        }
        if constexpr (E == TransportElement::Tensor)
          ps.advance_frame(i, s.F);
        ps.set_cell_id(i, s.R.id);
        ++n_accepted_loc;
        continue;
      }
      // At a wall, a failed region or projection, or past the cap: RK4 substeps from the start
      ++c.fallbacks;
      if (rk4_substeps<E>(intp, ps, i, t, dt))
        ++n_accepted_loc;
      else {
        outside_nodes_loc.push_back(i);
        ++n_declined_loc;
      }
    }
    #pragma omp critical
    {
      outside_nodes.insert(outside_nodes.end(), outside_nodes_loc.begin(), outside_nodes_loc.end());
      n_accepted += n_accepted_loc;
      n_declined += n_declined_loc;
      counts += c;
    }
  }
  counts_ += counts;
  particle_steps_ += ps.N();
  // Sort: threads merge out of order
  std::sort(outside_nodes.begin(), outside_nodes.end());
  return outside_nodes;
}

#endif
