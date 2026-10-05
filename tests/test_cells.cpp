// RK4cells (CellsIntegrator.hpp, steps_cells_impl.hpp) against RK4 on cases
// dolfin writes here, through the stepping layer the apps call.
//
// A major step of RK4cells is RK4 held in the cell the particle is in, cut
// where the path is predicted to meet a facet. A step that meets none is then
// RK4's step with every stage evaluated from the same cell's nodes, and the
// held evaluation is the loader's own bit for bit: the two schemes must give
// the same bits, in periodic boxes (positions unwrapped by whole periods),
// near walls under the near-wall P2 rule, between stamps and for the
// frequency loaders. A particle started on a facet, in the cell it is leaving,
// crosses before it steps (an immediate crossing); one stored in a cell that
// no longer holds it is relocated first. Both then step as RK4 does, the
// field being continuous across the facet, to round-off. A wall facet never
// makes a particle cross before it steps: one resting within the band of a
// wall steps held and ends inside; one leaving through it falls back. The
// same on the divergence-free split, whose regions are sub-cells (planes
// inside a cell as well as its facets), and on the felbm lattice, whose
// regions are floor cells and, next to solids, sub-cubes, with walls at the
// sub-cubes' mid-planes before a solid node. And the rare paths: a landing
// that ends the major step, or would land past it; a crossing into a narrower
// cell from within the band; the cap on steps a major step; an end below two
// walls at an acute angle, the joint projection onto walls and its failures;
// the lattice's cells through a checkpoint. Fourth order where steps land on
// predicted facets: split triangles, and line elements and frames on XDMF triangles.
#ifdef USE_DOLFIN
#include <catch2/catch.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <vector>

#include "ParticleSet.hpp"
#include "SimplexInterpol.hpp"
#include "SplitInterpol.hpp"
#include "StructuredInterpol.hpp"
#include "interpolators.hpp"
#include "XDMFInterpol.hpp"
#include "region_cases.hpp"
#include "stepping.hpp"
#include "steps_cells_impl.hpp"

namespace {

// Points well inside regions of the loader (every level above lo), some a few
// periods off along the periodic axes
template<typename I>
std::vector<Vector3d> inner_points(I& intp, const int D, const std::vector<bool>& per, const double t,
                                   const std::size_t n, const double lo){
  std::vector<Vector3d> out;
  for (const Vector3d& x : box_points(intp, D, per, 20*n)){
    CellPos pos;
    if (!intp.locate(x, t, pos)) continue;
    typename I::region_type R;
    typename I::levels_type lev;
    if (!intp.region_of(pos.id, x, 0., R, lev)) continue;
    bool inner = true;
    for (int k = 0; k < I::n_levels; ++k) inner = inner && lev[k] > lo;
    if (!inner) continue;
    out.push_back(x);
    if (out.size() == n) break;
  }
  return out;
}

// A particle set on intp carrying E from the points, with random unit line
// elements and random frames; cells given, or unknown
template<TransportElement E>
std::unique_ptr<ParticleSet> particles(std::shared_ptr<Interpol> intp, const std::vector<Vector3d>& pts,
                                       const std::vector<int>& cells = {}){
  auto ps = std::make_unique<ParticleSet>(intp, Uint(pts.size()));
  ps->carry(E);
  ps->add(pts, 0);
  std::mt19937 rng(23);
  std::normal_distribution<double> g(0., 1.);
  for (Uint i = 0; i < ps->N(); ++i){
    if (!cells.empty()) ps->set_cell_id(i, cells[i]);
    if constexpr (E == TransportElement::Vector) ps->set_rhohat(i, Vector3d(g(rng), g(rng), g(rng)).normalized());
    if constexpr (E == TransportElement::Tensor){
      Matrix3d F;
      for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) F(a, b) = (a == b) + 0.3*g(rng);
      ps->set_F(i, F);
    }
  }
  return ps;
}

// The largest difference of position and carried element between two sets, and whether every bit agrees
template<TransportElement E>
double largest_difference(const ParticleSet& a, const ParticleSet& b, bool& same){
  double d = 0.;
  same = true;
  for (Uint i = 0; i < a.N(); ++i){
    d = std::max(d, (a.x(i) - b.x(i)).norm());
    // RK4 keeps no cell on a lattice
    same = same && same_bits(a.x(i), b.x(i)) && (a.get_cell_id(i) == b.get_cell_id(i) || a.get_cell_id(i) < 0);
    if constexpr (E == TransportElement::Vector){
      d = std::max(d, (a.rhohat(i) - b.rhohat(i)).norm());
      same = same && same_bits(a.rhohat(i), b.rhohat(i)) && a.w(i) == b.w(i);
    }
    if constexpr (E == TransportElement::Tensor){
      d = std::max(d, (a.F(i) - b.F(i)).norm());
      same = same && same_bits(a.F(i), b.F(i));
    }
  }
  return d;
}

// One step of RK4 and one of RK4cells from the same state; RK4cells' counts
template<TransportElement E>
CellsCounts step_both(std::shared_ptr<Interpol> intp, ParticleSet& a, ParticleSet& b, const double t,
                      const double dt){
  RK4Integrator rk4;
  RK4CellsIntegrator cells;
  REQUIRE(rk4_step<E>(rk4, *intp, a, t, dt).empty());
  REQUIRE(cells_step<E>(cells, *intp, b, t, dt).empty());
  return cells.counts();
}

// Steps too short to reach a facet from points well inside: the same bits
template<TransportElement E, typename I>
void check_no_crossing(std::shared_ptr<I> intp, const int D, const std::vector<bool>& per, const double t){
  intp->set_needs_gradient(true);
  intp->update(t);
  const std::vector<Vector3d> pts = inner_points(*intp, D, per, t, 300, 0.15);
  REQUIRE(pts.size() == 300);
  auto a = particles<E>(intp, pts), b = particles<E>(intp, pts);
  // Twice: the second from the cells the first step stored
  for (int s = 0; s < 2; ++s){
    const CellsCounts c = step_both<E>(intp, *a, *b, t, 1e-4);
    REQUIRE(c.steps == pts.size());
    REQUIRE(c.crossings + c.relocations + c.rejects + c.fallbacks == 0);
    bool same = false;
    largest_difference<E>(*a, *b, same);
    REQUIRE(same);
  }
}

template<typename I>
void check_no_crossing_all(std::shared_ptr<I> intp, const int D, const std::vector<bool>& per, const double t){
  check_no_crossing<TransportElement::Point>(intp, D, per, t);
  check_no_crossing<TransportElement::Vector>(intp, D, per, t);
  check_no_crossing<TransportElement::Tensor>(intp, D, per, t);
}

// Particles on facets, stored in the cell they leave; inside a cell, stored
// in a neighbour; and beside vertices, stored in a cell they are just out of.
// Crossed or relocated, then RK4's step: to round-off (the field is
// continuous) for the first two; beside a vertex the path crosses the cells
// at it within the step, kinks RK4 steps over, so to that error
template<typename I>
void check_facets_and_vertices(std::shared_ptr<I> intp, const int D, const std::vector<bool>& per, const double t){
  constexpr int nv = I::n_levels;
  intp->set_needs_gradient(true);
  intp->update(t);
  std::array<std::vector<Vector3d>, 3> pts;
  std::array<std::vector<int>, 3> cells;
  std::size_t n = 0;
  for (const Vector3d& x : inner_points(*intp, D, per, t, 600, 0.05)){
    CellPos pos;
    REQUIRE(intp->locate(x, t, pos));
    typename I::region_type R = intp->region_in(pos.id, x);
    typename I::levels_type lev, rate;
    const Vector3d xe = intp->region_point(R, x, lev);
    PointValues pv(1.);
    intp->evaluate_motion(x, t, pos, pv);
    intp->level_rates(R, pv.U, rate);
    std::array<Vector3d, nv> g;
    for (int d = 0; d < 3; ++d){
      typename I::levels_type r;
      intp->level_rates(R, Vector3d::Unit(d), r);
      for (int k = 0; k < nv; ++k) g[k][d] = r[k];
    }
    const int kind = int(n % 3), k = int((n/3) % std::size_t(nv));
    ++n;
    typename I::region_type next;
    if (kind == 0){
      // Onto facet k along its normal, inside the facet
      const Vector3d xf = xe - R.offset - lev[k]/g[k].squaredNorm()*g[k];
      typename I::levels_type l2;
      intp->region_point(R, xf, l2);
      bool clear = true;
      for (int j = 0; j < nv; ++j) if (j != k && l2[j] < 1e-3) clear = false;
      if (!clear || intp->across(R, k, xf, next) < 0) continue;
      pts[0].push_back(xf);
      cells[0].push_back(rate[k] < 0. ? R.id : next.id);
    }
    else if (kind == 1){
      if (intp->across(R, k, x, next) < 0) continue;
      pts[1].push_back(x);
      cells[1].push_back(next.id);
    }
    else {
      // Vertex k: the other levels to zero along their gradients
      Matrix3d G = Matrix3d::Zero();
      Vector3d rhs = Vector3d::Zero();
      int row = 0;
      for (int j = 0; j < nv; ++j){
        if (j == k) continue;
        G.row(row) = g[j];
        rhs[row] = -lev[j];
        ++row;
      }
      if (D == 2) G(2, 2) = 1.;
      const Vector3d v = xe - R.offset + G.colPivHouseholderQr().solve(rhs);
      // 1e-8 out of the cell across a facet at it, in the fluid
      const int j = (k + 1) % nv;
      const Vector3d y = v - 1e-8*g[j].normalized();
      // Off the walls: a particle within the band of one, moving out, falls back
      bool wall = false;
      for (int d = 0; d < D; ++d)
        wall = wall || (!per[std::size_t(d)] && (std::abs(v[d]) < 1e-6 || std::abs(v[d] - 1.) < 1e-6));
      CellPos py;
      if (wall || !intp->locate(y, t, py)) continue;
      pts[2].push_back(y);
      cells[2].push_back(R.id);
    }
  }
  for (const int kind : {0, 1, 2}){
    const std::vector<Vector3d>& p = pts[std::size_t(kind)];
    REQUIRE(p.size() > 50);
    auto a = particles<TransportElement::Point>(intp, p, cells[std::size_t(kind)]);
    auto b = particles<TransportElement::Point>(intp, p, cells[std::size_t(kind)]);
    const CellsCounts c = step_both<TransportElement::Point>(intp, *a, *b, t, 1e-3);
    INFO((kind == 0 ? "on facets" : kind == 1 ? "stored beside" : "beside vertices") << ": " << p.size()
         << " particles, immediate " << c.immediate << ", relocations " << c.relocations << ", crossings "
         << c.crossings << ", rejects " << c.rejects);
    REQUIRE(c.fallbacks == 0);
    if (kind == 0) REQUIRE(c.immediate >= p.size()*9/10);
    else REQUIRE(c.relocations == p.size());
    bool same = false;
    const double d = largest_difference<TransportElement::Point>(*a, *b, same);
    if (kind == 0) REQUIRE(d < 1e-13);
    if (kind == 1){
      REQUIRE(c.crossings == 0);
      REQUIRE(same);
    }
    if (kind == 2) REQUIRE(d < 1e-9);
  }
}

// Particles in cells at a wall of the unit box, placed on a wall facet at the
// given wall levels. rest: the facets the field at rest on the walls moves
// towards; else those the flow leaves through at a level rate below -2
template<typename I>
std::vector<std::vector<Vector3d>> wall_points(I& intp, const int D, const double t, const bool rest,
                                               const std::vector<double>& levels,
                                               std::vector<std::vector<int>>& cells){
  constexpr int nv = I::n_levels;
  std::vector<std::vector<Vector3d>> pts(levels.size());
  cells.assign(levels.size(), {});
  std::size_t n = 0;
  for (const Vector3d& x : inner_points(intp, D, {false, false, false}, t, 4000, 0.)){
    CellPos pos;
    REQUIRE(intp.locate(x, t, pos));
    const typename I::region_type R = intp.region_in(pos.id, x);
    typename I::levels_type lev;
    const Vector3d xe = intp.region_point(R, x, lev);
    for (int k = 0; k < nv; ++k){
      typename I::region_type next;
      if (intp.across(R, k, x, next) >= 0) continue;
      const std::size_t kind = n % levels.size();
      std::array<Vector3d, nv> g;
      for (int d = 0; d < 3; ++d){
        typename I::levels_type r;
        intp.level_rates(R, Vector3d::Unit(d), r);
        for (int j = 0; j < nv; ++j) g[j][d] = r[j];
      }
      const Vector3d xf = xe - R.offset - (lev[k] - levels[kind])/g[k].squaredNorm()*g[k];
      typename I::levels_type l2;
      intp.region_point(R, xf, l2);
      bool clear = true;
      for (int j = 0; j < nv; ++j) if (j != k && l2[j] < 1e-3) clear = false;
      if (!clear) continue;
      // The level's rate, from the point inside (rest: the field's sign holds to the wall)
      CellPos q;
      q.id = R.id;
      const Vector3d xr = rest ? x : Vector3d(xe - R.offset - (lev[k] - 1e-3)/g[k].squaredNorm()*g[k]);
      REQUIRE(intp.locate(xr, t, q));
      PointValues pv(1.);
      intp.evaluate_motion(xr, t, q, pv);
      typename I::levels_type rate;
      intp.level_rates(R, pv.U, rate);
      if (!(rest ? rate[k] < 0. : rate[k] < -2.)) continue;
      pts[kind].push_back(xf);
      cells[kind].push_back(R.id);
      ++n;
    }
  }
  return pts;
}

// Within the band of a wall, at rest on it and moving towards it: the held
// step, no immediate crossing and no fallback, RK4's step from above the
// wall; from below it, moved onto the facet and nudged inside. Every end in
// its cell with the wall's level at or above zero, so located
template<typename I>
void check_resting(std::shared_ptr<I> intp, const int D, const double t){
  constexpr int nv = I::n_levels;
  intp->set_needs_gradient(true);
  intp->update(t);
  std::vector<std::vector<int>> cells;
  const auto pts = wall_points(*intp, D, t, true, {5e-10, -5e-10}, cells);
  for (const std::size_t kind : {std::size_t(0), std::size_t(1)}){
    const std::vector<Vector3d>& p = pts[kind];
    REQUIRE(p.size() > 20);
    auto a = particles<TransportElement::Point>(intp, p, cells[kind]);
    auto b = particles<TransportElement::Point>(intp, p, cells[kind]);
    CellsCounts c;
    if (kind == 0){
      c = step_both<TransportElement::Point>(intp, *a, *b, t, 0.05);
      bool same = false;
      REQUIRE(largest_difference<TransportElement::Point>(*a, *b, same) < 1e-13);
    }
    else {
      RK4CellsIntegrator cs;
      REQUIRE(cells_step<TransportElement::Point>(cs, *intp, *b, t, 0.05).empty());
      c = cs.counts();
    }
    INFO((kind == 0 ? "above" : "below") << " the wall: " << p.size() << " particles, immediate " << c.immediate
         << ", relocations " << c.relocations << ", rejects " << c.rejects << ", fallbacks " << c.fallbacks);
    REQUIRE(c.immediate == 0);
    REQUIRE(c.fallbacks == 0);
    REQUIRE(c.relocations == 0);
    for (Uint i = 0; i < b->N(); ++i){
      const typename I::region_type R = intp->region_in(b->get_cell_id(i), b->x(i));
      typename I::levels_type lev;
      intp->region_point(R, b->x(i), lev);
      double wall = std::numeric_limits<double>::infinity();
      for (int k = 0; k < nv; ++k){
        typename I::region_type next;
        if (intp->across(R, k, b->x(i), next) < 0) wall = std::min(wall, lev[k]);
      }
      REQUIRE(wall >= 0.);
      // Onto the facet
      if (kind == 1) REQUIRE(wall < 1e-12);
      CellPos q;
      q.id = b->get_cell_id(i);
      REQUIRE(intp->locate(b->x(i), t, q));
    }
  }
}

// Leaving through an open boundary, from within its band or a step before
// it: the fallback, declined as under RK4
template<typename I>
void check_leaving(std::shared_ptr<I> intp, const int D, const double t){
  intp->set_needs_gradient(true);
  intp->update(t);
  std::vector<std::vector<int>> cells;
  const auto pts = wall_points(*intp, D, t, false, {5e-10, 0.05}, cells);
  for (const std::size_t kind : {std::size_t(0), std::size_t(1)}){
    const std::vector<Vector3d>& p = pts[kind];
    REQUIRE(p.size() > 10);
    auto a = particles<TransportElement::Point>(intp, p, cells[kind]);
    auto b = particles<TransportElement::Point>(intp, p, cells[kind]);
    RK4Integrator rk4;
    RK4CellsIntegrator cs;
    const std::vector<Uint> oa = rk4_step<TransportElement::Point>(rk4, *intp, *a, t, 0.1);
    const std::vector<Uint> ob = cells_step<TransportElement::Point>(cs, *intp, *b, t, 0.1);
    INFO((kind == 0 ? "within the band" : "a step before") << ": " << p.size() << " particles, declined "
         << oa.size() << " and " << ob.size());
    REQUIRE(cs.counts().fallbacks == p.size());
    REQUIRE(oa == ob);
    REQUIRE(oa.size() > p.size()/2);
    bool same = false;
    largest_difference<TransportElement::Point>(*a, *b, same);
    REQUIRE(same);
  }
}

// The split: particles on an internal plane, in either sub-cell (the start
// takes the smaller barycentric's), and on a macro facet, stored in the cell
// they leave, cross it first; particles stored in a neighbour cell are
// relocated. Then RK4's step, to round-off (the split field is continuous;
// its gradient is not, so the positions alone)
template<typename Cell>
void check_split_planes(std::shared_ptr<SplitInterpol<Cell>> intp, const int D, const std::vector<bool>& per,
                        const double t){
  constexpr int nv = Cell::n_verts;
  intp->set_needs_gradient(true);
  intp->update(t);
  std::array<std::vector<Vector3d>, 3> pts;
  std::array<std::vector<int>, 3> cells;
  std::size_t n = 0;
  for (const Vector3d& x : inner_points(*intp, D, per, t, 900, 0.05)){
    CellPos pos;
    REQUIRE(intp->locate(x, t, pos));
    const SubRegion R = intp->region_in(pos.id, x);
    std::array<double, 4> lev, rate, l2;
    const Vector3d xe = intp->region_point(R, x, lev);
    PointValues pv(1.);
    intp->evaluate_motion(x, t, pos, pv);
    intp->level_rates(R, pv.U, rate);
    const int kind = int(n % 3), k = kind == 1 ? nv - 1 : int((n/3) % std::size_t(nv - 1));
    ++n;
    SubRegion next;
    if (kind < 2){
      const Vector3d g = intp->level_grad(R, k);
      const Vector3d xf = xe - R.offset - lev[k]/g.squaredNorm()*g;
      intp->region_point(R, xf, l2);
      bool clear = true;
      for (int j = 0; j < nv; ++j) if (j != k && l2[j] < 1e-3) clear = false;
      // Through the plane, not along it: a path back across within the step meets the gradient's jump
      if (!clear || std::abs(rate[k]) < 0.2*g.norm()*pv.U.norm() || intp->across(R, k, xf, next) < 0) continue;
      pts[std::size_t(kind)].push_back(xf);
      cells[std::size_t(kind)].push_back(rate[k] < 0. ? R.id : next.id);
    }
    else {
      if (intp->across(R, nv - 1, x, next) < 0) continue;
      pts[2].push_back(x);
      cells[2].push_back(next.id);
    }
  }
  for (const int kind : {0, 1, 2}){
    const std::vector<Vector3d>& p = pts[std::size_t(kind)];
    REQUIRE(p.size() > 50);
    auto a = particles<TransportElement::Point>(intp, p, cells[std::size_t(kind)]);
    auto b = particles<TransportElement::Point>(intp, p, cells[std::size_t(kind)]);
    const CellsCounts c = step_both<TransportElement::Point>(intp, *a, *b, t, kind == 2 ? 1e-4 : 1e-3);
    INFO((kind == 0 ? "on internal planes" : kind == 1 ? "on macro facets" : "stored beside") << ": " << p.size()
         << " particles, immediate " << c.immediate << ", relocations " << c.relocations << ", crossings "
         << c.crossings << ", rejects " << c.rejects);
    REQUIRE(c.fallbacks == 0);
    bool same = false;
    const double d = largest_difference<TransportElement::Point>(*a, *b, same);
    if (kind == 0) REQUIRE(c.immediate > 0);
    if (kind == 1) REQUIRE(c.immediate >= p.size()*9/10);
    if (kind < 2) REQUIRE(d < 1e-12);
    else {
      REQUIRE(c.relocations == p.size());
      REQUIRE(c.crossings == 0);
      REQUIRE(same);
    }
  }
}

// The lattice: particles on every kind of face (between floor cells, a floor
// cell and a sub-cube, two sub-cubes, a sub-cube's mid-plane), moving through
// it, cross it first or step in the region beyond as RK4 does, to round-off
void check_lattice_planes(std::shared_ptr<StructuredInterpol> intp, const int n, const double t){
  intp->set_needs_gradient(true);
  intp->update(t);
  std::array<std::vector<Vector3d>, 4> pts;
  std::size_t m = 0;
  for (const Vector3d& x : lattice_points(*intp, n, 4000)){
    LatticeRegion R;
    std::array<double, 6> lev, rate, l2;
    REQUIRE(intp->region_of(-1, x, 0., R, lev));
    CellPos pos;
    REQUIRE(intp->locate(x, t, pos));
    PointValues pv(1.);
    intp->evaluate(x, t, pos, pv);
    intp->level_rates(R, pv.U, rate);
    const int k = int(m++ % 6);
    const Vector3d g = intp->level_grad(R, k);
    const Vector3d xf = x - lev[k]/g.squaredNorm()*g;
    intp->region_point(R, xf, l2);
    bool clear = true;
    for (int j = 0; j < 6; ++j) if (j != k && l2[j] < 2e-2) clear = false;
    LatticeRegion next;
    if (!clear || std::abs(rate[k]) < 0.2*g.norm()*pv.U.norm() || intp->across(R, k, xf, next) < 0) continue;
    pts[next.f == R.f ? 3 : std::size_t(R.sub >= 0) + std::size_t(next.sub >= 0)].push_back(xf);
  }
  for (std::size_t kind = 0; kind < 4; ++kind){
    const std::vector<Vector3d>& p = pts[kind];
    REQUIRE(p.size() > 20);
    auto a = particles<TransportElement::Point>(intp, p);
    auto b = particles<TransportElement::Point>(intp, p);
    const CellsCounts c = step_both<TransportElement::Point>(intp, *a, *b, t, 1e-3);
    INFO((kind == 0 ? "between floor cells" : kind == 1 ? "floor cell and sub-cube" : kind == 2 ? "between sub-cubes"
          : "mid-planes") << ": " << p.size() << " particles, immediate " << c.immediate << ", crossings "
         << c.crossings << ", rejects " << c.rejects);
    REQUIRE(c.fallbacks == 0);
    REQUIRE(c.crossings > 0);
    bool same = false;
    REQUIRE(largest_difference<TransportElement::Point>(*a, *b, same) < 1e-12);
  }
}

// At a wall of the lattice (a sub-cube's mid-plane before a solid node),
// where the flow comes to rest: particles within the band of it, above and
// below, step held, never cross before a step, never fall back, and end in
// the fluid, a point below the wall moved onto it and inside
void check_lattice_resting(std::shared_ptr<StructuredInterpol> intp, const int n, const double t){
  intp->set_needs_gradient(true);
  intp->update(t);
  std::array<std::vector<Vector3d>, 2> pts;
  std::array<std::vector<int>, 2> ids;
  std::size_t m = 0;
  for (const Vector3d& x : lattice_points(*intp, n, 20000)){
    LatticeRegion R;
    std::array<double, 6> lev, rate, l2;
    REQUIRE(intp->region_of(-1, x, 0., R, lev));
    for (int k = 0; k < 6; ++k){
      if (!intp->is_wall(R, k)) continue;
      const Vector3d g = intp->level_grad(R, k);
      // Moving towards the wall a little inside
      CellPos pos;
      const Vector3d xi = x - (lev[k] - 1e-3)/g.squaredNorm()*g;
      if (!intp->locate(xi, t, pos)) continue;
      PointValues pv(1.);
      intp->evaluate(xi, t, pos, pv);
      intp->level_rates(R, pv.U, rate);
      if (!(rate[k] < 0.)) continue;
      const std::size_t kind = m++ % 2;
      const Vector3d xf = x - (lev[k] - (kind ? -5e-10 : 5e-10))/g.squaredNorm()*g;
      intp->region_point(R, xf, l2);
      bool clear = true;
      for (int j = 0; j < 6; ++j) if (j != k && l2[j] < 1e-3) clear = false;
      if (!clear) continue;
      pts[kind].push_back(xf);
      ids[kind].push_back(R.id);
    }
  }
  for (const std::size_t kind : {std::size_t(0), std::size_t(1)}){
    const std::vector<Vector3d>& p = pts[kind];
    REQUIRE(p.size() > 20);
    auto a = particles<TransportElement::Point>(intp, p, ids[kind]);
    auto b = particles<TransportElement::Point>(intp, p, ids[kind]);
    CellsCounts c;
    if (kind == 0){
      c = step_both<TransportElement::Point>(intp, *a, *b, t, 0.5);
      bool same = false;
      REQUIRE(largest_difference<TransportElement::Point>(*a, *b, same) < 1e-12);
    }
    else {
      RK4CellsIntegrator cs;
      REQUIRE(cells_step<TransportElement::Point>(cs, *intp, *b, t, 0.5).empty());
      c = cs.counts();
    }
    INFO((kind == 0 ? "above" : "below") << " the wall: " << p.size() << " particles, immediate " << c.immediate
         << ", relocations " << c.relocations << ", rejects " << c.rejects << ", fallbacks " << c.fallbacks);
    REQUIRE(c.immediate == 0);
    REQUIRE(c.fallbacks == 0);
    REQUIRE(c.relocations == 0);
    for (Uint i = 0; i < b->N(); ++i){
      CellPos q;
      REQUIRE(intp->locate(b->x(i), t, q));
      const LatticeRegion R = intp->region_in(b->get_cell_id(i), b->x(i));
      std::array<double, 6> lev;
      intp->region_point(R, b->x(i), lev);
      for (int k = 0; k < 6; ++k) if (intp->is_wall(R, k)) REQUIRE(lev[k] > 0.);
    }
  }
}


// The first event of a major step of dt from x stored in id at t, by the
// kernel's own functions: the facet predicted, its time, and, landed on it
// with the major step's end out of the way, the time taken and whether the
// region entered holds the landing
struct Landing {
  int k = -1;
  double tau = 0., taken = 0.;
  bool landed = false, held = false;
};

template<typename I>
Landing landing_of(I& intp, const Vector3d& x, const int id, const double t, const double dt){
  PointValues pv(intp.get_U0());
  HeldEval<I> ev(intp, pv);
  ev.band = cells::band;
  cells::State<TransportElement::Point, I> s;
  s.x = x;
  s.t = t;
  s.r = dt;
  s.R = intp.region_in(id, x);
  ev.bind(s.R);
  Landing L;
  if (!cells::first_stage(intp, ev, pv, s)) return L;
  L.tau = cells::first_crossing<I::n_levels>(intp, s.R, s.lev, s.d1, s.d2, dt, -1, cells::t_eps*dt, L.k);
  if (!(L.tau > 0. && L.tau < cells::inf)) return L;
  CellsCounts c;
  s.r = 4.*dt;
  if (!cells::step_to_facet(intp, ev, pv, s, L.k, L.tau, cells::t_eps*dt, c) || c.crossings != 1) return L;
  L.landed = true;
  L.taken = s.t - t;
  L.held = c.relocations == 0;
  return L;
}

// Landings that end the major step: a major step of exactly the time a
// landing takes stores a region that holds the end, the entered one or, where
// the landing left it near an edge, the one relocated to. A major step ending
// between the predicted crossing and the landing: the landing is rejected and
// the step taken to the end without the crossing, RK4's step to round-off
template<typename I>
void check_landings(std::shared_ptr<I> intp, const int D, const std::vector<bool>& per, const double t,
                    const double dt){
  constexpr int nv = I::n_levels;
  intp->set_needs_gradient(true);
  intp->update(t);
  std::vector<Vector3d> out_pts, past_pts;
  std::vector<int> out_ids, past_ids;
  std::vector<double> out_dt, past_dt;
  for (const Vector3d& x : inner_points(*intp, D, per, t, 6000, 0.01)){
    CellPos pos;
    REQUIRE(intp->locate(x, t, pos));
    const Landing L = landing_of(*intp, x, pos.id, t, dt);
    if (!L.landed || !(L.taken > L.tau*(1. + 1e-9))) continue;
    if (!L.held){
      out_pts.push_back(x);
      out_ids.push_back(pos.id);
      out_dt.push_back(L.taken);
    }
    // Well before the facet at the end
    else if (L.taken > L.tau*1.001 && past_pts.size() < 40){
      past_pts.push_back(x);
      past_ids.push_back(pos.id);
      past_dt.push_back(L.tau + 0.5*(L.taken - L.tau));
    }
  }
  INFO(out_pts.size() << " landings out of the region entered, " << past_pts.size() << " past the end");
  REQUIRE(out_pts.size() >= 3);
  REQUIRE(past_pts.size() >= 20);
  for (std::size_t i = 0; i < out_pts.size(); ++i){
    auto b = particles<TransportElement::Point>(intp, {out_pts[i]}, {out_ids[i]});
    RK4CellsIntegrator cs;
    REQUIRE(cells_step<TransportElement::Point>(cs, *intp, *b, t, out_dt[i]).empty());
    const CellsCounts c = cs.counts();
    INFO("ending on a landing: crossings " << c.crossings << ", relocations " << c.relocations << ", rejects "
         << c.rejects << ", fallbacks " << c.fallbacks);
    REQUIRE(c.crossings == 1);
    REQUIRE(c.fallbacks == 0);
    typename I::levels_type lev;
    intp->region_point(intp->region_in(b->get_cell_id(0), b->x(0)), b->x(0), lev);
    REQUIRE(min_level<nv>(lev) >= -cells::band);
  }
  // One step to the end after the rejection: RK4's
  std::size_t one = 0;
  for (std::size_t i = 0; i < past_pts.size(); ++i){
    auto a = particles<TransportElement::Point>(intp, {past_pts[i]}, {past_ids[i]});
    auto b = particles<TransportElement::Point>(intp, {past_pts[i]}, {past_ids[i]});
    const CellsCounts c = step_both<TransportElement::Point>(intp, *a, *b, t, past_dt[i]);
    INFO("landing past the end: steps " << c.steps << ", crossings " << c.crossings << ", rejects " << c.rejects
         << ", fallbacks " << c.fallbacks);
    REQUIRE(c.crossings == 0);
    REQUIRE(c.rejects >= 1);
    REQUIRE(c.fallbacks == 0);
    REQUIRE(b->get_cell_id(0) == past_ids[i]);
    if (c.steps != 2) continue;
    ++one;
    bool same = false;
    REQUIRE(largest_difference<TransportElement::Point>(*a, *b, same) < 1e-12);
  }
  REQUIRE(one >= past_pts.size()/2);
}

// x -> 1 - (1 - x)^2: columns narrowing towards x = 1
inline Vector3d graded_x(const Vector3d& x){ return Vector3d(1. - (1. - x[0])*(1. - x[0]), x[1], x[2]); }

// Particles within the band of a facet into a narrower cell, stored in the
// wider one they leave: the entered facet's level is below the band there.
// Crossed once each, then RK4's step to round-off; never back and forth
template<typename I>
void check_into_narrower(std::shared_ptr<I> intp, const int D, const double t){
  constexpr int nv = I::n_levels;
  intp->set_needs_gradient(true);
  intp->update(t);
  std::vector<Vector3d> pts;
  std::vector<int> ids;
  for (const Vector3d& x : inner_points(*intp, D, {false, false, false}, t, 2000, 0.05)){
    CellPos pos;
    REQUIRE(intp->locate(x, t, pos));
    const typename I::region_type R = intp->region_in(pos.id, x);
    typename I::levels_type lev, rate, l2, l3;
    const Vector3d xe = intp->region_point(R, x, lev);
    PointValues pv(1.);
    intp->evaluate_motion(x, t, pos, pv);
    intp->level_rates(R, pv.U, rate);
    for (int k = 0; k < nv; ++k){
      if (!(rate[k] < 0.)) continue;
      const Vector3d g = intp->level_grad(R, k);
      const Vector3d xf = xe - R.offset - (lev[k] - 0.9*cells::band)/g.squaredNorm()*g;
      typename I::region_type next;
      const int e = intp->across(R, k, xf, next);
      if (e < 0) continue;
      intp->region_point(R, xf, l2);
      intp->region_point(next, xf, l3);
      bool clear = true;
      for (int j = 0; j < nv; ++j) if (j != k && l2[j] < 1e-3) clear = false;
      if (!clear || !(l3[e] < -cells::band)) continue;
      pts.push_back(xf);
      ids.push_back(R.id);
    }
  }
  REQUIRE(pts.size() > 20);
  auto a = particles<TransportElement::Point>(intp, pts, ids);
  auto b = particles<TransportElement::Point>(intp, pts, ids);
  const CellsCounts c = step_both<TransportElement::Point>(intp, *a, *b, t, 1e-3);
  INFO(pts.size() << " particles, immediate " << c.immediate << ", relocations " << c.relocations
       << ", fallbacks " << c.fallbacks);
  REQUIRE(c.fallbacks == 0);
  REQUIRE(c.relocations == 0);
  REQUIRE(c.immediate == pts.size());
  bool same = false;
  REQUIRE(largest_difference<TransportElement::Point>(*a, *b, same) < 1e-13);
}

// Past the cap on steps a major step: the fallback, RK4's substeps from the
// start. A path through many narrow cells into a wide one passes the wide
// cell's cap
template<typename I>
void check_cap(std::shared_ptr<I> intp, const double t, const double dt){
  intp->set_needs_gradient(true);
  intp->update(t);
  std::mt19937 rng(5);
  std::uniform_real_distribution<double> ux(0.001, 0.02), uy(0.3, 0.7);
  std::vector<Vector3d> pts;
  while (pts.size() < 30){
    const Vector3d x(ux(rng), uy(rng), 0.);
    CellPos pos;
    if (intp->locate(x, t, pos)) pts.push_back(x);
  }
  auto a = particles<TransportElement::Point>(intp, pts);
  auto b = particles<TransportElement::Point>(intp, pts);
  RK4CellsIntegrator cs;
  REQUIRE(cells_step<TransportElement::Point>(cs, *intp, *b, t, dt).empty());
  for (Uint i = 0; i < a->N(); ++i)
    REQUIRE(rk4_substeps<TransportElement::Point>(*intp, *a, i, t, dt));
  INFO("crossings " << cs.counts().crossings << ", fallbacks " << cs.counts().fallbacks);
  REQUIRE(cs.counts().fallbacks == pts.size());
  bool same = false;
  largest_difference<TransportElement::Point>(*a, *b, same);
  REQUIRE(same);
}

// Two walls of one cell at an acute angle (the unit square sheared): ends
// below both, or on one and below the other, are moved onto both and inside,
// a displacement of the band's order, and located
template<typename I>
void check_acute_corner(std::shared_ptr<I> intp, const double t){
  constexpr int nv = I::n_levels;
  intp->set_needs_gradient(true);
  intp->update(t);
  const std::vector<std::array<double, 2>> targets = {{-5e-10, -5e-10}, {0., -5e-10}, {-5e-10, 0.}, {-9e-10, -1e-10}};
  std::vector<Vector3d> pts;
  std::vector<int> ids;
  std::vector<std::array<int, 2>> walls;
  std::vector<int> seen;
  for (const Vector3d& x : inner_points(*intp, 2, {false, false, false}, t, 4000, 0.)){
    CellPos pos;
    REQUIRE(intp->locate(x, t, pos));
    if (std::find(seen.begin(), seen.end(), pos.id) != seen.end()) continue;
    const typename I::region_type R = intp->region_in(pos.id, x);
    std::vector<int> w;
    for (int k = 0; k < nv; ++k) if (intp->is_wall(R, k)) w.push_back(k);
    if (w.size() != 2 || !(intp->level_grad(R, w[0]).dot(intp->level_grad(R, w[1])) < 0.)) continue;
    seen.push_back(pos.id);
    typename I::levels_type lev;
    intp->region_point(R, x, lev);
    Eigen::Matrix2d G;
    for (int j = 0; j < 2; ++j) G.row(j) = intp->level_grad(R, w[std::size_t(j)]).template head<2>();
    for (const auto& tg : targets){
      const Eigen::Vector2d d = G.inverse()*Eigen::Vector2d(tg[0] - lev[w[0]], tg[1] - lev[w[1]]);
      pts.push_back(x + Vector3d(d[0], d[1], 0.));
      ids.push_back(pos.id);
      walls.push_back({w[0], w[1]});
    }
  }
  REQUIRE(pts.size() >= 8);
  auto b = particles<TransportElement::Point>(intp, pts, ids);
  RK4CellsIntegrator cs;
  REQUIRE(cells_step<TransportElement::Point>(cs, *intp, *b, t, 1e-3).empty());
  INFO(pts.size() << " particles, fallbacks " << cs.counts().fallbacks);
  REQUIRE(cs.counts().fallbacks == 0);
  for (Uint i = 0; i < b->N(); ++i){
    REQUIRE(b->get_cell_id(i) == ids[i]);
    typename I::levels_type lev;
    intp->region_point(intp->region_in(ids[i], b->x(i)), b->x(i), lev);
    INFO("particle " << i << ": wall levels " << lev[walls[i][0]] << ", " << lev[walls[i][1]]);
    REQUIRE(lev[walls[i][0]] >= 0.);
    REQUIRE(lev[walls[i][1]] >= 0.);
    REQUIRE((b->x(i) - pts[i]).norm() < 1e-8);
    CellPos q;
    q.id = ids[i];
    REQUIRE(intp->locate(b->x(i), t, q));
  }
}

// A checkpoint on the lattice keeps the regions' node ids, every one of the
// lattice's, and drops the rest. Particles on a sub-cube's mid-plane, stored
// in the lower sub-cube they move into, start there, not in the upper one a
// point on the plane is located in: resumed from the checkpoint, they step as
// they would have, the same events and bits
void check_lattice_restart(std::shared_ptr<StructuredInterpol> intp, const int n, const double t){
  intp->set_needs_gradient(true);
  intp->update(t);
  std::vector<Vector3d> pts;
  std::vector<int> ids;
  std::size_t m = 0;
  for (Vector3d x : lattice_points(*intp, n, 20000)){
    LatticeRegion R;
    std::array<double, 6> lev;
    REQUIRE(intp->region_of(-1, x, 0., R, lev));
    if (R.sub < 0) continue;
    const int a = int(m++ % 3);
    x[a] = R.f[a] + 0.5;
    // The lower sub-cube's node: the region of a point just below the plane
    Vector3d y = x;
    y[a] -= 1e-6;
    LatticeRegion below;
    if (!intp->region_of(-1, y, 0., below, lev) || below.sub != (R.sub & ~(1 << a)) || below.f != R.f) continue;
    CellPos pos;
    if (!intp->locate(y, t, pos)) continue;
    PointValues pv(1.);
    intp->evaluate(x, t, pos, pv);
    if (!(pv.U[a] < -0.1)) continue;
    pts.push_back(x);
    ids.push_back(below.id);
    if (pts.size() == 60) break;
  }
  REQUIRE(pts.size() >= 30);
  auto a = particles<TransportElement::Point>(intp, pts, ids);
  // The last particle at the largest node id, written and read alone
  auto edge = particles<TransportElement::Point>(intp, {pts[0], pts[1]}, {n*n*n - 1, n*n*n});
  CaseDir dir("cells_lattice_restart");
  const auto round_trip = [&](const ParticleSet& ps){
    {
      H5::H5File h5(dir.file("checkpoint.h5"), H5F_ACC_TRUNC);
      ps.write_checkpoint(h5, true);
    }
    auto back = std::make_unique<ParticleSet>(intp, ps.N());
    H5::H5File h5(dir.file("checkpoint.h5"), H5F_ACC_RDONLY);
    back->read_checkpoint(h5, true);
    return back;
  };
  auto e = round_trip(*edge);
  REQUIRE(e->get_cell_id(0) == n*n*n - 1);
  REQUIRE(e->get_cell_id(1) == -1);
  auto b = round_trip(*a);
  for (Uint i = 0; i < b->N(); ++i) REQUIRE(b->get_cell_id(i) == ids[i]);
  RK4CellsIntegrator ca, cb;
  REQUIRE(cells_step<TransportElement::Point>(ca, *intp, *a, t, 0.5).empty());
  REQUIRE(cells_step<TransportElement::Point>(cb, *intp, *b, t, 0.5).empty());
  INFO("immediate crossings " << ca.counts().immediate << " and " << cb.counts().immediate << ", fallbacks "
       << ca.counts().fallbacks << " and " << cb.counts().fallbacks);
  REQUIRE(ca.counts().immediate == cb.counts().immediate);
  REQUIRE(ca.counts().crossings == cb.counts().crossings);
  REQUIRE(ca.counts().fallbacks == cb.counts().fallbacks);
  for (Uint i = 0; i < a->N(); ++i){
    REQUIRE(same_bits(a->x(i), b->x(i)));
    REQUIRE(a->get_cell_id(i) == b->get_cell_id(i));
  }
}

// The 90th percentile
double p90(std::vector<double> v){
  std::sort(v.begin(), v.end());
  return v[std::size_t(0.9*double(v.size() - 1))];
}

// Per particle, the difference of position and of the carried element between two sets
template<TransportElement E>
void differences(const ParticleSet& a, const ParticleSet& b, std::vector<double>& dx, std::vector<double>& de){
  dx.clear();
  de.clear();
  for (Uint i = 0; i < a.N(); ++i){
    dx.push_back((a.x(i) - b.x(i)).norm());
    if constexpr (E == TransportElement::Vector)
      de.push_back((a.rhohat(i) - b.rhohat(i)).norm() + std::abs(a.w(i) - b.w(i)));
    if constexpr (E == TransportElement::Tensor) de.push_back((a.F(i) - b.F(i)).norm());
  }
}

// Points well inside regions and inside the box [lo, hi]^D
template<typename I>
std::vector<Vector3d> central_points(I& intp, const int D, const double t, const std::size_t n, const double lo,
                                     const double hi){
  std::vector<Vector3d> out;
  for (const Vector3d& x : inner_points(intp, D, {false, false, false}, t, 20*n, 0.01)){
    bool in = true;
    for (int d = 0; d < D; ++d) in = in && x[d] > lo && x[d] < hi;
    if (in) out.push_back(x);
    if (out.size() == n) break;
  }
  return out;
}

// RK4cells over T at dt, dt/2, dt/4 and dt/8 from the same points: steps
// that land on facets predicted ahead, none falling back, and the changes of
// position and element falling by about 16 a halving (p90)
template<TransportElement E, typename I>
void check_order(std::shared_ptr<I> intp, const std::vector<Vector3d>& pts, const double t, const double T,
                 const double dt){
  intp->set_needs_gradient(true);
  intp->update(t);
  std::vector<std::unique_ptr<ParticleSet>> runs;
  for (int k = 0; k < 4; ++k){
    auto ps = particles<E>(intp, pts);
    RK4CellsIntegrator cs;
    const double h = dt/double(1 << k);
    const int n = int(std::lround(T/h));
    for (int j = 0; j < n; ++j) REQUIRE(cells_step<E>(cs, *intp, *ps, t + j*h, h).empty());
    const CellsCounts c = cs.counts();
    INFO("dt " << h << ": " << pts.size() << " particles, crossings " << c.crossings << ", immediate "
         << c.immediate << ", rejects " << c.rejects << ", relocations " << c.relocations);
    REQUIRE(c.fallbacks == 0);
    REQUIRE(c.crossings - c.immediate > pts.size());
    REQUIRE(c.rejects < c.crossings/10);
    runs.push_back(std::move(ps));
  }
  std::vector<double> rx, re;
  for (int k = 0; k < 2; ++k){
    std::vector<double> dx0, de0, dx1, de1;
    differences<E>(*runs[k], *runs[k + 1], dx0, de0);
    differences<E>(*runs[k + 1], *runs[k + 2], dx1, de1);
    rx.push_back(p90(dx0)/p90(dx1));
    if (!de0.empty()) re.push_back(p90(de0)/p90(de1));
  }
  INFO("ratios of the position's changes " << rx[0] << ", " << rx[1]);
  for (const double r : rx) REQUIRE(r > 12.);
  INFO("ratios of the element's changes " << (re.empty() ? 0. : re[0]) << ", " << (re.empty() ? 0. : re[1]));
  for (const double r : re) REQUIRE(r > 12.);
}

// Of the points, the first in a cell with two walls at an acute angle: its cell, the point and the walls
template<typename I>
bool acute_cell(I& intp, const std::vector<Vector3d>& pts, const double t, int& id, Vector3d& x,
                std::array<int, 2>& w){
  for (const Vector3d& y : pts){
    CellPos pos;
    REQUIRE(intp.locate(y, t, pos));
    const typename I::region_type R = intp.region_in(pos.id, y);
    std::vector<int> ws;
    for (int k = 0; k < I::n_levels; ++k) if (intp.is_wall(R, k)) ws.push_back(k);
    if (ws.size() != 2 || !(intp.level_grad(R, ws[0]).dot(intp.level_grad(R, ws[1])) < 0.)) continue;
    id = pos.id;
    x = y;
    w = {ws[0], ws[1]};
    return true;
  }
  return false;
}

// The point nearest x of cell id where levels a and b take the values la and lb
template<typename I>
Vector3d at_levels(I& intp, const int id, const Vector3d& x, const int a, const int b, const double la,
                   const double lb){
  const typename I::region_type R = intp.region_in(id, x);
  typename I::levels_type lev;
  intp.region_point(R, x, lev);
  Eigen::Matrix<double, 2, 3> G;
  G.row(0) = intp.level_grad(R, a);
  G.row(1) = intp.level_grad(R, b);
  return x + G.transpose()*(G*G.transpose()).inverse()*Eigen::Vector2d(la - lev[a], lb - lev[b]);
}

// onto_walls on a state at y in cell id; its result, the state's point and levels after
template<typename I>
bool onto(I& intp, const int id, const Vector3d& y, const int mask, Vector3d& x, typename I::levels_type& lev){
  cells::State<TransportElement::Point, I> s;
  s.R = intp.region_in(id, y);
  s.x = y;
  intp.region_point(s.R, s.x, lev);
  const bool ok = cells::onto_walls(intp, s, mask, lev);
  x = s.x;
  return ok;
}

// Points just below two walls at an acute angle are moved onto both (an edge
// of a tet) and inside, often only on a second try with a nudge (a level
// rounded below zero). It fails where the walls do not meet (all of a cell's
// facets), where the end is left below another facet by more than the band,
// and where a level moves more than proj_max; the first before moving the point
template<typename I>
void check_onto_walls(std::shared_ptr<I> intp, const int D, const double t){
  constexpr int nv = I::n_levels;
  intp->update(t);
  int id = -1;
  Vector3d x0, x;
  std::array<int, 2> w = {-1, -1};
  REQUIRE(acute_cell(*intp, inner_points(*intp, D, {false, false, false}, t, 4000, 0.), t, id, x0, w));
  const int mask = (1 << w[0]) | (1 << w[1]);
  typename I::levels_type lev;
  std::mt19937 rng(3);
  std::uniform_real_distribution<double> below(-cells::band, 0.);
  for (int i = 0; i < 400; ++i){
    const Vector3d y = at_levels(*intp, id, x0, w[0], w[1], below(rng), below(rng));
    INFO("point " << i);
    REQUIRE(onto(*intp, id, y, mask, x, lev));
    for (const int k : w){
      REQUIRE(lev[k] >= 0.);
      REQUIRE(lev[k] <= cells::band);
    }
    REQUIRE((x - y).norm() < 1e-8);
  }
  // All facets of a cell do not meet: untouched
  const Vector3d y = at_levels(*intp, id, x0, w[0], w[1], -5e-10, -5e-10);
  REQUIRE(!onto(*intp, id, y, (1 << nv) - 1, x, lev));
  REQUIRE(same_bits(x, y));
  // Onto one wall, the other left 5e-9 below
  REQUIRE(!onto(*intp, id, at_levels(*intp, id, x0, w[0], w[1], -5e-10, -5e-9), 1 << w[0], x, lev));
  REQUIRE(lev[w[0]] >= 0.);
  REQUIRE(lev[w[1]] < -cells::band);
  // A wall 2e-6 below
  REQUIRE(!onto(*intp, id, at_levels(*intp, id, x0, w[0], w[1], -2e-6, 0.1), 1 << w[0], x, lev));
  REQUIRE(lev[w[0]] >= 0.);
  REQUIRE(min_level<nv>(lev) >= 0.);
}

// A mesh 1e8 from the origin: a point below a facet by less than a step of
// the levels between neighbouring doubles is not moved by any nudge up to the
// band, so the projection gives up, the point untouched
template<typename I>
void check_onto_walls_far(std::shared_ptr<I> intp, const double t){
  constexpr int nv = I::n_levels;
  intp->update(t);
  const std::vector<Vector3d> pts = inner_points(*intp, 2, {false, false, false}, t, 1, 0.1);
  REQUIRE(pts.size() == 1);
  CellPos pos;
  REQUIRE(intp->locate(pts[0], t, pos));
  const typename I::region_type R = intp->region_in(pos.id, pts[0]);
  // A slanted facet
  int k = -1;
  for (int m = 0; m < nv && k < 0; ++m){
    const Vector3d g = intp->level_grad(R, m);
    if (std::abs(g[0]) > 0.2*g.norm() && std::abs(g[1]) > 0.2*g.norm()) k = m;
  }
  REQUIRE(k >= 0);
  const Vector3d g = intp->level_grad(R, k);
  typename I::levels_type lev;
  intp->region_point(R, pts[0], lev);
  const Vector3d xf = pts[0] - lev[k]/g.squaredNorm()*g;
  // Doubles about the facet: one below it that no nudge up to 4 bands moves
  bool found = false;
  Vector3d y;
  for (int i = -20; i <= 20 && !found; ++i)
    for (int j = -20; j <= 20 && !found; ++j){
      y = xf;
      for (int a = 0; a < 2; ++a){
        const int n = a ? j : i;
        for (int q = 0; q < std::abs(n); ++q) y[a] = std::nextafter(y[a], n > 0 ? 1e300 : -1e300);
      }
      intp->region_point(R, y, lev);
      if (!(lev[k] < 0.)) continue;
      found = true;
      for (int a = 0; a < 2; ++a){
        const double ulp = std::nextafter(y[a], 1e300) - y[a];
        found = found && (4.*cells::band - lev[k])*std::abs(g[a])/g.squaredNorm() < 0.25*ulp;
      }
    }
  REQUIRE(found);
  Vector3d x;
  REQUIRE(!onto(*intp, pos.id, y, 1 << k, x, lev));
  REQUIRE(same_bits(x, y));
}

// Ends below both walls of a corner 1e-4 rad wide: their normals are parallel
// to round-off, so the walls do not meet for the joint projection, which
// fails, and the major step falls back, RK4's substeps from the start
template<typename I>
void check_projection_fallback(std::shared_ptr<I> intp, const std::vector<Vector3d>& corner, const double t){
  intp->set_needs_gradient(true);
  intp->update(t);
  int id = -1;
  Vector3d x0;
  std::array<int, 2> w = {-1, -1};
  REQUIRE(acute_cell(*intp, corner, t, id, x0, w));
  std::vector<Vector3d> pts;
  for (const auto& lv : std::vector<std::array<double, 2>>{{0., -5e-10}, {-5e-10, 0.}, {-5e-10, -5e-10}, {-9e-10, -1e-10}}){
    pts.push_back(at_levels(*intp, id, x0, w[0], w[1], lv[0], lv[1]));
    Vector3d x;
    typename I::levels_type lev;
    REQUIRE(!onto(*intp, id, pts.back(), (1 << w[0]) | (1 << w[1]), x, lev));
    REQUIRE(same_bits(x, pts.back()));
  }
  auto a = particles<TransportElement::Point>(intp, pts, std::vector<int>(pts.size(), id));
  auto b = particles<TransportElement::Point>(intp, pts, std::vector<int>(pts.size(), id));
  RK4CellsIntegrator cs;
  const std::vector<Uint> out = cells_step<TransportElement::Point>(cs, *intp, *b, t, 1e-3);
  std::vector<Uint> out_rk4;
  for (Uint i = 0; i < a->N(); ++i)
    if (!rk4_substeps<TransportElement::Point>(*intp, *a, i, t, 1e-3)) out_rk4.push_back(i);
  INFO("fallbacks " << cs.counts().fallbacks << ", declined " << out.size());
  REQUIRE(cs.counts().fallbacks == pts.size());
  REQUIRE(out == out_rk4);
  // Every start below the walls: declined, and left at the start, not at the step's end
  REQUIRE(out.size() == pts.size());
  for (Uint i = 0; i < b->N(); ++i)
    REQUIRE(same_bits(b->x(i), pts[i]));
  bool same = false;
  largest_difference<TransportElement::Point>(*a, *b, same);
  REQUIRE(same);
}

}  // namespace

TEST_CASE("At a wall: resting particles step held and end inside, leaving ones fall back", "[cells]"){
  SECTION("P2 triangles at rest on the walls"){
    CaseDir c("cells_rest_tri");
    write_stamped<Triangle>(c, 6, "P2", {false, false, false}, true);
    check_resting(std::make_shared<SimplexInterpol<Triangle>>(c.file("h5_params.dat")), 2, 0.3);
  }
  SECTION("P2 tets at rest on the walls"){
    CaseDir c("cells_rest_tet");
    write_stamped<Tet>(c, 3, "P2", {false, false, false}, true);
    check_resting(std::make_shared<SimplexInterpol<Tet>>(c.file("h5_params.dat")), 3, 0.3);
  }
  SECTION("P2 triangles, open boundaries"){
    CaseDir c("cells_open_tri");
    write_stamped<Triangle>(c, 6, "P2", {false, false, false}, false);
    check_leaving(std::make_shared<SimplexInterpol<Triangle>>(c.file("h5_params.dat")), 2, 0.3);
  }
}

TEST_CASE("A step that meets no facet is RK4's step bit for bit", "[cells]"){
  SECTION("P2 tets, periodic, between stamps"){
    CaseDir c("cells_p2_tet");
    write_stamped<Tet>(c, 3, "P2", {true, true, true}, false);
    check_no_crossing_all(std::make_shared<SimplexInterpol<Tet>>(c.file("h5_params.dat")), 3, {true, true, true}, 0.3);
  }
  SECTION("P2 triangles, periodic in x"){
    CaseDir c("cells_p2_tri");
    write_stamped<Triangle>(c, 6, "P2", {true, false, false}, false);
    check_no_crossing_all(std::make_shared<SimplexInterpol<Triangle>>(c.file("h5_params.dat")), 2,
                          {true, false, false}, 0.7);
  }
  SECTION("near-wall P2 tets, XDMF"){
    CaseDir c("cells_wall_tet");
    write_stamped<Tet>(c, 3, "P1", {false, false, false}, true);
    std::ofstream(c.file("xdmf_params.dat"), std::ios::app) << "wall_p2=edge\n";
    check_no_crossing_all(std::make_shared<XDMFInterpol<Tet>>(c.file("xdmf_params.dat")), 3,
                          {false, false, false}, 0.3);
  }
  SECTION("near-wall P2 triangles, XDMF"){
    CaseDir c("cells_wall_tri");
    write_stamped<Triangle>(c, 6, "P1", {false, false, false}, true);
    std::ofstream(c.file("xdmf_params.dat"), std::ios::app) << "wall_p2=edge\n";
    check_no_crossing_all(std::make_shared<XDMFInterpol<Triangle>>(c.file("xdmf_params.dat")), 2,
                          {false, false, false}, 0.3);
  }
  SECTION("frequency tets"){
    CaseDir c("cells_freq_tet");
    write_freq<Tet>(c, 3);
    check_no_crossing_all(std::make_shared<TetFreqInterpol>(c.params()), 3, {true, true, true}, 0.25);
  }
  SECTION("frequency triangles"){
    CaseDir c("cells_freq_tri");
    write_freq<Triangle>(c, 6);
    check_no_crossing_all(std::make_shared<TriangleFreqInterpol>(c.params()), 2, {true, true, false}, 0.25);
  }
}

TEST_CASE("The split: a step in one sub-cell is RK4's, planes are crossed, stale cells relocated", "[cells][split]"){
  SECTION("tets, periodic in x, between stamps"){
    CaseDir c("cells_split_tet");
    write_split<Tet>(c, 3, true);
    auto intp = std::make_shared<SplitInterpol<Tet>>(c.file("h5_params.dat"));
    check_no_crossing_all(intp, 3, {true, false, false}, 0.3);
    check_split_planes(intp, 3, {true, false, false}, 0.3);
  }
  SECTION("triangles, walls"){
    CaseDir c("cells_split_tri");
    write_split<Triangle>(c, 6, false);
    auto intp = std::make_shared<SplitInterpol<Triangle>>(c.file("h5_params.dat"));
    check_no_crossing_all(intp, 2, {false, false, false}, 0.7);
    check_split_planes(intp, 2, {false, false, false}, 0.7);
  }
}

TEST_CASE("The lattice: a step in one region is RK4's, faces are crossed, walls are not", "[cells][lattice]"){
  const int n = 12;
  CaseDir c("cells_lattice");
  write_lattice(c, n);
  auto intp = std::make_shared<StructuredInterpol>(c.file("felbm_params.dat"));
  check_no_crossing_all(intp, 3, {true, true, true}, 0.3);
  check_lattice_planes(intp, n, 0.3);
  check_lattice_resting(intp, n, 0.7);
  check_lattice_restart(intp, n, 0.3);
}

TEST_CASE("A particle on a facet crosses first, one beside its cell is relocated", "[cells]"){
  SECTION("P2 tets, periodic"){
    CaseDir c("cells_facet_tet");
    write_stamped<Tet>(c, 3, "P2", {true, true, true}, false);
    check_facets_and_vertices(std::make_shared<SimplexInterpol<Tet>>(c.file("h5_params.dat")), 3,
                              {true, true, true}, 0.3);
  }
  SECTION("P2 triangles, at rest on walls"){
    CaseDir c("cells_facet_tri");
    write_stamped<Triangle>(c, 6, "P2", {false, false, false}, true);
    check_facets_and_vertices(std::make_shared<SimplexInterpol<Triangle>>(c.file("h5_params.dat")), 2,
                              {false, false, false}, 0.3);
  }
}

TEST_CASE("Moved particles lose their cell, inserted ones take the cell located, under every scheme", "[cells]"){
  CaseDir c("cells_keep_tri");
  write_stamped<Triangle>(c, 6, "P2", {false, false, false}, true);
  auto intp = std::make_shared<SimplexInterpol<Triangle>>(c.file("h5_params.dat"));
  intp->update(0.);
  const std::vector<Vector3d> pts = inner_points(*intp, 2, {false, false, false}, 0., 4, 0.05);
  REQUIRE(pts.size() == 4);
  // No scheme asks for it
  ParticleSet ps(intp, 10);
  ps.add(pts, 0);
  for (Uint i = 0; i < ps.N(); ++i){
    CellPos pos;
    REQUIRE(intp->locate(pts[i], 0., pos));
    ps.set_cell_id(i, pos.id);
  }
  ps.move(0, pts[3]);
  REQUIRE(ps.get_cell_id(0) == -1);
  REQUIRE(ps.insert_node_between(1, 2, true));
  REQUIRE(ps.insert_node_between(1, 3, false));
  REQUIRE(ps.N() == 6);
  const int id = ps.get_cell_id(4);
  REQUIRE(id >= 0);
  std::array<double, 4> lev;
  intp->region_point(intp->region_in(id, ps.x(4)), ps.x(4), lev);
  REQUIRE(min_level<3>(lev) >= 0.);
  REQUIRE(ps.get_cell_id(5) == -1);
}

TEST_CASE("A landing that ends the major step stores a region holding the end; one past the end is not taken", "[cells]"){
  SECTION("P2 triangles, periodic in x"){
    CaseDir c("cells_land_tri");
    write_stamped<Triangle>(c, 6, "P2", {true, false, false}, false);
    check_landings(std::make_shared<SimplexInterpol<Triangle>>(c.file("h5_params.dat")), 2, {true, false, false},
                   0.3, 0.3);
  }
  SECTION("P2 tets, periodic"){
    CaseDir c("cells_land_tet");
    write_stamped<Tet>(c, 3, "P2", {true, true, true}, false);
    check_landings(std::make_shared<SimplexInterpol<Tet>>(c.file("h5_params.dat")), 3, {true, true, true}, 0.3, 0.3);
  }
}

TEST_CASE("Into a narrower cell from within the band: crossed once, not back and forth", "[cells]"){
  CaseDir c("cells_graded_tri");
  write_stamped<Triangle>(c, 6, "P2", {false, false, false}, true, graded_x);
  check_into_narrower(std::make_shared<SimplexInterpol<Triangle>>(c.file("h5_params.dat")), 2, 0.3);
}

TEST_CASE("Past the cap on steps a major step: RK4's substeps from the start", "[cells]"){
  CaseDir c("cells_cap_tri");
  // 39 columns over x < 0.3, one beyond
  write_stamped<Triangle>(c, 40, "P2", {false, false, false}, false,
                          [](const Vector3d& x){
                            const double i = std::round(40.*x[0]);
                            return Vector3d(i < 40. ? 0.3*i/39. : 1., x[1], 0.); },
                          [](const Vector3d&, const int){ return Vector3d(1., 0.05, 0.); });
  check_cap(std::make_shared<SimplexInterpol<Triangle>>(c.file("h5_params.dat")), 0.3, 0.34);
}

TEST_CASE("An end below two walls at an acute angle: onto both and inside", "[cells]"){
  CaseDir c("cells_acute_tri");
  write_stamped<Triangle>(c, 4, "P2", {false, false, false}, true,
                          [](const Vector3d& x){ return Vector3d(x[0] - 1.5*x[1], x[1], 0.); });
  check_acute_corner(std::make_shared<SimplexInterpol<Triangle>>(c.file("h5_params.dat")), 0.3);
}

TEST_CASE("Split triangles: steps landing on predicted planes converge at fourth order", "[cells][split]"){
  CaseDir c("cells_split_order_tri");
  write_split<Triangle>(c, 6, false);
  auto intp = std::make_shared<SplitInterpol<Triangle>>(c.file("h5_params.dat"));
  intp->update(0.3);
  const std::vector<Vector3d> pts = central_points(*intp, 2, 0.3, 100, 0.35, 0.65);
  REQUIRE(pts.size() == 100);
  check_order<TransportElement::Tensor>(intp, pts, 0.3, 0.1, 0.005);
}

TEST_CASE("XDMF triangles: line elements and frames across facets converge at fourth order", "[cells]"){
  CaseDir c("cells_xdmf_order_tri");
  write_stamped<Triangle>(c, 6, "P1", {false, false, false}, true);
  std::ofstream(c.file("xdmf_params.dat"), std::ios::app) << "wall_p2=edge\n";
  auto intp = std::make_shared<XDMFInterpol<Triangle>>(c.file("xdmf_params.dat"));
  intp->update(0.3);
  const std::vector<Vector3d> pts = central_points(*intp, 2, 0.3, 100, 0.3, 0.7);
  REQUIRE(pts.size() == 100);
  check_order<TransportElement::Vector>(intp, pts, 0.3, 0.4, 0.04);
  check_order<TransportElement::Tensor>(intp, pts, 0.3, 0.4, 0.04);
}

TEST_CASE("onto_walls: onto two walls, nudged until inside; it fails where they do not meet or a level moves too far", "[cells]"){
  SECTION("tets, an acute edge"){
    CaseDir c("cells_onto_tet");
    write_stamped<Tet>(c, 2, "P2", {false, false, false}, true,
                       [](const Vector3d& x){ return Vector3d(x[0] - 1.5*x[1], x[1], x[2]); });
    check_onto_walls(std::make_shared<SimplexInterpol<Tet>>(c.file("h5_params.dat")), 3, 0.3);
  }
  SECTION("a mesh far from the origin"){
    CaseDir c("cells_onto_far_tri");
    write_stamped<Triangle>(c, 2, "P2", {false, false, false}, true,
                            [](const Vector3d& x){ return Vector3d(x[0] + 0.3*x[1] + 1e8, x[1] + 1e8, 0.); },
                            [](const Vector3d&, const int){ return Vector3d::Zero(); });
    check_onto_walls_far(std::make_shared<SimplexInterpol<Triangle>>(c.file("h5_params.dat")), 0.3);
  }
}

TEST_CASE("An end below the walls of a sliver corner: the projection fails, RK4's substeps from the start", "[cells]"){
  CaseDir c("cells_sliver_tri");
  // Drifting slowly down: the fallback starts again from the step's start
  write_stamped<Triangle>(c, 2, "P2", {false, false, false}, true,
                          [](const Vector3d& x){ return Vector3d(x[0] - 1e4*x[1], x[1], 0.); },
                          [](const Vector3d&, const int k){ return Vector3d(0., -1e-8*(1. + k), 0.); });
  // Beside the corner at (1, 0), 1e-4 rad wide
  check_projection_fallback(std::make_shared<SimplexInterpol<Triangle>>(c.file("h5_params.dat")),
                            {Vector3d(0.75, 1e-6, 0.)}, 0.3);
}

TEST_CASE("fall_time: the first downward root of l + b tau + a tau^2; none for a level that never reaches zero", "[cells]"){
  const double inf = cells::inf;
  REQUIRE(cells::fall_time(0.5, -2., 0.) == Approx(0.25));
  REQUIRE(cells::fall_time(0.5, 2., 0.) == inf);
  // (1 - tau)(1 - 2 tau): falls through at 1/2, rises through at 1
  REQUIRE(cells::fall_time(1., -3., 2.) == Approx(0.5));
  REQUIRE(cells::fall_time(1., 1., -1.) == Approx(0.5*(1. + std::sqrt(5.))));
  // Negative discriminant
  REQUIRE(cells::fall_time(1., -1., 1.) == inf);
  REQUIRE(cells::fall_time(0.35813077214585504, -1.4005753639509075, 1.369340128434507) == inf);
  // first_crossing reaches it when the level at the horizon rounds to zero (no FMA)
  struct NoWalls {
    using region_type = int;
    bool is_wall(const int, const int) const { return false; }
  };
  const std::array<double, 4> lev = {0.35813077214585504, 1., 1., 1.};
  const std::array<double, 4> d1 = {-1.4005753639509075, 0., 0., 0.};
  const std::array<double, 4> d2 = {2.*1.369340128434507, 0., 0., 0.};
  int k = -1;
  REQUIRE(cells::first_crossing<3>(NoWalls{}, 0, lev, d1, d2, 0.5114052144050252, -1, 1e-12, k) == inf);
  REQUIRE(k == -1);
}

#endif
