// Regions of the mesh loaders (regions.hpp) and the held evaluator
// (held_eval.hpp), on dolfin's unit meshes and on cases dolfin writes here.
//
// A region is a cell and the offset carrying an unwrapped point into that
// cell's periodic image; its levels are the cell's barycentrics at the
// evaluation point. What the stepping across facets relies on, checked here:
// the levels are contains' barycentrics (and, where the cell holds the wrapped
// point, locate's, bit for bit); across lands in the cell the walk finds for a
// point just beyond the facet, periodic partners and the boxes one cell thick
// (where a neighbour meets a cell across two facets) included; the held
// evaluator, which gathers a cell's nodes once, gives the loader's own
// evaluation bit for bit; and the velocity from the two sides of a facet
// agrees to round-off on every kind of facet the loaders have (P1, P2, the
// near-wall P2 rule, OpenFOAM's split, the frequency loaders), since the
// scheme assumes a continuous field and carries no jump term.
//
// The divergence-free split's regions are the sub-cells of a cell's
// barycentric split, its levels a sub-cell's barycentrics; the lattice's are
// the floor cells whose nodes are all fluid and the eight sub-cubes of the
// others, its levels the evaluation's own 1D weights. For both: the levels
// and their rates, across on every kind of plane (a sub-cell's internal planes
// and its macro facet; a lattice's faces between floor cells and sub-cubes,
// the sub-cubes' mid-planes, the walls before a solid node), the held
// evaluator bit for bit, and continuity.
#ifdef USE_DOLFIN
#include <catch2/catch.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <dolfin.h>

#include "MeshCore.hpp"
#include "OpenFoamInterpol.hpp"
#include "SimplexInterpol.hpp"
#include "SplitInterpol.hpp"
#include "StructuredInterpol.hpp"
#include "interpolators.hpp"
#include "XDMFInterpol.hpp"
#include "case_dir.hpp"
#include "dolfin_ref.hpp"
#include "held_eval.hpp"
#include "openfoam_load.hpp"
#include "region_cases.hpp"
#include "regions.hpp"

namespace {

// The cell and facet tables of a unit mesh, without fields or files
template<typename Cell>
struct BareMesh : public MeshCore<Cell> {
  std::shared_ptr<dolfin::Mesh> mesh;   // the dolfin cells refer to it
  std::vector<dolfin::Cell> dolfin_cells_;
  BareMesh(std::shared_ptr<dolfin::Mesh> m, const std::vector<bool>& periodic) : MeshCore<Cell>(""), mesh(m) {
    this->periodic = periodic;
    this->dim = m->geometry().dim();
    this->x_min = Vector3d::Zero();
    this->x_max = Vector3d(1., 1., this->dim == 3 ? 1. : 0.);
    for (dolfin::CellIterator c(*m); !c.end(); ++c){
      dolfin_cells_.push_back(*c);
      this->cells_.push_back(Cell(*c));
    }
    build_facet_neighbours(this->facet_neigh_, m, dolfin_cells_, nullptr, this->periodic,
                           this->x_min, this->x_max, this->dim, this->periodic_tol);
    this->set_period();
  }
  const Cell& cell(const int id) const { return this->cells_[std::size_t(id)]; }
  std::int32_t facet(const int id, const int k) const { return this->facet_neigh_[std::size_t(id)*Cell::n_verts + k]; }
  const std::vector<Cell>& cells() const { return this->cells_; }
  const std::vector<std::int32_t>& facets() const { return this->facet_neigh_; }
  Vector3d wrap(const Vector3d& x) const { return this->_modx(x); }
  Vector3d vertex(const int id, const int k) const {
    for (dolfin::VertexIterator v(dolfin_cells_[std::size_t(id)]); !v.end(); ++v)
      if (int(v.pos()) == k){
        Vector3d x = Vector3d::Zero();
        for (Uint d = 0; d < this->dim; ++d) x[d] = v->x(d);
        return x;
      }
    return Vector3d::Zero();
  }
  // A point of the cell from barycentric weights
  Vector3d point(const int id, const std::array<double, 4>& w) const {
    Vector3d x = Vector3d::Zero();
    for (int k = 0; k < Cell::n_verts; ++k) x += w[k]*vertex(id, k);
    return x;
  }
  bool locate(const Vector3d& x, const double, CellPos& pos){
    return walk_to_cell(this->cells_, this->facet_neigh_, this->_modx(x), pos);
  }
  using MeshCore<Cell>::locate;
  void update(const double) {}
  void evaluate(const Vector3d&, const double, const CellPos&, PointValues&) {}
  double get_t_min() { return 0.; }
  double get_t_max() { return 1.; }
};

template<typename Cell>
void check_levels(){
  constexpr int nv = Cell::n_verts;
  BareMesh<Cell> m(unit_mesh<Cell>(4), {true, true, false});
  std::mt19937 rng(3);
  std::uniform_int_distribution<int> any_cell(0, int(m.cells().size()) - 1), image(-2, 2);
  std::uniform_real_distribution<double> uni(-1., 1.);
  for (int s = 0; s < 2000; ++s){
    const int id = any_cell(rng);
    // Unwrapped: a few periods off along the periodic axes
    const Vector3d x = m.point(id, random_weights<Cell>(rng)) + Vector3d(image(rng), image(rng), 0.);
    Region R;
    std::array<double, 4> lev, b;
    REQUIRE(m.region_of(id, x, 0., R, lev));
    REQUIRE(R.id == id);
    // Where the cell holds the wrapped point: contains' barycentrics there, as locate's
    REQUIRE(m.cell(id).contains(m.wrap(x), b));
    for (int k = 0; k < nv; ++k) REQUIRE(lev[k] == b[k]);
    std::array<double, 4> l2;
    REQUIRE(m.region_point(R, x, l2) == m.wrap(x));
    for (int k = 0; k < nv; ++k) REQUIRE(l2[k] == b[k]);
    REQUIRE((x + R.offset - m.wrap(x)).norm() < 1e-14);
    // Beyond the cell: x + offset, extrapolated
    const Vector3d y = x + 0.5*Vector3d(uni(rng), uni(rng), dim_of<Cell> == 3 ? uni(rng) : 0.);
    if (!m.cell(id).contains(m.wrap(y), b)){
      REQUIRE(m.region_point(R, y, l2) == y + R.offset);
      m.cell(id).contains(y + R.offset, b);
      for (int k = 0; k < nv; ++k) REQUIRE(l2[k] == b[k]);
    }
    // Rates along v from the barycentric gradients
    const Vector3d v(uni(rng), uni(rng), dim_of<Cell> == 3 ? uni(rng) : 0.);
    std::array<double, 4> rate;
    m.level_rates(R, v, rate);
    for (int k = 0; k < nv; ++k) REQUIRE(rate[k] == m.cell(id).bary_grad(k).dot(v));
    // A stale cell beside it: found by the walk
    const std::int32_t nb = m.facet(id, s % nv);
    if (nb >= 0){
      Region S;
      REQUIRE(m.region_of(nb, x, 0., S, l2));
      REQUIRE(S.id == id);
      REQUIRE(S.offset == R.offset);
    }
  }
}

// A point on the periodic face x = 1, stored in a cell on that side: the cell
// is kept with the offset moved by a period, though the wrapped point is at x = 0
template<typename Cell>
void check_face_point(){
  constexpr int nv = Cell::n_verts;
  BareMesh<Cell> m(unit_mesh<Cell>(4), {true, false, false});
  std::size_t kept = 0;
  for (int id = 0; id < int(m.cells().size()); ++id)
    for (int k = 0; k < nv; ++k){
      if (m.facet(id, k) > facet_wall) continue;
      Vector3d c = Vector3d::Zero();
      for (int j = 0; j < nv; ++j) if (j != k) c += m.vertex(id, j)/double(nv - 1);
      if (c[0] != 1.) continue;
      Region R;
      std::array<double, 4> lev;
      REQUIRE(m.wrap(c)[0] == 0.);
      REQUIRE(m.region_of(id, c, 1e-12, R, lev));
      REQUIRE(R.id == id);
      REQUIRE(R.offset.norm() == 0.);
      REQUIRE(min_level<nv>(lev) >= -1e-12);
      ++kept;
    }
  REQUIRE(kept > 0);
}

// across against the walk, from every facet of every cell: the neighbour
// holds a point just beyond the facet, the entered facet names the cell
// left, and across from it comes back
template<typename Cell>
void check_across(const std::size_t n, const std::vector<bool>& periodic, const bool expect_twice){
  constexpr int nv = Cell::n_verts;
  BareMesh<Cell> m(unit_mesh<Cell>(n), periodic);
  std::size_t plain = 0, partner = 0, walls = 0, twice = 0;
  for (int id = 0; id < int(m.cells().size()); ++id)
    for (int k = 0; k < nv; ++k){
      Vector3d c = Vector3d::Zero();
      for (int j = 0; j < nv; ++j) if (j != k) c += m.vertex(id, j)/double(nv - 1);
      const Vector3d g = m.cell(id).bary_grad(k);
      const double h = 1e-6/double(n);
      const Vector3d out = c - h*g.normalized(), in = c + h*g.normalized();
      const Region R{id, Vector3d::Zero()};
      Region next;
      const int e = m.across(R, k, out, next);
      const std::int32_t a = m.facet(id, k);
      if (a == facet_wall){
        REQUIRE(e == -1);
        ++walls;
        continue;
      }
      REQUIRE(e >= 0);
      (a >= 0 ? plain : partner)++;
      // The walk's cell
      CellPos pos;
      pos.id = id;
      REQUIRE(walk_to_cell(m.cells(), m.facets(), m.wrap(out), pos));
      REQUIRE(pos.id == next.id);
      REQUIRE((out + next.offset - m.wrap(out)).norm() < 1e-14);
      // Held by it, the entered facet nearest
      std::array<double, 4> lev;
      m.region_point(next, out, lev);
      REQUIRE(min_level<nv>(lev) >= 0.);
      REQUIRE(m.facet(next.id, e) == (a >= 0 ? id : facet_periodic(id)));
      for (int j = 0; j < nv; ++j) REQUIRE(lev[e] <= lev[j]);
      int matches = 0;
      for (int j = 0; j < nv; ++j) matches += m.facet(next.id, j) == m.facet(next.id, e);
      if (matches > 1) ++twice;
      // And back
      Region back;
      REQUIRE(m.across(next, e, in, back) == k);
      REQUIRE(back.id == id);
      REQUIRE(back.offset == R.offset);
    }
  REQUIRE(plain > 0);
  if (!periodic[1]) REQUIRE(walls > 0);
  if (periodic[0]) REQUIRE(partner > 0);
  if (expect_twice) REQUIRE(twice > 0);
}

bool same_region(const Region& a, const Region& b){ return a.id == b.id; }
bool same_region(const SubRegion& a, const SubRegion& b){ return a.id == b.id && a.sub == b.sub; }

// The held evaluator against the loader's evaluation in the located cell, at
// the point and at a second point of the same region from the same gather
template<typename I>
std::size_t check_held_is_located(I& intp, const std::vector<Vector3d>& pts, const std::vector<double>& times){
  std::mt19937 rng(5);
  std::uniform_real_distribution<double> uni(-1., 1.);
  std::size_t n = 0;
  for (const int order : {1, 2}){
    intp.set_int_order(order);
    for (const double t : times){
      intp.update(t);
      for (const Vector3d& x : pts){
        CellPos pos;
        if (!intp.locate(x, t, pos)) continue;
        typename I::region_type R;
        typename I::levels_type lev;
        REQUIRE(intp.region_of(pos.id, x, 0., R, lev));
        REQUIRE(R.id == pos.id);
        PointValues b(1.);
        HeldEval<I> ev(intp, b);
        ev.bind(R);
        const Vector3d step = 1e-3*Vector3d(uni(rng), uni(rng), uni(rng));
        for (const Vector3d& y : {x, Vector3d(x + step)}){
          CellPos p2;
          p2.id = pos.id;
          if (!intp.locate(y, t, p2) || p2.id != pos.id) continue;
          typename I::region_type Ry;
          REQUIRE(intp.region_of(pos.id, y, 0., Ry, lev));
          if (!same_region(Ry, R)) continue;
          PointValues a(1.);
          intp.evaluate_motion(y, t, p2, a);
          REQUIRE(ev(y, t));
          REQUIRE(same_bits(a.U, b.U));
          REQUIRE(same_bits(a.A, b.A));
          if (order == 2){
            REQUIRE(same_bits(a.gradU, b.gradU));
            REQUIRE(same_bits(a.gradA, b.gradA));
          }
          ++n;
        }
      }
    }
  }
  return n;
}

// The largest jump of u across a facet over the facets of the located cells,
// relative to the largest |u|; counts the facets and the periodic ones
template<typename I>
double facet_jump(I& intp, const std::vector<Vector3d>& pts, const double t,
                  std::size_t& n_facets, std::size_t& n_periodic){
  constexpr int nv = I::n_levels;
  intp.set_int_order(1);
  intp.update(t);
  double jump = 0., umax = 0.;
  n_facets = n_periodic = 0;
  for (const Vector3d& x : pts){
    CellPos pos;
    if (!intp.locate(x, t, pos)) continue;
    typename I::region_type R;
    typename I::levels_type lev;
    REQUIRE(intp.region_of(pos.id, x, 0., R, lev));
    // The levels' gradients, as rates along the axes
    std::array<Vector3d, nv> g;
    for (int d = 0; d < 3; ++d){
      typename I::levels_type r;
      intp.level_rates(R, Vector3d::Unit(d), r);
      for (int k = 0; k < nv; ++k) g[k][d] = r[k];
    }
    for (int k = 0; k < nv; ++k){
      // Onto facet k along its normal
      const Vector3d xe = intp.region_point(R, x, lev);
      const Vector3d xf = xe - R.offset - lev[k]/g[k].squaredNorm()*g[k];
      typename I::levels_type l2;
      intp.region_point(R, xf, l2);
      bool clear = true;
      for (int j = 0; j < nv; ++j) if (j != k && l2[j] < 1e-3) clear = false;
      if (!clear) continue;
      typename I::region_type R2;
      if (intp.across(R, k, xf, R2) < 0) continue;
      PointValues p1(1.), p2(1.);
      HeldEval<I> e1(intp, p1), e2(intp, p2);
      e1.bind(R);
      e2.bind(R2);
      e1.at(xf, t);
      e2.at(xf, t);
      jump = std::max(jump, (p1.U - p2.U).norm());
      umax = std::max(umax, p1.U.norm());
      ++n_facets;
      if (R2.offset != R.offset) ++n_periodic;
    }
  }
  REQUIRE(umax > 0.);
  return jump/umax;
}

template<typename I>
void check_loader(I& intp, const int D, const std::vector<bool>& per, const std::vector<double>& times,
                  const double t_jump){
  const std::vector<Vector3d> pts = box_points(intp, D, per, 300);
  REQUIRE(check_held_is_located(intp, pts, times) > pts.size());
  std::size_t n_facets = 0, n_periodic = 0;
  const double j = facet_jump(intp, pts, t_jump, n_facets, n_periodic);
  INFO("relative jump " << j << " over " << n_facets << " facets, " << n_periodic << " periodic");
  REQUIRE(n_facets > pts.size());
  if (per[0]) REQUIRE(n_periodic > 0);
  REQUIRE(j < 1e-13);
}

template<typename Cell>
void check_stamped(const std::string& tag, const std::size_t n, const std::string& u_el,
                   const std::vector<bool>& per, const bool rest, const std::string& extra){
  constexpr int D = dim_of<Cell>;
  CaseDir c(tag);
  write_stamped<Cell>(c, n, u_el, per, rest);
  if (!extra.empty()) std::ofstream(c.file("h5_params.dat"), std::ios::app) << extra;
  {
    SimplexInterpol<Cell> intp(c.file("h5_params.dat"));
    check_loader(intp, D, per, {0., 0.3, 1.}, 0.3);
  }
  if (u_el == "P1"){
    if (!extra.empty()) std::ofstream(c.file("xdmf_params.dat"), std::ios::app) << extra;
    XDMFInterpol<Cell> intp(c.file("xdmf_params.dat"));
    check_loader(intp, D, per, {0., 0.3, 1.}, 0.3);
  }
}

template<typename Cell, typename I>
void check_freq(const std::string& tag, const std::size_t n){
  constexpr int D = dim_of<Cell>;
  CaseDir c(tag);
  write_freq<Cell>(c, n);
  I intp(c.params());
  check_loader(intp, D, {true, true, D == 3}, {0., 0.25, 0.7}, 0.25);
}

// The split's levels at located points: the sub-cell of the smallest
// barycentric, its barycentrics mu as the evaluation computes them, and their
// rates along v the gradients' (the levels are affine in x)
template<typename Cell>
void check_split_levels(SplitInterpol<Cell>& intp, const std::vector<Vector3d>& pts){
  constexpr int nv = Cell::n_verts;
  std::mt19937 rng(9);
  std::uniform_real_distribution<double> uni(-1., 1.);
  std::size_t n = 0;
  for (const Vector3d& x : pts){
    CellPos pos;
    if (!intp.locate(x, 0., pos)) continue;
    SubRegion R;
    std::array<double, 4> lev, l2, mu;
    REQUIRE(intp.region_of(pos.id, x, 0., R, lev));
    REQUIRE(R.id == pos.id);
    REQUIRE(R.sub == split_eval::sub_cell<nv>(pos.bary, mu.data()));
    for (int k = 0; k < nv; ++k) REQUIRE(lev[k] == mu[k]);
    REQUIRE(min_level<nv>(lev) >= 0.);
    const SubRegion S = intp.region_in(pos.id, x);
    REQUIRE((S.id == R.id && S.sub == R.sub && S.offset == R.offset));
    REQUIRE(intp.region_point(R, x, l2) == x + R.offset);
    for (int k = 0; k < nv; ++k) REQUIRE(l2[k] == mu[k]);
    const Vector3d v(uni(rng), uni(rng), dim_of<Cell> == 3 ? uni(rng) : 0.);
    std::array<double, 4> rate, l3;
    intp.level_rates(R, v, rate);
    intp.region_point(R, x + 1e-3*v, l3);
    for (int k = 0; k < nv; ++k){
      const double g = intp.level_grad(R, k).dot(v);
      REQUIRE(std::abs(rate[k] - g) <= 1e-12*(1. + std::abs(g)));
      REQUIRE(std::abs((l3[k] - lev[k])/1e-3 - rate[k]) <= 1e-9*(1. + std::abs(g)));
    }
    ++n;
  }
  REQUIRE(n > pts.size()/2);
}

// across from every plane of the located sub-cells, against locate: an
// internal plane leads to the sub-cell of the same cell beyond it, the macro
// facet to the neighbour's sub-cell on it, periodic partners included; a point
// just beyond is held there with the entered level the lowest, locate puts it
// in that cell and sub-cell, and across from the entered plane comes back
template<typename Cell>
void check_split_across(SplitInterpol<Cell>& intp, const std::vector<Vector3d>& pts, const bool periodic){
  constexpr int nv = Cell::n_verts;
  std::size_t internal = 0, macro = 0, partner = 0, walls = 0;
  for (const Vector3d& x : pts){
    CellPos pos;
    if (!intp.locate(x, 0., pos)) continue;
    SubRegion R;
    std::array<double, 4> lev;
    REQUIRE(intp.region_of(pos.id, x, 0., R, lev));
    for (int k = 0; k < nv; ++k){
      const Vector3d g = intp.level_grad(R, k);
      const Vector3d xe = intp.region_point(R, x, lev);
      const Vector3d xf = xe - R.offset - lev[k]/g.squaredNorm()*g;
      std::array<double, 4> l2;
      intp.region_point(R, xf, l2);
      bool clear = true;
      for (int j = 0; j < nv; ++j) if (j != k && l2[j] < 1e-2) clear = false;
      if (!clear) continue;
      const double h = 1e-7/g.norm();
      const Vector3d out = xf - h*g.normalized(), in = xf + h*g.normalized();
      SubRegion next;
      const int e = intp.across(R, k, out, next);
      if (intp.is_wall(R, k)){
        REQUIRE(e == -1);
        ++walls;
        continue;
      }
      REQUIRE(e >= 0);
      if (k < nv - 1){
        REQUIRE(next.id == R.id);
        REQUIRE(next.sub == split_eval::sub_beyond(R.sub, k));
        REQUIRE(e < nv - 1);
        ++internal;
      }
      else {
        REQUIRE(e == nv - 1);
        ++macro;
        if (next.offset != R.offset) ++partner;
      }
      intp.region_point(next, out, l2);
      REQUIRE(min_level<nv>(l2) >= 0.);
      for (int j = 0; j < nv; ++j) REQUIRE(l2[e] <= l2[j]);
      CellPos q;
      q.id = R.id;
      REQUIRE(intp.locate(out, 0., q));
      std::array<double, 4> mu;
      REQUIRE(q.id == next.id);
      REQUIRE(split_eval::sub_cell<nv>(q.bary, mu.data()) == next.sub);
      SubRegion back;
      REQUIRE(intp.across(next, e, in, back) == k);
      REQUIRE((back.id == R.id && back.sub == R.sub && back.offset == R.offset));
    }
  }
  REQUIRE(internal > pts.size());
  REQUIRE(macro > pts.size()/4);
  if (periodic) REQUIRE(partner > 0);
  else REQUIRE(walls > 0);
}

template<typename Cell>
void check_split(const std::string& tag, const std::size_t n, const bool periodic){
  constexpr int D = dim_of<Cell>;
  CaseDir c(tag);
  write_split<Cell>(c, n, periodic);
  SplitInterpol<Cell> intp(c.file("h5_params.dat"));
  const std::vector<bool> per = {periodic, false, false};
  intp.update(0.);
  const std::vector<Vector3d> pts = box_points(intp, D, per, 300);
  check_split_levels(intp, pts);
  check_split_across(intp, pts, periodic);
  check_loader(intp, D, per, {0., 0.3, 1.}, 0.3);
}

// The lattice's regions at fluid points: the floor cell where its nodes are
// all fluid, else the sub-cube of the point's halves; levels the evaluation's
// own weights (lower face, upper face, along each axis), affine in x
void check_lattice_levels(StructuredInterpol& intp, const int n, const std::vector<Vector3d>& pts,
                          std::size_t& bulk, std::size_t& near){
  std::mt19937 rng(9);
  std::uniform_real_distribution<double> uni(-1., 1.);
  bulk = near = 0;
  for (const Vector3d& x : pts){
    LatticeRegion R;
    std::array<double, 6> lev, l2, l3, rate;
    REQUIRE(intp.region_of(-1, x, 0., R, lev));
    bool fluid = true;
    for (int c = 0; c < 8; ++c)
      fluid = fluid && !lattice_solid(n, int(imodulo(R.f[0] + (c & 1), n)), int(imodulo(R.f[1] + (c >> 1 & 1), n)),
                                      int(imodulo(R.f[2] + (c >> 2 & 1), n)));
    int sub = 0;
    for (int a = 0; a < 3; ++a){
      REQUIRE(R.f[a] == int(std::floor(x[a])));
      if (x[a] - R.f[a] >= 0.5) sub |= 1 << a;
    }
    REQUIRE(R.sub == (fluid ? -1 : sub));
    (fluid ? bulk : near)++;
    REQUIRE(intp.region_point(R, x, l2) == x);
    for (int k = 0; k < 6; ++k){
      REQUIRE(l2[k] == lev[k]);
      REQUIRE(lev[k] >= 0.);
      REQUIRE(lev[k] <= 1.);
    }
    const LatticeRegion S = intp.region_in(R.id, x);
    REQUIRE((S.f == R.f && S.sub == R.sub && S.id == R.id));
    const Vector3d v(uni(rng), uni(rng), uni(rng));
    intp.level_rates(R, v, rate);
    intp.region_point(R, x + 1e-3*v, l3);
    for (int k = 0; k < 6; ++k){
      REQUIRE(std::abs(rate[k] - intp.level_grad(R, k).dot(v)) <= 1e-15*std::abs(rate[k]));
      REQUIRE(std::abs((l3[k] - lev[k])/1e-3 - rate[k]) <= 1e-9);
    }
  }
}

// across from every face of the regions at the points, against the region a
// point just beyond is in: a floor cell's faces and a sub-cube's outer ones
// lead to the next floor cell or its sub-cube on the face, a sub-cube's
// mid-planes to the next sub-cube of the cell, or to a wall where that one's
// node is solid (locate refuses the point); the entered level the lowest, and
// across from the entered face comes back
void check_lattice_across(StructuredInterpol& intp, const std::vector<Vector3d>& pts,
                          std::array<std::size_t, 6>& kinds){
  kinds.fill(0);   // bulk-bulk, bulk-near, near-bulk, near-near, mid-plane, wall
  for (const Vector3d& x : pts){
    LatticeRegion R;
    std::array<double, 6> lev, l2;
    REQUIRE(intp.region_of(-1, x, 0., R, lev));
    for (int k = 0; k < 6; ++k){
      const Vector3d g = intp.level_grad(R, k);
      const Vector3d xf = x - lev[k]/g.squaredNorm()*g;
      intp.region_point(R, xf, l2);
      bool clear = true;
      for (int j = 0; j < 6; ++j) if (j != k && l2[j] < 1e-3) clear = false;
      if (!clear) continue;
      const Vector3d out = xf - 1e-7*g.normalized(), in = xf + 1e-7*g.normalized();
      LatticeRegion next;
      const int e = intp.across(R, k, out, next);
      CellPos q;
      if (e < 0){
        REQUIRE(intp.is_wall(R, k));
        REQUIRE(!intp.locate(out, 0., q));
        REQUIRE(intp.locate(in, 0., q));
        ++kinds[5];
        continue;
      }
      REQUIRE(!intp.is_wall(R, k));
      REQUIRE(e == (k ^ 1));
      ++kinds[next.f == R.f ? 4 : 2*(R.sub >= 0) + (next.sub >= 0)];
      intp.region_point(next, out, l2);
      for (int j = 0; j < 6; ++j){
        REQUIRE(l2[j] >= 0.);
        REQUIRE(l2[e] <= l2[j]);
      }
      LatticeRegion S;
      REQUIRE(intp.region_of(-1, out, 0., S, l2));
      REQUIRE((S.f == next.f && S.sub == next.sub && S.id == next.id));
      LatticeRegion back;
      REQUIRE(intp.across(next, e, in, back) == k);
      REQUIRE((back.f == R.f && back.sub == R.sub && back.id == R.id));
    }
  }
}

// The held evaluator bit for bit the lattice's evaluation in the same region,
// at the point and at a second point; u from the two sides of every face
// agrees to round-off; J is the gradient of u in the region
void check_lattice_fields(StructuredInterpol& intp, const std::vector<Vector3d>& pts, const std::vector<double>& times){
  std::mt19937 rng(5);
  std::uniform_real_distribution<double> uni(-1., 1.);
  std::size_t same = 0, faces = 0, grads = 0;
  double jump = 0., umax = 0., gerr = 0., gmax = 0.;
  for (const double t : times){
    intp.update(t);
    for (const Vector3d& x : pts){
      LatticeRegion R;
      std::array<double, 6> lev, l2;
      REQUIRE(intp.region_of(-1, x, 0., R, lev));
      PointValues b(1.);
      HeldEval<StructuredInterpol> ev(intp, b);
      ev.bind(R);
      const Vector3d step = 0.1*Vector3d(uni(rng), uni(rng), uni(rng));
      for (const Vector3d& y : {x, Vector3d(x + step)}){
        LatticeRegion S;
        CellPos pos;
        if (!intp.region_of(-1, y, 0., S, l2) || S.f != R.f || S.sub != R.sub || !intp.locate(y, t, pos)) continue;
        PointValues a(1.);
        intp.evaluate(y, t, pos, a);
        REQUIRE(ev(y, t));
        REQUIRE(same_bits(a.U, b.U));
        REQUIRE(same_bits(a.A, b.A));
        REQUIRE(same_bits(a.gradU, b.gradU));
        PointValues c(1.);
        HeldEval<StructuredInterpol> ev2(intp, c);
        ev2.bind(R);
        REQUIRE(ev2.velocity(y, t));
        REQUIRE(same_bits(a.U, c.U));
        ++same;
      }
      // J against central differences of u in the region
      if (*std::min_element(lev.begin(), lev.end()) > 1e-3){
        const double h = 1e-5;
        ev.at(x, t);
        const Matrix3d J = b.gradU;
        for (int d = 0; d < 3; ++d){
          ev.at(x + h*Vector3d::Unit(d), t);
          const Vector3d up = b.U;
          ev.at(x - h*Vector3d::Unit(d), t);
          gerr = std::max(gerr, (J.col(d) - (up - b.U)/(2.*h)).norm());
        }
        gmax = std::max(gmax, J.norm());
        ++grads;
      }
      for (int k = 0; k < 6; ++k){
        const Vector3d g = intp.level_grad(R, k);
        const Vector3d xf = x - lev[k]/g.squaredNorm()*g;
        intp.region_point(R, xf, l2);
        bool clear = true;
        for (int j = 0; j < 6; ++j) if (j != k && l2[j] < 1e-3) clear = false;
        LatticeRegion next;
        if (!clear || intp.across(R, k, xf, next) < 0) continue;
        PointValues p1(1.), p2(1.);
        HeldEval<StructuredInterpol> e1(intp, p1), e2(intp, p2);
        e1.bind(R);
        e2.bind(next);
        e1.at(xf, t);
        e2.at(xf, t);
        jump = std::max(jump, (p1.U - p2.U).norm());
        umax = std::max(umax, p1.U.norm());
        ++faces;
      }
    }
  }
  INFO("bit for bit at " << same << " points; relative jump " << jump/umax << " over " << faces
       << " faces; J's error " << gerr/gmax << " at " << grads << " points");
  REQUIRE(same > pts.size());
  REQUIRE(faces > pts.size());
  REQUIRE(jump <= 1e-13*umax);
  REQUIRE(grads > pts.size()/2);
  REQUIRE(gerr <= 1e-6*gmax);
}

}  // namespace

TEST_CASE("Levels are contains' barycentrics at the evaluation point", "[regions]"){
  SECTION("triangle"){ check_levels<Triangle>(); }
  SECTION("tet"){ check_levels<Tet>(); }
}

TEST_CASE("A point on a periodic face keeps its stored cell", "[regions]"){
  SECTION("triangle"){ check_face_point<Triangle>(); }
  SECTION("tet"){ check_face_point<Tet>(); }
}

TEST_CASE("across lands where the walk finds a point just beyond the facet", "[regions]"){
  SECTION("triangle, walls"){ check_across<Triangle>(4, {false, false, false}, false); }
  SECTION("tet, walls"){ check_across<Tet>(3, {false, false, false}, false); }
  SECTION("triangle, periodic in x"){ check_across<Triangle>(4, {true, false, false}, false); }
  SECTION("tet, periodic in x and y"){ check_across<Tet>(3, {true, true, false}, false); }
  // One cell thick: the triangles meet a neighbour across two facets, told
  // apart by the levels (dolfin's six tets of a cube never do)
  SECTION("triangle, one cell, periodic in x and y"){ check_across<Triangle>(1, {true, true, false}, true); }
  SECTION("tet, one cell, periodic in x, y and z"){ check_across<Tet>(1, {true, true, true}, false); }
}

TEST_CASE("The held evaluator is the loader's evaluation, and u is continuous across facets", "[regions]"){
  SECTION("P1 triangles, periodic"){ check_stamped<Triangle>("reg_p1_tri", 6, "P1", {true, true, false}, false, ""); }
  SECTION("P2 triangles, periodic"){ check_stamped<Triangle>("reg_p2_tri", 6, "P2", {true, true, false}, false, ""); }
  SECTION("P1 tets, periodic"){ check_stamped<Tet>("reg_p1_tet", 3, "P1", {true, true, true}, false, ""); }
  SECTION("P2 tets, periodic"){ check_stamped<Tet>("reg_p2_tet", 3, "P2", {true, false, false}, false, ""); }
  SECTION("near-wall P2 triangles"){ check_stamped<Triangle>("reg_wall_tri", 6, "P1", {false, false, false}, true, "wall_p2=edge\n"); }
  SECTION("near-wall P2 tets"){ check_stamped<Tet>("reg_wall_tet", 3, "P1", {false, false, false}, true, "wall_p2=edge\n"); }
  SECTION("frequency triangles"){ check_freq<Triangle, TriangleFreqInterpol>("reg_freq_tri", 6); }
  SECTION("frequency tets"){ check_freq<Tet, TetFreqInterpol>("reg_freq_tet", 3); }
}

TEST_CASE("The divergence-free split: sub-cells, their planes, the held evaluator and continuity", "[regions][split]"){
  SECTION("triangles, periodic in x"){ check_split<Triangle>("reg_split_tri", 6, true); }
  SECTION("triangles, walls"){ check_split<Triangle>("reg_split_tri_w", 6, false); }
  SECTION("tets, periodic in x"){ check_split<Tet>("reg_split_tet", 3, true); }
  SECTION("tets, walls"){ check_split<Tet>("reg_split_tet_w", 3, false); }
}

TEST_CASE("The lattice: floor cells and sub-cubes, their faces, the held evaluator and continuity", "[regions][lattice]"){
  const int n = 12;
  CaseDir c("reg_lattice");
  write_lattice(c, n);
  StructuredInterpol intp(c.file("felbm_params.dat"));
  intp.update(0.);
  const std::vector<Vector3d> pts = lattice_points(intp, n, 600);
  std::size_t bulk = 0, near = 0;
  check_lattice_levels(intp, n, pts, bulk, near);
  INFO(bulk << " points in floor cells, " << near << " in sub-cubes");
  REQUIRE(bulk > 100);
  REQUIRE(near > 100);
  std::array<std::size_t, 6> kinds;
  check_lattice_across(intp, pts, kinds);
  INFO("faces: bulk-bulk " << kinds[0] << ", bulk-near " << kinds[1] << ", near-bulk " << kinds[2]
       << ", near-near " << kinds[3] << ", mid-planes " << kinds[4] << ", walls " << kinds[5]);
  for (const std::size_t k : kinds) REQUIRE(k > 20);
  check_lattice_fields(intp, pts, {0., 0.3, 1.});
}

TEST_CASE("OpenFOAM's split: the held evaluator and continuity", "[regions][openfoam]"){
  if (!openfoam_load::available()){
    WARN("no OpenFOAM reader in this build");
    return;
  }
  const std::string data = std::string(PARTRAC_SOURCE_DIR) + "/data_example/";
  SECTION("tets, cyclic in two directions"){
    OpenFoamInterpol<Tet> intp(data + "openfoam_channel3d/partrac_params.dat");
    const double t = intp.get_t_min();
    check_loader(intp, 3, {true, true, false}, {t}, t);
  }
  SECTION("triangles"){
    OpenFoamInterpol<Triangle> intp(data + "openfoam_cavity/partrac_params.dat");
    const double t0 = intp.get_t_min(), t1 = intp.get_t_max();
    check_loader(intp, 2, {false, false, false}, {t0, 0.5*(t0 + t1), t1}, 0.5*(t0 + t1));
  }
}

#endif
