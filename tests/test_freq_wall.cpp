// wall_upgrade in the frequency loaders (mode=trianglefreq and tetfreq), P1
// and P2, on dolfin's unit square and cube with every side a wall.
//
// A wall vertex is at rest only if every component is at rest there, within
// one tolerance for the sum, so the near-wall rule raises the summed field the
// way the stamped loader raises a stamp. That gives the references here: the
// stamped loader (mode=triangle and tet) reading the same files as two
// stamps, or reading one field that is the components' sum. Lines of
// `omega phi a file` chosen so that at t = 0.4 the weights a cos(omega t + phi)
// are the stamps' blend, 0.6 and 0.4, and their rates are the stamps' rate, -1
// and 1, make the two loaders' velocity, acceleration and both gradients the
// same field.
#ifdef USE_DOLFIN
#include <catch2/catch.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "SimplexInterpol.hpp"
#include "interpolators.hpp"
#include "case_dir.hpp"
#include "region_cases.hpp"

namespace {

constexpr double T0 = 0.4;
const std::vector<bool> walls = {false, false, false};

// Zero on every side of the unit box
double bubble(const Vector3d& x, const int D){
  double b = 1.;
  for (int d = 0; d < D; ++d) b *= x[d]*(1. - x[d]);
  return b;
}

// Two fields at rest on every wall
template<int D>
Vector3d rest_a(const Vector3d& x){
  return bubble(x, D)*Vector3d(4.*(1. + x[1]), 1. + x[0], D == 3 ? 1. + x[0] + x[1] : 0.);
}

template<int D>
Vector3d rest_b(const Vector3d& x){
  return bubble(x, D)*Vector3d(1. - x[0]*x[1], 2. + std::sin(3.*x[0]), D == 3 ? 0.5 + x[2] : 0.);
}

// Moving on the side x = 0 alone
template<int D>
Vector3d moving(const Vector3d& x){
  double b = 1. - x[0];
  for (int d = 1; d < D; ++d) b *= x[d]*(1. - x[d]);
  return b*Vector3d(1., 0.5, D == 3 ? 0.3 : 0.);
}

using Comp = Vector3d (*)(const Vector3d&);

// The lines whose weights at T0 are the blend (0.6, 0.4) and whose rates are
// (-1, 1); swapped: up_1's line first
std::string blend_lines(const bool swapped = false){
  std::ostringstream s;
  s.precision(17);
  const double c[2] = {1. - T0, T0}, r[2] = {-1., 1.};
  for (int i = 0; i < 2; ++i){
    const int k = swapped ? 1 - i : i;
    // a cos(t + phi) = c, a sin(t + phi) = -r
    const double a = std::hypot(c[k], r[k]), phi = std::atan2(-r[k], c[k]) - T0;
    s << "1 " << phi << " " << a << " up_" << k << ".h5\n";
  }
  return s.str();
}

// Both components a stamp, the stamped loader's parameters (h5_params.dat)
// and the frequency loader's (freq_params.dat), with wall_upgrade
template<typename Cell>
void write_components(const CaseDir& c, const std::size_t n, const std::string& u_el, const Comp f0,
                      const Comp f1, const std::string& lines){
  write_stamped<Cell>(c, n, u_el, walls, true, {}, [&](const Vector3d& x, const int k){
    return k == 0 ? f0(x) : f1(x);
  });
  std::ofstream(c.file("h5_params.dat"), std::ios::app) << "wall_upgrade=true\n";
  write_freq_of_stamps(c, u_el, walls, lines, "wall_upgrade=true\n");
}

// One field, w0 f0 + w1 f1, at both stamps
template<typename Cell>
void write_sum(const CaseDir& c, const std::size_t n, const std::string& u_el, const Comp f0, const Comp f1,
               const double w0, const double w1){
  write_stamped<Cell>(c, n, u_el, walls, true, {}, [&](const Vector3d& x, int){
    return Vector3d(w0*f0(x) + w1*f1(x));
  });
  std::ofstream(c.file("h5_params.dat"), std::ios::app) << "wall_upgrade=true\n";
}

// evaluate_motion, or evaluate, with the gradients at points of the box
template<typename I>
std::vector<PointValues> sample(I& intp, const int D, const double t, const bool motion = true){
  intp.set_int_order(2);
  intp.update(t);
  std::vector<PointValues> out;
  for (const Vector3d& x : box_points(intp, D, walls, 2000)){
    CellPos pos;
    REQUIRE(intp.locate(x, t, pos));
    PointValues v(1.);
    if (motion) intp.evaluate_motion(x, t, pos, v);
    else intp.evaluate(x, t, pos, v);
    out.push_back(v);
  }
  return out;
}

// The largest difference of a field over the points, relative to its largest value
template<typename F>
double gap(const std::vector<PointValues>& a, const std::vector<PointValues>& b, const F& f){
  double d = 0., m = 0.;
  for (std::size_t i = 0; i < a.size(); ++i){
    d = std::max(d, (f(a[i]) - f(b[i])).cwiseAbs().maxCoeff());
    m = std::max(m, f(a[i]).cwiseAbs().maxCoeff());
  }
  REQUIRE(m > 0.);
  return d/m;
}

const auto U = [](const PointValues& v){ return v.U; };
const auto A = [](const PointValues& v){ return v.A; };
const auto gradU = [](const PointValues& v){ return v.gradU; };
const auto gradA = [](const PointValues& v){ return v.gradA; };

// Components at rest on every wall: the frequency loader is the stamped one
template<typename Cell>
void same_as_stamps(const std::size_t n, const std::string& u_el){
  constexpr int D = Cell::n_verts - 1;
  CaseDir c("fwall_same");
  write_components<Cell>(c, n, u_el, rest_a<D>, rest_b<D>, blend_lines());
  SimplexInterpol<Cell> stamped(c.file("h5_params.dat"));
  SimplexFreqInterpol<Cell> freq(c.file("freq_params.dat"));
  const auto a = sample(freq, D, T0), b = sample(stamped, D, T0);
  CHECK(gap(a, b, U) < 1e-12);
  CHECK(gap(a, b, A) < 1e-12);
  CHECK(gap(a, b, gradU) < 1e-12);
  CHECK(gap(a, b, gradA) < 1e-12);
  // The rule is in play: without it the wall cells read otherwise
  write_freq_of_stamps(c, u_el, walls, blend_lines(), "wall_upgrade=false\n");
  SimplexFreqInterpol<Cell> plain(c.file("freq_params.dat"));
  const auto p = sample(plain, D, T0);
  CHECK(gap(a, p, gradU) > 1e-3);
  CHECK(gap(a, p, gradA) > 1e-3);
}

// Components that differ on a wall: the frequency loader is the stamped one
// on their sum, with the components in either order, and not the stamped one
// on the components
template<typename Cell>
void same_as_sum(const std::size_t n, const std::string& u_el, const Comp f0, const Comp f1){
  constexpr int D = Cell::n_verts - 1;
  CaseDir c("fwall_comp"), s("fwall_sum");
  write_components<Cell>(c, n, u_el, f0, f1, blend_lines());
  write_sum<Cell>(s, n, u_el, f0, f1, 1. - T0, T0);
  SimplexInterpol<Cell> sum(s.file("h5_params.dat")), stamps(c.file("h5_params.dat"));
  const auto ref = sample(sum, D, T0);
  for (const bool swapped : {false, true}){
    write_freq_of_stamps(c, u_el, walls, blend_lines(swapped), "wall_upgrade=true\n");
    SimplexFreqInterpol<Cell> freq(c.file("freq_params.dat"));
    const auto a = sample(freq, D, T0);
    CHECK(gap(a, ref, U) < 1e-12);
    CHECK(gap(a, ref, gradU) < 1e-12);
    if (!swapped) CHECK(gap(a, sample(stamps, D, T0), gradU) > 1e-9);
  }
}

// Frozen at a time where the moving component's weight vanishes
template<typename Cell>
void frozen_is_live(const std::size_t n, const std::string& u_el){
  constexpr int D = Cell::n_verts - 1;
  const double t_f = 0.3;
  CaseDir c("fwall_frozen"), s("fwall_frozen_sum");
  std::ostringstream lines;
  lines.precision(17);
  lines << "1 0.1 0.7 up_0.h5\n1 " << M_PI/2 - t_f << " 1.3 up_1.h5\n";
  write_components<Cell>(c, n, u_el, rest_a<D>, moving<D>, lines.str());
  SimplexFreqInterpol<Cell> freq(c.file("freq_params.dat"));
  const auto live = sample(freq, D, t_f);
  freq.freeze(t_f);
  for (const double t : {t_f, t_f + 0.3}){
    const auto frozen = sample(freq, D, t);
    CHECK(gap(live, frozen, U) < 1e-13);
    CHECK(gap(live, frozen, gradU) < 1e-13);
  }
  // The sum's own values leave that wall at rest: a stamp of it raises it
  write_sum<Cell>(s, n, u_el, rest_a<D>, moving<D>, 0.7*std::cos(t_f + 0.1), 1.3*std::cos(M_PI/2));
  SimplexInterpol<Cell> sum(s.file("h5_params.dat"));
  CHECK(gap(live, sample(sum, D, t_f), gradU) > 1e-6);
}

// ignore_pressure: the velocity and its derivatives as with the pressure, in
// the cells the rule raises and the others, and no pressure
template<typename Cell>
void without_pressure(const std::size_t n, const std::string& u_el){
  constexpr int D = Cell::n_verts - 1;
  CaseDir c("fwall_nop");
  write_components<Cell>(c, n, u_el, rest_a<D>, moving<D>, blend_lines());
  SimplexFreqInterpol<Cell> with(c.file("freq_params.dat"));
  write_freq_of_stamps(c, u_el, walls, blend_lines(), "wall_upgrade=true\nignore_pressure=true\n");
  SimplexFreqInterpol<Cell> without(c.file("freq_params.dat"));
  // evaluate, which reads the pressure
  const auto a = sample(with, D, T0, false), b = sample(without, D, T0, false);
  CHECK(gap(a, b, U) == 0.);
  CHECK(gap(a, b, A) == 0.);
  CHECK(gap(a, b, gradU) == 0.);
  CHECK(gap(a, b, gradA) == 0.);
  double p = 0.;
  for (std::size_t i = 0; i < a.size(); ++i){
    p = std::max(p, std::abs(a[i].get_p()));
    CHECK(b[i].get_p() == 0.);
  }
  CHECK(p > 0.);
}

}  // namespace

TEST_CASE("Components at rest on the same walls are raised as the stamps they are", "[freq_wall]"){
  SECTION("P1 triangles"){ same_as_stamps<Triangle>(6, "P1"); }
  SECTION("P2 triangles"){ same_as_stamps<Triangle>(6, "P2"); }
  SECTION("P1 tets"){ same_as_stamps<Tet>(3, "P1"); }
  SECTION("P2 tets"){ same_as_stamps<Tet>(3, "P2"); }
}

TEST_CASE("A wall is at rest only where every component is", "[freq_wall]"){
  SECTION("P1 triangles"){ same_as_sum<Triangle>(6, "P1", rest_a<2>, moving<2>); }
  SECTION("P2 triangles"){ same_as_sum<Triangle>(6, "P2", rest_a<2>, moving<2>); }
  SECTION("P1 tets"){ same_as_sum<Tet>(3, "P1", rest_a<3>, moving<3>); }
  SECTION("P2 tets"){ same_as_sum<Tet>(3, "P2", rest_a<3>, moving<3>); }
}

namespace {
// A weak component whose wall values are a strong one's round-off
template<int D>
Vector3d weak(const Vector3d& x){ return 1e-6*rest_b<D>(x) + Vector3d::Constant(1e-17); }
}

TEST_CASE("Rest is tested with one tolerance for the sum", "[freq_wall]"){
  SECTION("P1 triangles"){ same_as_sum<Triangle>(6, "P1", rest_a<2>, weak<2>); }
  SECTION("P2 triangles"){ same_as_sum<Triangle>(6, "P2", rest_a<2>, weak<2>); }
  SECTION("P1 tets"){ same_as_sum<Tet>(3, "P1", rest_a<3>, weak<3>); }
  SECTION("P2 tets"){ same_as_sum<Tet>(3, "P2", rest_a<3>, weak<3>); }
}

TEST_CASE("Frozen fields keep the rest of the components and read as at the freezing time", "[freq_wall]"){
  SECTION("P1 triangles"){ frozen_is_live<Triangle>(6, "P1"); }
  SECTION("P2 triangles"){ frozen_is_live<Triangle>(6, "P2"); }
  SECTION("P1 tets"){ frozen_is_live<Tet>(3, "P1"); }
  SECTION("P2 tets"){ frozen_is_live<Tet>(3, "P2"); }
}

TEST_CASE("ignore_pressure leaves the velocity as it is and gives no pressure", "[freq_wall]"){
  SECTION("P1 triangles"){ without_pressure<Triangle>(6, "P1"); }
  SECTION("P2 triangles"){ without_pressure<Triangle>(6, "P2"); }
  SECTION("P1 tets"){ without_pressure<Tet>(3, "P1"); }
  SECTION("P2 tets"){ without_pressure<Tet>(3, "P2"); }
}

#endif
