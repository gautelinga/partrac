// The compile-time Butcher tableau and its stage loop (tableau.hpp).
//
// RK4's table must round exactly as RK4 written out stage by stage, kept here
// as the reference, and both are run on random states through a mock
// interpolator whose fluid is a ball, so that stages and ends fall outside
// too, with and without Stop. This unit is built with the production flags,
// which is where the bits matter.
//
// Then each table's order on linear ODEs with known solutions: x' = M x
// with its matrix exponential, carrying n' = M n and F' = M F as the vector
// and tensor elements do, and a quadrature x' = p(t) that a scheme of order
// q integrates exactly for deg p < q and not at deg p = q, which checks the
// nodes c. Beside RK4: Kutta's 3/8 rule (negative numerators, a nonzero
// a_31) and Bogacki-Shampine 3(2) (embedded weights, FSAL), tables defined
// here only, since the code ships RK4 alone.
#include <catch2/catch.hpp>

#include <cmath>
#include <cstring>
#include <random>

#include "steps_impl.hpp"
#include "tableau.hpp"

namespace {

struct Kutta38 {
  static constexpr int S = 4;
  static constexpr int a[S][S] = {{0, 0, 0, 0}, {1, 0, 0, 0}, {-1, 3, 0, 0}, {1, -1, 1, 0}};
  static constexpr int c[S] = {0, 1, 2, 1};
  static constexpr int den[S] = {1, 3, 3, 1};
  static constexpr int b[S] = {1, 3, 3, 1};
  static constexpr int b_den = 8;
  static constexpr int bh[S] = {0, 0, 0, 0};
  static constexpr int bh_den = 0;
  static constexpr bool fsal = false;
};

struct BS32 {
  static constexpr int S = 4;
  static constexpr int a[S][S] = {{0, 0, 0, 0}, {1, 0, 0, 0}, {0, 3, 0, 0}, {2, 3, 4, 0}};
  static constexpr int c[S] = {0, 1, 3, 9};
  static constexpr int den[S] = {1, 2, 4, 9};
  static constexpr int b[S] = {2, 3, 4, 0};
  static constexpr int b_den = 9;
  static constexpr int bh[S] = {7, 6, 8, 3};
  static constexpr int bh_den = 24;
  static constexpr bool fsal = true;
};

// A smooth unsteady field, fluid inside a ball of radius 1; the cell id a
// hash of the point, so a located end can be compared too
struct BallInterp {
  double get_U0() const { return 1.3; }
  bool locate(const Vector3d& x, const double t, CellPos& pos){
    if (x.squaredNorm() >= 1.) return false;
    pos.id = int(std::floor(8*x[0])) + 16*int(std::floor(8*x[1])) + 256*int(std::floor(8*x[2] + t));
    return true;
  }
  void evaluate(const Vector3d& x, const double t, const CellPos& pos, PointValues& pv){
    const double s = std::sin(x[1] + t), c = std::cos(x[2] - 0.5*t);
    pv.U = {s + 0.3*x[2], c*x[0] - 0.2, x[0]*x[1]*(1 + t)};
    pv.gradU << 0., std::cos(x[1] + t), 0.3,
                c, 0., -std::sin(x[2] - 0.5*t)*x[0],
                x[1]*(1 + t), x[0]*(1 + t), 0.;
  }
};

// RK4 written out stage by stage: the reference
template<TransportElement E, bool Stop, typename Interp>
inline bool rk4_stages_reference(Interp& intp, PointValues& ptvals, CellPos& pos, const Vector3d& x,
                                 [[maybe_unused]] const Vector3d& n, [[maybe_unused]] const Matrix3d& F,
                                 const double t, const double h,
                                 Vector3d& dx, [[maybe_unused]] Vector3d& el, [[maybe_unused]] Matrix3d& dF){
  Vector3d k1, k2, k3, k4;
  [[maybe_unused]] Vector3d F1, F2, F3, F4;
  [[maybe_unused]] Matrix3d dF1, dF2, dF3, dF4;
  bool inside = true;
  k1 = k2 = k3 = k4 = Vector3d::Zero();
  if constexpr (E == TransportElement::Vector) F1 = F2 = F3 = F4 = Vector3d::Zero();
  if constexpr (E == TransportElement::Tensor) dF1 = dF2 = dF3 = dF4 = Matrix3d::Zero();
  if (intp.locate(x, t, pos)){
    evaluate_motion(intp, x, t, pos, ptvals);
    k1 = ptvals.get_u();
    if constexpr (E == TransportElement::Vector){ const Matrix3d J = ptvals.get_J(); F1 = J * n; }
    if constexpr (E == TransportElement::Tensor){ const Matrix3d J = ptvals.get_J(); dF1 = J * F; }
  } else inside = false;
  if ((!Stop || inside) && intp.locate(x + k1 * h/2, t + h/2, pos)){
    evaluate_motion(intp, x + k1 * h/2, t + h/2, pos, ptvals);
    k2 = ptvals.get_u();
    if constexpr (E == TransportElement::Vector){ const Matrix3d J = ptvals.get_J(); const Vector3d n2 = n + F1 * h/2; F2 = J * n2; }
    if constexpr (E == TransportElement::Tensor){ const Matrix3d J = ptvals.get_J(); const Matrix3d F2 = F + dF1 * h/2; dF2 = J * F2; }
  } else inside = false;
  if ((!Stop || inside) && intp.locate(x + k2 * h/2, t + h/2, pos)){
    evaluate_motion(intp, x + k2 * h/2, t + h/2, pos, ptvals);
    k3 = ptvals.get_u();
    if constexpr (E == TransportElement::Vector){ const Matrix3d J = ptvals.get_J(); const Vector3d n3 = n + F2 * h/2; F3 = J * n3; }
    if constexpr (E == TransportElement::Tensor){ const Matrix3d J = ptvals.get_J(); const Matrix3d F3 = F + dF2 * h/2; dF3 = J * F3; }
  } else inside = false;
  if ((!Stop || inside) && intp.locate(x + k3 * h, t + h, pos)){
    evaluate_motion(intp, x + k3 * h, t + h, pos, ptvals);
    k4 = ptvals.get_u();
    if constexpr (E == TransportElement::Vector){ const Matrix3d J = ptvals.get_J(); const Vector3d n4 = n + F3 * h; F4 = J * n4; }
    if constexpr (E == TransportElement::Tensor){ const Matrix3d J = ptvals.get_J(); const Matrix3d F4 = F + dF3 * h; dF4 = J * F4; }
  } else inside = false;
  dx = (k1 + 2*k2 + 2*k3 + k4) * h/6;
  if constexpr (E == TransportElement::Vector) el = n + (F1 + 2*F2 + 2*F3 + F4) * h/6;
  if constexpr (E == TransportElement::Tensor) dF = (dF1 + 2*dF2 + 2*dF3 + dF4) * h/6;
  return inside && intp.locate(x + dx, t + h, pos);
}

template<class M>
bool same_bits(const M& a, const M& b){ return std::memcmp(a.data(), b.data(), sizeof(double)*a.size()) == 0; }

// Random states near the ball's surface, steps long enough to leave it; counts
// the steps that ended outside, to show both paths ran
template<TransportElement E, bool Stop>
int compare_with_reference(const int n_states){
  std::mt19937 rng(12345 + int(E) + 7*Stop);
  std::uniform_real_distribution<double> uni(-1., 1.);
  BallInterp intp;
  int n_out = 0;
  for (int s = 0; s < n_states; ++s){
    const Vector3d x = 0.98*Vector3d(uni(rng), uni(rng), uni(rng)).normalized()*std::cbrt(0.5*(1 + uni(rng)));
    const Vector3d n = Vector3d(uni(rng), uni(rng), uni(rng)).normalized();
    Matrix3d F = Matrix3d::Identity();
    for (int i = 0; i < 9; ++i) F.data()[i] += 0.3*uni(rng);
    const double t = 2*uni(rng), h = 0.1*(1.05 + uni(rng));
    PointValues pv_ref(intp.get_U0()), pv_table(intp.get_U0());
    CellPos pos_ref, pos_table;
    pos_ref.id = pos_table.id = 3;
    Vector3d dx_ref, dx_table, el_ref = Vector3d::Zero(), el_table = Vector3d::Zero();
    Matrix3d dF_ref = Matrix3d::Zero(), dF_table = Matrix3d::Zero();
    const bool in_ref = rk4_stages_reference<E, Stop>(intp, pv_ref, pos_ref, x, n, F, t, h, dx_ref, el_ref, dF_ref);
    LocatedEval<BallInterp> ev{intp, pv_table, pos_table};
    const bool in_table = rk_stages<RK4Tableau, E, Stop>(ev, x, n, F, t, h, dx_table, el_table, dF_table);
    REQUIRE(in_ref == in_table);
    REQUIRE(pos_ref.id == pos_table.id);
    REQUIRE(same_bits(dx_ref, dx_table));
    if constexpr (E == TransportElement::Vector) REQUIRE(same_bits(el_ref, el_table));
    if constexpr (E == TransportElement::Tensor) REQUIRE(same_bits(dF_ref, dF_table));
    n_out += !in_table;
  }
  return n_out;
}

// x' = M x + p(t), always inside
struct LinearEval {
  Matrix3d M;
  int deg = -1;   // p = (deg + 1) t^deg e, none below 0
  Vector3d u_ = Vector3d::Zero();
  bool operator()(const Vector3d& x, const double t){
    u_ = M*x;
    if (deg >= 0) u_ += (deg + 1)*std::pow(t, deg)*Vector3d(1., -2., 0.5);
    return true;
  }
  Vector3d u(){ return u_; }
  Matrix3d J(){ return M; }
  bool end(const Vector3d&, const double){ return true; }
};

// Rotation at w in the xy plane, growth at g along z
Matrix3d rot_growth(const double w, const double g){
  Matrix3d M;
  M << 0., -w, 0., w, 0., 0., 0., 0., g;
  return M;
}
Matrix3d expm_rot_growth(const double w, const double g, const double t){
  Matrix3d E;
  E << std::cos(w*t), -std::sin(w*t), 0., std::sin(w*t), std::cos(w*t), 0., 0., 0., std::exp(g*t);
  return E;
}

// Errors of x, n and F at T = 1 in N steps of x' = M x
template<class Tab>
std::array<double, 3> linear_errors(const int N, const bool embedded = false){
  const double w = 2.1, g = -0.7;
  LinearEval ev{rot_growth(w, g)};
  const Vector3d x0(0.3, -0.5, 0.8);
  const Vector3d n0 = Vector3d(1., 2., -1.).normalized();
  Vector3d x = x0, n = n0;
  Matrix3d F = Matrix3d::Identity();
  const double h = 1./N;
  for (int s = 0; s < N; ++s){
    Vector3d dx, dxh, el, unused;
    Matrix3d dF, unusedF;
    if constexpr (Tab::bh_den != 0){
      rk_stages<Tab, TransportElement::Point, false>(ev, x, n, F, s*h, h, dx, unused, unusedF, dxh);
      if (embedded) dx = dxh;
    } else
      rk_stages<Tab, TransportElement::Point, false>(ev, x, n, F, s*h, h, dx, unused, unusedF);
    rk_stages<Tab, TransportElement::Vector, false>(ev, x, n, F, s*h, h, unused, el, unusedF);
    rk_stages<Tab, TransportElement::Tensor, false>(ev, x, n, F, s*h, h, unused, unused, dF);
    x += dx;
    n = el;
    F += dF;
  }
  const Matrix3d E = expm_rot_growth(w, g, 1.);
  return {(x - E*x0).norm(), (n - E*n0).norm(), (F - E).norm()};
}

// Error at T = 1 of x' = p(t), deg p = deg, x(0) = 0
template<class Tab>
double quadrature_error(const int deg, const int N){
  LinearEval ev{Matrix3d::Zero(), deg};
  Vector3d x = Vector3d::Zero(), el;
  Matrix3d dF;
  const double h = 1./N;
  for (int s = 0; s < N; ++s){
    Vector3d dx;
    rk_stages<Tab, TransportElement::Point, false>(ev, x, x, Matrix3d::Identity(), s*h, h, dx, el, dF);
    x += dx;
  }
  return (x - Vector3d(1., -2., 0.5)).norm();
}

template<class Tab>
void check_order(const int q){
  const auto e1 = linear_errors<Tab>(20), e2 = linear_errors<Tab>(40);
  for (int c = 0; c < 3; ++c){
    INFO("component " << c << ": " << e1[c] << " -> " << e2[c]);
    const double rate = std::log2(e1[c]/e2[c]);
    CHECK(rate > q - 0.15);
    CHECK(rate < q + 0.15);
  }
  for (int deg = 0; deg < q; ++deg)
    CHECK(quadrature_error<Tab>(deg, 7) < 1e-13);
  CHECK(quadrature_error<Tab>(q, 7) > 1e-8);
}

} // namespace

TEST_CASE("RK4's table rounds as the hand-written block, Stop or not", "[tableau]"){
  const int N = 20000;
  int out = 0;
  out += compare_with_reference<TransportElement::Point, false>(N);
  out += compare_with_reference<TransportElement::Point, true>(N);
  out += compare_with_reference<TransportElement::Vector, false>(N);
  out += compare_with_reference<TransportElement::Vector, true>(N);
  out += compare_with_reference<TransportElement::Tensor, false>(N);
  out += compare_with_reference<TransportElement::Tensor, true>(N);
  // both paths ran
  CHECK(out > 6*N/20);
  CHECK(out < 6*N*19/20);
}

TEST_CASE("Each table has its order on linear ODEs", "[tableau]"){
  SECTION("RK4"){ check_order<RK4Tableau>(4); }
  SECTION("Kutta's 3/8"){ check_order<Kutta38>(4); }
  SECTION("Bogacki-Shampine 3(2)"){ check_order<BS32>(3); }
}

TEST_CASE("Embedded weights have their own order", "[tableau]"){
  const auto e1 = linear_errors<BS32>(20, true), e2 = linear_errors<BS32>(40, true);
  const double rate = std::log2(e1[0]/e2[0]);
  CHECK(rate > 1.85);
  CHECK(rate < 2.15);
}

TEST_CASE("A first stage passed in is the stage evaluated", "[tableau]"){
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> uni(-1., 1.);
  BallInterp intp;
  for (int s = 0; s < 200; ++s){
    const Vector3d x = 0.5*Vector3d(uni(rng), uni(rng), uni(rng));
    const Vector3d n = Vector3d(uni(rng), uni(rng), 1.).normalized();
    const Matrix3d F = Matrix3d::Identity() + 0.2*Matrix3d::Random();
    const double t = uni(rng), h = 0.05*(1 + uni(rng));
    PointValues pv(intp.get_U0());
    CellPos pos;
    LocatedEval<BallInterp> ev{intp, pv, pos};
    Vector3d dx1, dx2, el1, el2;
    Matrix3d dF1, dF2;
    const bool in1 = rk_stages<RK4Tableau, TransportElement::Tensor, false>(ev, x, n, F, t, h, dx1, el1, dF1);
    REQUIRE(ev(x, t));
    const Vector3d u1 = ev.u();
    const Matrix3d J1 = ev.J();
    const bool in2 = rk_stages_from<RK4Tableau, TransportElement::Tensor, false>(ev, u1, J1, x, n, F, t, h, dx2, el2, dF2);
    REQUIRE(in1 == in2);
    REQUIRE(same_bits(dx1, dx2));
    REQUIRE(same_bits(dF1, dF2));
    rk_stages<RK4Tableau, TransportElement::Vector, false>(ev, x, n, F, t, h, dx1, el1, dF1);
    rk_stages_from<RK4Tableau, TransportElement::Vector, false>(ev, u1, J1, x, n, F, t, h, dx2, el2, dF2);
    REQUIRE(same_bits(el1, el2));
  }
}

TEST_CASE("An FSAL table's last stage is its end", "[tableau]"){
  // the last stage point is x + dx bit for bit, so its values seed the next step
  LinearEval ev{rot_growth(1.3, 0.4)};
  const Vector3d x(0.2, 0.1, -0.4);
  Vector3d dx, el, dxh;
  Matrix3d dF;
  rk_stages<BS32, TransportElement::Point, false>(ev, x, x, Matrix3d::Identity(), 0., 0.1, dx, el, dF, dxh);
  const Vector3d u_last = ev.u();
  REQUIRE(ev(x + dx, 0.1));
  CHECK(same_bits(u_last, ev.u()));
}
