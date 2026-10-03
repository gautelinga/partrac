#ifndef __TABLEAU_HPP
#define __TABLEAU_HPP

// Explicit Runge-Kutta steps from a compile-time Butcher tableau

#include <array>
#include <type_traits>
#include <utility>
#include "typedefs.hpp"
#include "TransportElement.hpp"

// A tableau is a type: S stages; row i of a and c over den[i]; b over b_den;
// embedded weights bh over bh_den (0: none); fsal: the last stage is the end
struct RK4Tableau {
  static constexpr int S = 4;
  static constexpr int a[S][S] = {{0, 0, 0, 0}, {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}};
  static constexpr int c[S] = {0, 1, 1, 1};
  static constexpr int den[S] = {1, 2, 2, 1};
  static constexpr int b[S] = {1, 2, 2, 1};
  static constexpr int b_den = 6;
  static constexpr int bh[S] = {0, 0, 0, 0};
  static constexpr int bh_den = 0;
  static constexpr bool fsal = false;
};

namespace rk_detail {

// Explicit, rows summing to c, an FSAL last row equal to b at c = 1
template<class Tab>
constexpr bool valid(){
  for (int i = 0; i < Tab::S; ++i){
    if (Tab::den[i] <= 0) return false;
    int sum = 0;
    for (int j = 0; j < Tab::S; ++j){
      if (j >= i && Tab::a[i][j] != 0) return false;
      sum += Tab::a[i][j];
    }
    if (sum != Tab::c[i]) return false;
  }
  if (Tab::c[0] != 0 || Tab::b_den <= 0 || Tab::bh_den < 0) return false;
  if (Tab::fsal){
    const int l = Tab::S - 1;
    if (Tab::c[l] != Tab::den[l] || Tab::den[l] != Tab::b_den) return false;
    for (int j = 0; j < Tab::S; ++j)
      if (Tab::a[l][j] != Tab::b[j]) return false;
  }
  return true;
}

// Indices of a row's nonzero numerators
template<int S>
struct Nonzero { std::array<int, S> idx{}; int n = 0; };

// Row I of a (I < S), the weights b (I == S) or bh (I == S + 1)
template<class Tab, int I>
constexpr int row_num(const int j){
  if constexpr (I < Tab::S) return Tab::a[I][j];
  else if constexpr (I == Tab::S) return Tab::b[j];
  else return Tab::bh[j];
}

template<class Tab, int I>
constexpr int row_den(){
  if constexpr (I < Tab::S) return Tab::den[I];
  else if constexpr (I == Tab::S) return Tab::b_den;
  else return Tab::bh_den;
}

template<class Tab, int I>
constexpr Nonzero<Tab::S> row_nonzero(){
  Nonzero<Tab::S> r;
  for (int j = 0; j < Tab::S; ++j)
    if (row_num<Tab, I>(j) != 0) r.idx[r.n++] = j;
  return r;
}

template<class Tab, int I, std::size_t... P>
constexpr auto row_seq(std::index_sequence<P...>){
  constexpr Nonzero<Tab::S> nz = row_nonzero<Tab, I>();
  return std::integer_sequence<int, nz.idx[P]...>{};
}

template<class Tab, int I>
struct Row {
  static constexpr int num(const int j){ return row_num<Tab, I>(j); }
  static constexpr int den(){ return row_den<Tab, I>(); }
  static constexpr Nonzero<Tab::S> nz = row_nonzero<Tab, I>();
  using seq = decltype(row_seq<Tab, I>(std::make_index_sequence<row_nonzero<Tab, I>().n>{}));
};

// num k, a unit numerator dropped
template<int Num, class K>
inline decltype(auto) term(const K& k){
  if constexpr (Num == 1) return (k);
  else return double(Num) * k;
}

// Sum of num_j k_j, left to right in ascending j
template<class R, class Arr, int... Js>
inline decltype(auto) wsum(const Arr& k, std::integer_sequence<int, Js...>){
  return (... + term<R::num(Js)>(k[Js]));
}

// (sum num_j k_j) h / den, a unit denominator dropped
template<class R, class Arr>
inline auto increment(const Arr& k, const double h){
  if constexpr (R::den() == 1) return wsum<R>(k, typename R::seq{}) * h;
  else return wsum<R>(k, typename R::seq{}) * h / double(R::den());
}

// y + the row's increment
template<class R, class Y, class Arr>
inline decltype(auto) point(const Y& y, const Arr& k, const double h){
  if constexpr (R::nz.n == 0) return (y);
  else return y + increment<R>(k, h);
}

// t + c h / den
template<class Tab, int I>
inline double stage_time(const double t, const double h){
  constexpr int c = Tab::c[I], d = Tab::den[I];
  if constexpr (c == 0) return t;
  else if constexpr (c == d) return t + h;
  else if constexpr (c == 1) return t + h/double(d);
  else return t + double(c)*h/double(d);
}

template<class Fn, int... Is>
inline void each_stage(Fn&& f, std::integer_sequence<int, Is...>){
  (f(std::integral_constant<int, Is + 1>{}), ...);
}

// The step; Given: the first stage's u1 and J1 at (x, t) instead of an evaluation
template<class Tab, TransportElement E, bool Stop, bool Given, class Eval>
inline bool stages(Eval& ev, [[maybe_unused]] const Vector3d& u1, [[maybe_unused]] const Matrix3d& J1,
                   const Vector3d& x, [[maybe_unused]] const Vector3d& n, [[maybe_unused]] const Matrix3d& F,
                   const double t, const double h,
                   Vector3d& dx, [[maybe_unused]] Vector3d& el, [[maybe_unused]] Matrix3d& dF,
                   [[maybe_unused]] Vector3d& dxh){
  static_assert(valid<Tab>(), "not an explicit consistent tableau");
  constexpr int S = Tab::S;
  std::array<Vector3d, S> k;
  [[maybe_unused]] std::array<Vector3d, S> Fk;
  [[maybe_unused]] std::array<Matrix3d, S> dFk;
  bool inside = true;
  // Outside: zero rates
  auto zero = [&](const int i){
    k[i] = Vector3d::Zero();
    if constexpr (E == TransportElement::Vector) Fk[i] = Vector3d::Zero();
    if constexpr (E == TransportElement::Tensor) dFk[i] = Matrix3d::Zero();
  };
  if constexpr (Given){
    k[0] = u1;
    if constexpr (E == TransportElement::Vector) Fk[0] = J1 * n;
    if constexpr (E == TransportElement::Tensor) dFk[0] = J1 * F;
  }
  else if (ev(x, t)){
    k[0] = ev.u();
    if constexpr (E == TransportElement::Vector){ const Matrix3d J = ev.J(); Fk[0] = J * n; }
    if constexpr (E == TransportElement::Tensor){ const Matrix3d J = ev.J(); dFk[0] = J * F; }
  }
  else { inside = false; zero(0); }
  each_stage([&](auto ic){
    constexpr int I = decltype(ic)::value;
    using R = Row<Tab, I>;
    if ((!Stop || inside) && ev(point<R>(x, k, h), stage_time<Tab, I>(t, h))){
      k[I] = ev.u();
      if constexpr (E == TransportElement::Vector){ const Matrix3d J = ev.J(); const Vector3d ni = point<R>(n, Fk, h); Fk[I] = J * ni; }
      if constexpr (E == TransportElement::Tensor){ const Matrix3d J = ev.J(); const Matrix3d Fi = point<R>(F, dFk, h); dFk[I] = J * Fi; }
    }
    else { inside = false; zero(I); }
  }, std::make_integer_sequence<int, S - 1>{});
  using B = Row<Tab, S>;
  dx = increment<B>(k, h);
  if constexpr (E == TransportElement::Vector) el = n + increment<B>(Fk, h);
  if constexpr (E == TransportElement::Tensor) dF = increment<B>(dFk, h);
  if constexpr (Tab::bh_den != 0) dxh = increment<Row<Tab, S + 1>>(k, h);
  // FSAL: the last stage was the end
  if constexpr (Tab::fsal) return inside;
  else return inside && ev.end(x + dx, t + h);
}

} // namespace rk_detail

// One step of h by Tab from (x, t): dx, the element's el (n plus its
// increment, unnormalised) or its increment dF; false if a stage or the end
// is outside, with zero rates for an outside stage. Stop: no stage after an
// outside one. Eval: ev(x, t) evaluates a stage and says it is inside, ev.u()
// and ev.J() its values, ev.end(x, t) tests the end
template<class Tab, TransportElement E, bool Stop, class Eval>
inline bool rk_stages(Eval& ev, const Vector3d& x, const Vector3d& n, const Matrix3d& F,
                      const double t, const double h, Vector3d& dx, Vector3d& el, Matrix3d& dF){
  Vector3d dxh;
  return rk_detail::stages<Tab, E, Stop, false>(ev, x, Matrix3d(), x, n, F, t, h, dx, el, dF, dxh);
}

// The same with the embedded increment dxh
template<class Tab, TransportElement E, bool Stop, class Eval>
inline bool rk_stages(Eval& ev, const Vector3d& x, const Vector3d& n, const Matrix3d& F,
                      const double t, const double h, Vector3d& dx, Vector3d& el, Matrix3d& dF, Vector3d& dxh){
  static_assert(Tab::bh_den != 0, "no embedded weights");
  return rk_detail::stages<Tab, E, Stop, false>(ev, x, Matrix3d(), x, n, F, t, h, dx, el, dF, dxh);
}

// From a first stage evaluated at (x, t), inside: its u1 and J1
template<class Tab, TransportElement E, bool Stop, class Eval>
inline bool rk_stages_from(Eval& ev, const Vector3d& u1, const Matrix3d& J1, const Vector3d& x,
                           const Vector3d& n, const Matrix3d& F, const double t, const double h,
                           Vector3d& dx, Vector3d& el, Matrix3d& dF){
  Vector3d dxh;
  return rk_detail::stages<Tab, E, Stop, true>(ev, u1, J1, x, n, F, t, h, dx, el, dF, dxh);
}

#endif
