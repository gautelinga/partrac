#ifndef __HELD_EVAL_HPP
#define __HELD_EVAL_HPP

// Evaluations held in one region of a mesh loader: bound to the region at the
// stamps in play, it gathers the cell's nodes once; each point then costs its
// levels, its basis and their products

#include <array>
#include <type_traits>
#include <utility>
#include "PointValues.hpp"
#include "typedefs.hpp"

// The interpolators with regions: the mesh loaders with a held evaluation
template<typename Interp, typename = void>
struct has_regions : std::false_type {};
template<typename Interp>
struct has_regions<Interp, std::void_t<typename Interp::Held>> : std::true_type {};

// The loaders with a held evaluation without the gradients
template<typename Interp, typename = void>
struct has_held_velocity : std::false_type {};
template<typename Interp>
struct has_held_velocity<Interp, std::void_t<decltype(std::declval<Interp&>().held_velocity(
  0, std::declval<const typename Interp::levels_type&>(), 0., std::declval<const typename Interp::Held&>(),
  std::declval<PointValues&>()))>> : std::true_type {};

// The gradients as the run asks for them (wants_gradient)
template<typename Interp>
class HeldEval {
public:
  static constexpr int n_levels = Interp::n_levels;
  using region_type = typename Interp::region_type;
  using levels_type = typename Interp::levels_type;
  HeldEval(Interp& intp, PointValues& pv) : intp_(intp), pv_(pv) {}
  void bind(const region_type& R){
    R_ = R;
    intp_.hold(R, held_);
  }
  // The values at (x, t), inside or not
  bool at(const Vector3d& x, const double t){
    intp_.region_point(R_, x, lev_);
    intp_.held_motion(R_.id, lev_, t, held_, pv_);
    return inside();
  }
  // A stage: the values if inside
  bool operator()(const Vector3d& x, const double t){
    intp_.region_point(R_, x, lev_);
    if (!inside())
      return false;
    intp_.held_motion(R_.id, lev_, t, held_, pv_);
    return true;
  }
  // A stage's velocity and its rate, if inside; the gradients where the loader can leave them out
  bool velocity(const Vector3d& x, const double t){
    if constexpr (!has_held_velocity<Interp>::value)
      return (*this)(x, t);
    else {
      intp_.region_point(R_, x, lev_);
      if (!inside())
        return false;
      intp_.held_velocity(R_.id, lev_, t, held_, pv_);
      return true;
    }
  }
  Vector3d u(){ return pv_.get_u(); }
  Matrix3d J(){ return pv_.get_J(); }
  bool end(const Vector3d& x, const double){
    intp_.region_point(R_, x, lev_);
    return inside();
  }
  // At the last point
  const levels_type& levels() const { return lev_; }
  double band = 0.;   // levels down to -band are inside
  int skip = -1;      // a facet whose level is not tested
private:
  bool inside() const {
    for (int k = 0; k < n_levels; ++k)
      if (k != skip && lev_[k] < -band) return false;
    return true;
  }
  Interp& intp_;
  PointValues& pv_;
  typename Interp::Held held_;
  region_type R_;
  levels_type lev_;
};

#endif
