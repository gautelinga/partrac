#ifndef __STAMPEDINTERPOL_WALL_HPP
#define __STAMPEDINTERPOL_WALL_HPP

// StampedInterpol's P2 wall cells, in a unit of their own per format: apart
// from the evaluations, whose inlining they would change

#include <array>
#include <cstdint>
#include "StampedInterpol_eval.hpp"
#include "wall_split.hpp"

// A P2 wall cell made ready from both stamps
template<typename Cell, typename Format>
void StampedInterpol<Cell, Format>::wall_hold(const int id, Held& h) const
{
  wall_split::ready<Cell>(cells_[id], split_, split_.cell[std::size_t(wall_index_[id])], h.prev.data(),
                          h.next.data(), h.ready);
}

// A P2 wall cell: P2 plus the wall correction at each stamp
template<typename Cell, typename Format>
template<bool Velocity, bool Hold>
__attribute__((noinline))
void StampedInterpol<Cell, Format>::upgrade_at(const int id, const std::array<double, 4>& bary,
                                               const double alpha_t, const Held* h, PointValues& fields) const
{
  const Cell& cell = cells_[id];
  const wall_split::WallCell& c = split_.cell[std::size_t(wall_index_[id])];
  const wall_split::Ready<Cell>* r;
  if constexpr (Hold) r = &h->ready;
  else {
    WallMemo& m = wall_memo();
    if (m.key != stamps_in_play_ || m.id != id){
      std::array<double, Cell::n_dofs_max*3> u_prev, u_next;
      gather_stamps_nodes<Cell::n_verts, Cell::n_dofs_max, D>(u_dofs_[id], u_dofs_.stride(),
                          u_prev_, u_next_, u_prev.data(), u_next.data());
      wall_split::ready<Cell>(cell, split_, c, u_prev.data(), u_next.data(), m.ready);
      m.key = stamps_in_play_;
      m.id = id;
    }
    r = &m.ready;
  }
  const int sub = Hold ? h->sub : -1;
  Vector3d U_prev = Vector3d::Zero(), U_next = Vector3d::Zero();
  if (!Velocity && wants_gradient()){
    Matrix3d gradU_prev = Matrix3d::Zero(), gradU_next = Matrix3d::Zero();
    wall_split::add_ready<Cell, true>(split_, c, *r, bary, U_prev, U_next, gradU_prev, gradU_next, sub);
    fields.gradU = alpha_t * gradU_next + (1-alpha_t) * gradU_prev;
    fields.gradA = stamp_rate(gradU_next, gradU_prev, t_prev, t_next);
  }
  else {
    Matrix3d unused;
    wall_split::add_ready<Cell, false>(split_, c, *r, bary, U_prev, U_next, unused, unused, sub);
  }
  fields.U = alpha_t * U_next + (1-alpha_t) * U_prev;
  fields.A = stamp_rate(U_next, U_prev, t_prev, t_next);
}

// The P2 wall cells of one format's two cells
#define STAMPED_WALL_INSTANCES(Format) \
  template void StampedInterpol<Triangle, Format>::wall_hold(int, Held&) const; \
  template void StampedInterpol<Tet, Format>::wall_hold(int, Held&) const; \
  template void StampedInterpol<Triangle, Format>::upgrade_at<false, false>(int, const std::array<double, 4>&, \
    double, const Held*, PointValues&) const; \
  template void StampedInterpol<Triangle, Format>::upgrade_at<false, true>(int, const std::array<double, 4>&, \
    double, const Held*, PointValues&) const; \
  template void StampedInterpol<Triangle, Format>::upgrade_at<true, true>(int, const std::array<double, 4>&, \
    double, const Held*, PointValues&) const; \
  template void StampedInterpol<Tet, Format>::upgrade_at<false, false>(int, const std::array<double, 4>&, \
    double, const Held*, PointValues&) const; \
  template void StampedInterpol<Tet, Format>::upgrade_at<false, true>(int, const std::array<double, 4>&, \
    double, const Held*, PointValues&) const; \
  template void StampedInterpol<Tet, Format>::upgrade_at<true, true>(int, const std::array<double, 4>&, \
    double, const Held*, PointValues&) const;

#endif
