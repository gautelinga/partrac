#ifndef __STAMPEDINTERPOL_EVAL_HPP
#define __STAMPEDINTERPOL_EVAL_HPP

// StampedInterpol's evaluation, located and held, in units of their own per format

#include <array>
#include <cassert>
#include "StampedInterpol.hpp"
#include "p12_eval.hpp"

// A cell without the near-wall rule: velocity and its rate at bary from the gathered nodes
template<typename Cell, typename Format>
template<bool Velocity>
__attribute__((always_inline))
inline void StampedInterpol<Cell, Format>::plain_at(const int id, const std::array<double, 4>& bary,
                                                    const double alpha_t, const double* prev,
                                                    const double* next, PointValues& fields) const
{
  // Pk basis at bary
  std::array<double, Cell::n_dofs_max> _Nu_, _Nux_, _Nuy_, _Nuz_;   // _Nuz_ unused in 2D
  cell_basis(cells_[id], bary, ncoeffs_u, _Nu_.data(), "u");

  // Evaluate
  const Vector3d U_prev = block_value<D>(_Nu_.data(), prev, ncoeffs_u);
  const Vector3d U_next = block_value<D>(_Nu_.data(), next, ncoeffs_u);
  fields.U = alpha_t * U_next + (1-alpha_t) * U_prev;
  fields.A = stamp_rate(U_next, U_prev, t_prev, t_next);

  if (!Velocity && wants_gradient()){
    cell_deriv(cells_[id], bary, ncoeffs_u, _Nux_.data(), _Nuy_.data(), _Nuz_.data(), "u");
    const Matrix3d gradU_prev = block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(), prev, ncoeffs_u);
    const Matrix3d gradU_next = block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(), next, ncoeffs_u);
    fields.gradU = alpha_t * gradU_next + (1-alpha_t) * gradU_prev;
    fields.gradA = stamp_rate(gradU_next, gradU_prev, t_prev, t_next);
  }
}

// The P2 blocks of the stamps whose wall vertices are at rest
template<typename Cell, typename Format>
__attribute__((always_inline))
inline std::array<bool, 2> StampedInterpol<Cell, Format>::wall_blocks(const int id, const double* prev,
                                                                      const double* next, double* prev2,
                                                                      double* next2) const
{
  const WallEdges& w = wall_cells_[wall_index_[id]];
  const bool quad_prev = near_wall::wall_block<Cell>(prev, prev2, w, rest_tol_prev_);
  const bool quad_next = near_wall::wall_block<Cell>(next, next2, w, rest_tol_next_);
  return {quad_prev, quad_next};
}

// A wall cell: P2 for a stamp whose wall vertices are at rest, else P1
template<typename Cell, typename Format>
template<bool Velocity>
__attribute__((always_inline))
inline void StampedInterpol<Cell, Format>::wall_at(const int id, const std::array<double, 4>& bary,
                                                   const double alpha_t, const double* prev,
                                                   const double* next, const double* prev2,
                                                   const double* next2, const bool quad_prev,
                                                   const bool quad_next, PointValues& fields) const
{
  if (!(quad_prev || quad_next)){
    plain_at<Velocity>(id, bary, alpha_t, prev, next, fields);
    return;
  }
  // P1 only for a stamp not at rest
  const Uint n2 = Cell::n_dofs_max;
  const bool p1 = !(quad_prev && quad_next);
  std::array<double, Cell::n_dofs_max> _Nu_, _Nux_, _Nuy_, _Nuz_;       // _Nuz_ unused in 2D
  std::array<double, Cell::n_dofs_max> _Nu2_, _Nu2x_, _Nu2y_, _Nu2z_;   // _Nu2z_ unused in 2D
  if (p1) cell_basis(cells_[id], bary, ncoeffs_u, _Nu_.data(), "u");
  cell_basis(cells_[id], bary, n2, _Nu2_.data(), "u");

  const Vector3d U_prev = quad_prev ? block_value<D>(_Nu2_.data(), prev2, n2)
                                    : block_value<D>(_Nu_.data(), prev, ncoeffs_u);
  const Vector3d U_next = quad_next ? block_value<D>(_Nu2_.data(), next2, n2)
                                    : block_value<D>(_Nu_.data(), next, ncoeffs_u);
  fields.U = alpha_t * U_next + (1-alpha_t) * U_prev;
  fields.A = stamp_rate(U_next, U_prev, t_prev, t_next);

  if (!Velocity && wants_gradient()){
    if (p1) cell_deriv(cells_[id], bary, ncoeffs_u, _Nux_.data(), _Nuy_.data(), _Nuz_.data(), "u");
    cell_deriv(cells_[id], bary, n2, _Nu2x_.data(), _Nu2y_.data(), _Nu2z_.data(), "u");
    const Matrix3d gradU_prev = quad_prev
      ? block_gradient<D>(_Nu2x_.data(), _Nu2y_.data(), _Nu2z_.data(), prev2, n2)
      : block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(), prev, ncoeffs_u);
    const Matrix3d gradU_next = quad_next
      ? block_gradient<D>(_Nu2x_.data(), _Nu2y_.data(), _Nu2z_.data(), next2, n2)
      : block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(), next, ncoeffs_u);
    fields.gradU = alpha_t * gradU_next + (1-alpha_t) * gradU_prev;
    fields.gradA = stamp_rate(gradU_next, gradU_prev, t_prev, t_next);
  }
}

template<typename Cell, typename Format>
void StampedInterpol<Cell, Format>::wall_motion(const int id, const CellPos& pos, const double alpha_t,
                                                PointValues& fields) const
{
  std::array<double, Cell::n_dofs_max*3> u_prev_block, u_next_block;
  gather_stamps_nodes<Cell::n_verts, Cell::n_dofs_max, D>(u_dofs_[id], u_dofs_.stride(),
                      u_prev_, u_next_, u_prev_block.data(), u_next_block.data());
  std::array<double, D*Cell::n_dofs_max> u_prev_block_2, u_next_block_2;
  const std::array<bool, 2> quad = wall_blocks(id, u_prev_block.data(), u_next_block.data(),
                                               u_prev_block_2.data(), u_next_block_2.data());
  wall_at(id, pos.bary, alpha_t, u_prev_block.data(), u_next_block.data(), u_prev_block_2.data(),
          u_next_block_2.data(), quad[0], quad[1], fields);
}

template<typename Cell, typename Format>
void StampedInterpol<Cell, Format>::hold(const Region& R, Held& h) const
{
  const int id = R.id;
  h.wall = wall_p2_ == WallP2::Edge && wall_index_[id] >= 0;
  gather_stamps_nodes<Cell::n_verts, Cell::n_dofs_max, D>(u_dofs_[id], u_dofs_.stride(),
                      u_prev_, u_next_, h.prev.data(), h.next.data());
  if (h.wall){
    const std::array<bool, 2> quad = wall_blocks(id, h.prev.data(), h.next.data(), h.prev2.data(), h.next2.data());
    h.quad_prev = quad[0];
    h.quad_next = quad[1];
  }
}

template<typename Cell, typename Format>
template<bool Velocity>
__attribute__((always_inline))
inline void StampedInterpol<Cell, Format>::held_at(const int id, const std::array<double, 4>& lev,
                                                   const double t, const Held& h, PointValues& fields) const
{
  assert(in_bracket(t, t_prev, t_next, this->stamp_snap));
  const double alpha_t = stamp_weight(t, t_prev, t_next);
  if (h.wall)
    wall_at<Velocity>(id, lev, alpha_t, h.prev.data(), h.next.data(), h.prev2.data(), h.next2.data(),
                      h.quad_prev, h.quad_next, fields);
  else
    plain_at<Velocity>(id, lev, alpha_t, h.prev.data(), h.next.data(), fields);
}

template<typename Cell, typename Format>
void StampedInterpol<Cell, Format>::held_motion(const int id, const std::array<double, 4>& lev,
                                                const double t, const Held& h, PointValues& fields)
{
  held_at<false>(id, lev, t, h, fields);
}

template<typename Cell, typename Format>
void StampedInterpol<Cell, Format>::held_velocity(const int id, const std::array<double, 4>& lev,
                                                  const double t, const Held& h, PointValues& fields)
{
  held_at<true>(id, lev, t, h, fields);
}

template<typename Cell, typename Format>
void StampedInterpol<Cell, Format>::evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  evaluate_impl<true>(x, t, pos, fields);
}

template<typename Cell, typename Format>
void StampedInterpol<Cell, Format>::evaluate_motion(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  evaluate_impl<false>(x, t, pos, fields);
}

template<typename Cell, typename Format>
template<bool Scalars>
void StampedInterpol<Cell, Format>::evaluate_impl(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  // Assuming inside fluid
  assert(in_bracket(t, t_prev, t_next, this->stamp_snap));
  const double alpha_t = stamp_weight(t, t_prev, t_next);

  // Quadratic near walls, else P1
  const int id = pos.id;
  if (wall_p2_ == WallP2::Edge && wall_index_[id] >= 0){
    wall_motion(id, pos, alpha_t, fields);
  }
  else {
    // Gathered by node: the D components of a node are consecutive
    std::array<double, Cell::n_dofs_max*3> u_prev_block, u_next_block;
    gather_stamps_nodes<Cell::n_verts, Cell::n_dofs_max, D>(u_dofs_[id], u_dofs_.stride(),
                        u_prev_, u_next_, u_prev_block.data(), u_next_block.data());
    plain_at(id, pos.bary, alpha_t, u_prev_block.data(), u_next_block.data(), fields);
  }

  if constexpr (Scalars){
    std::array<double, Cell::n_dofs_max> _Np_;
    if (include_pressure){
      cell_basis(cells_[id], pos.bary, ncoeffs_p, _Np_.data(), "p");
      const CellDofs& pt = p_table();
      std::array<double, Cell::n_dofs_max> p_prev_block, p_next_block;
      gather_stamps_nodes<Cell::n_verts, Cell::n_dofs_max, 1>(pt[id], pt.stride(),
                          p_prev_, p_next_, p_prev_block.data(), p_next_block.data());
      const double P_prev = block_scalar(_Np_.data(), p_prev_block.data(), ncoeffs_p);
      const double P_next = block_scalar(_Np_.data(), p_next_block.data(), ncoeffs_p);
      fields.P = alpha_t * P_next + (1-alpha_t) * P_prev;
    }

    if (include_phi){
      // The pressure's basis where phi is in its space
      std::array<double, Cell::n_dofs_max> _Nphi_;
      const double* Nphi = _Np_.data();
      if (!include_pressure || ncoeffs_phi != ncoeffs_p){
        cell_basis(cells_[id], pos.bary, ncoeffs_phi, _Nphi_.data(), "phi");
        Nphi = _Nphi_.data();
      }
      const CellDofs& ft = phi_table();
      std::array<double, Cell::n_dofs_max> phi_prev_block, phi_next_block;
      gather_stamps_nodes<Cell::n_verts, Cell::n_dofs_max, 1>(ft[id], ft.stride(),
                          phi_prev_, phi_next_, phi_prev_block.data(), phi_next_block.data());
      const double Phi_prev = block_scalar(Nphi, phi_prev_block.data(), ncoeffs_phi);
      const double Phi_next = block_scalar(Nphi, phi_next_block.data(), ncoeffs_phi);
      fields.Phi = alpha_t * Phi_next + (1-alpha_t) * Phi_prev;
    }

    fields.cell_type = cell_type_[id];
  }
}

template<typename Cell, typename Format>
void StampedInterpol<Cell, Format>::evaluate_phase_gradient(const Vector3d&, const double t, const CellPos& pos,
                                                            Vector3d& g)
{
  g.setZero();
  if constexpr (!Format::vertex_fields) return;
  if (!phase_gradient_) return;
  // P1 on the velocity's vertex table, D a vertex
  const int id = pos.id;
  const double alpha_t = stamp_weight(t, t_prev, t_next);
  std::array<double, Cell::n_dofs_max> N;
  cell_basis(cells_[id], pos.bary, Uint(Cell::n_verts), N.data(), "phi");
  std::array<double, Cell::n_dofs_max*3> prev, next;
  gather_stamps_nodes<Cell::n_verts, Cell::n_dofs_max, D>(u_dofs_[id], u_dofs_.stride(), phi_prev_ + nverts_,
                                                          phi_next_ + nverts_, prev.data(), next.data());
  g = alpha_t*block_value<D>(N.data(), next.data(), Uint(Cell::n_verts))
    + (1 - alpha_t)*block_value<D>(N.data(), prev.data(), Uint(Cell::n_verts));
}

// The evaluation of one format's two cells
#define STAMPED_EVAL_INSTANCES(Format) \
  template void StampedInterpol<Triangle, Format>::evaluate(const Vector3d&, double, const CellPos&, PointValues&); \
  template void StampedInterpol<Triangle, Format>::evaluate_motion(const Vector3d&, double, const CellPos&, PointValues&); \
  template void StampedInterpol<Tet, Format>::evaluate(const Vector3d&, double, const CellPos&, PointValues&); \
  template void StampedInterpol<Tet, Format>::evaluate_motion(const Vector3d&, double, const CellPos&, PointValues&); \
  template void StampedInterpol<Triangle, Format>::evaluate_phase_gradient(const Vector3d&, double, const CellPos&, Vector3d&); \
  template void StampedInterpol<Tet, Format>::evaluate_phase_gradient(const Vector3d&, double, const CellPos&, Vector3d&);

// The held evaluation of one format's two cells, in a unit of its own: apart
// from the located evaluation, whose inlining it would change
#define STAMPED_HELD_INSTANCES(Format) \
  template void StampedInterpol<Triangle, Format>::hold(const Region&, Held&) const; \
  template void StampedInterpol<Tet, Format>::hold(const Region&, Held&) const; \
  template void StampedInterpol<Triangle, Format>::held_motion(int, const std::array<double, 4>&, double, const Held&, PointValues&); \
  template void StampedInterpol<Tet, Format>::held_motion(int, const std::array<double, 4>&, double, const Held&, PointValues&); \
  template void StampedInterpol<Triangle, Format>::held_velocity(int, const std::array<double, 4>&, double, const Held&, PointValues&); \
  template void StampedInterpol<Tet, Format>::held_velocity(int, const std::array<double, 4>&, double, const Held&, PointValues&);

#endif
