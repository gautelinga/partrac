#ifndef __STAMPEDINTERPOL_EVAL_HPP
#define __STAMPEDINTERPOL_EVAL_HPP

// StampedInterpol's evaluation, in a unit of its own per format

#include <array>
#include <cassert>
#include "StampedInterpol.hpp"
#include "p12_eval.hpp"

template<typename Cell, typename Format>
void StampedInterpol<Cell, Format>::wall_motion(const int id, const CellPos& pos, const double alpha_t,
                                                PointValues& fields) const
{
  // Compute P1 basis at x
  std::array<double, Cell::n_dofs_max> _Nu_, _Nux_, _Nuy_, _Nuz_;   // _Nuz_ unused in 2D
  cell_basis(cells_[id], pos.bary, ncoeffs_u, _Nu_.data(), "u");

  std::array<double, Cell::n_dofs_max*3> u_prev_block, u_next_block;
  gather_stamps_nodes<Cell::n_verts, Cell::n_dofs_max, D>(u_dofs_[id], u_dofs_.stride(),
                      u_prev_, u_next_, u_prev_block.data(), u_next_block.data());

  Vector3d U_prev = block_value<D>(_Nu_.data(), u_prev_block.data(), ncoeffs_u);
  Vector3d U_next = block_value<D>(_Nu_.data(), u_next_block.data(), ncoeffs_u);

  Matrix3d gradU_prev = Matrix3d::Zero(), gradU_next = Matrix3d::Zero();
  if (wants_gradient()){
    cell_deriv(cells_[id], pos.bary, ncoeffs_u, _Nux_.data(), _Nuy_.data(), _Nuz_.data(), "u");
    gradU_prev = block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(), u_prev_block.data(), ncoeffs_u);
    gradU_next = block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(), u_next_block.data(), ncoeffs_u);
  }

  // Walls: P2 blocks of the stamps whose wall vertices are at rest
  std::array<double, D*Cell::n_dofs_max> u_prev_block_2, u_next_block_2;
  const WallEdges& w = wall_cells_[wall_index_[id]];
  const bool quad_prev = near_wall::wall_block<Cell>(u_prev_block.data(), u_prev_block_2.data(), w, rest_tol_prev_);
  const bool quad_next = near_wall::wall_block<Cell>(u_next_block.data(), u_next_block_2.data(), w, rest_tol_next_);

  if (quad_prev || quad_next){
    const Uint n2 = Cell::n_dofs_max;
    std::array<double, Cell::n_dofs_max> _Nu2_, _Nu2x_, _Nu2y_, _Nu2z_;   // _Nu2z_ unused in 2D
    cell_basis(cells_[id], pos.bary, n2, _Nu2_.data(), "u");
    if (wants_gradient())
      cell_deriv(cells_[id], pos.bary, n2, _Nu2x_.data(), _Nu2y_.data(), _Nu2z_.data(), "u");

    if (quad_prev){
      U_prev = block_value<D>(_Nu2_.data(), u_prev_block_2.data(), n2);
      if (wants_gradient())
        gradU_prev = block_gradient<D>(_Nu2x_.data(), _Nu2y_.data(), _Nu2z_.data(), u_prev_block_2.data(), n2);
    }
    if (quad_next){
      U_next = block_value<D>(_Nu2_.data(), u_next_block_2.data(), n2);
      if (wants_gradient())
        gradU_next = block_gradient<D>(_Nu2x_.data(), _Nu2y_.data(), _Nu2z_.data(), u_next_block_2.data(), n2);
    }
  }

  fields.U = alpha_t * U_next + (1-alpha_t) * U_prev;
  fields.A = stamp_rate(U_next, U_prev, t_prev, t_next);
  if (wants_gradient()){
    fields.gradU = alpha_t * gradU_next + (1-alpha_t) * gradU_prev;
    fields.gradA = stamp_rate(gradU_next, gradU_prev, t_prev, t_next);
  }
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
  assert(t <= t_next && t >= t_prev);
  const double alpha_t = stamp_weight(t, t_prev, t_next);

  // Quadratic near walls, else P1
  const int id = pos.id;
  if (wall_p2_ == WallP2::Edge && wall_index_[id] >= 0){
    wall_motion(id, pos, alpha_t, fields);
  }
  else {
    // Compute Pk basis at x
    std::array<double, Cell::n_dofs_max> _Nu_, _Nux_, _Nuy_, _Nuz_;   // _Nuz_ unused in 2D
    cell_basis(cells_[id], pos.bary, ncoeffs_u, _Nu_.data(), "u");

    // Gathered by node: the D components of a node are consecutive
    std::array<double, Cell::n_dofs_max*3> u_prev_block, u_next_block;
    gather_stamps_nodes<Cell::n_verts, Cell::n_dofs_max, D>(u_dofs_[id], u_dofs_.stride(),
                        u_prev_, u_next_, u_prev_block.data(), u_next_block.data());

    // Evaluate
    const Vector3d U_prev = block_value<D>(_Nu_.data(), u_prev_block.data(), ncoeffs_u);
    const Vector3d U_next = block_value<D>(_Nu_.data(), u_next_block.data(), ncoeffs_u);
    fields.U = alpha_t * U_next + (1-alpha_t) * U_prev;
    fields.A = stamp_rate(U_next, U_prev, t_prev, t_next);

    if (wants_gradient()){
      cell_deriv(cells_[id], pos.bary, ncoeffs_u, _Nux_.data(), _Nuy_.data(), _Nuz_.data(), "u");
      const Matrix3d gradU_prev = block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(), u_prev_block.data(), ncoeffs_u);
      const Matrix3d gradU_next = block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(), u_next_block.data(), ncoeffs_u);
      fields.gradU = alpha_t * gradU_next + (1-alpha_t) * gradU_prev;
      fields.gradA = stamp_rate(gradU_next, gradU_prev, t_prev, t_next);
    }
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


// The evaluation of one format's two cells
#define STAMPED_EVAL_INSTANCES(Format) \
  template void StampedInterpol<Triangle, Format>::evaluate(const Vector3d&, double, const CellPos&, PointValues&); \
  template void StampedInterpol<Triangle, Format>::evaluate_motion(const Vector3d&, double, const CellPos&, PointValues&); \
  template void StampedInterpol<Tet, Format>::evaluate(const Vector3d&, double, const CellPos&, PointValues&); \
  template void StampedInterpol<Tet, Format>::evaluate_motion(const Vector3d&, double, const CellPos&, PointValues&);

#endif
