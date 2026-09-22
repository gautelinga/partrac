// SplitInterpol's evaluation, in a unit of its own

#include "SplitInterpol.hpp"
#include "p12_eval.hpp"
#include "split_eval.hpp"
#include <array>
#include <cassert>

template<typename Cell>
void SplitInterpol<Cell>::evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  evaluate_impl<true>(x, t, pos, fields);
}

template<typename Cell>
void SplitInterpol<Cell>::evaluate_motion(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  evaluate_impl<false>(x, t, pos, fields);
}

template<typename Cell>
template<bool Scalars>
void SplitInterpol<Cell>::evaluate_impl(const Vector3d &, const double t, const CellPos& pos, PointValues& fields)
{
  // Assuming inside fluid
  assert(t <= t_next && t >= t_prev);
  const double alpha_t = stamp_weight(t, t_prev, t_next);
  const int id = pos.id;
  constexpr std::size_t N = Split::n_sub;
  constexpr std::size_t ni = Split::n_int*std::size_t(D);

  // The sub-cell of the smallest barycentric, and the P2 basis in its own
  // barycentrics
  double mu[Cell::n_verts];
  const int i = split_eval::sub_cell<Cell::n_verts>(pos.bary, mu);
  std::array<double, N> _Nu_;
  std::array<double, N*3> u_prev_block, u_next_block;
  split_eval::sub_basis<Cell>(mu, _Nu_.data());
  // The macro P2 nodes of the facet opposite i, and this cell's interior block
  split_eval::gather_sub_stamps<Cell, D>(i, u_dofs_[std::size_t(id)], u_prev_, u_next_,
                                         int_prev_ + std::size_t(id)*ni,
                                         int_next_ + std::size_t(id)*ni,
                                         u_prev_block.data(), u_next_block.data());

  const Vector3d U_prev = block_value<D>(_Nu_.data(), u_prev_block.data(), N);
  const Vector3d U_next = block_value<D>(_Nu_.data(), u_next_block.data(), N);
  fields.U = alpha_t * U_next + (1-alpha_t) * U_prev;
  fields.A = stamp_rate(U_next, U_prev, t_prev, t_next);

  if (wants_gradient()){
    std::array<double, N> _Nux_, _Nuy_, _Nuz_;   // _Nuz_ unused in 2D
    split_eval::sub_deriv<Cell>(cells_[std::size_t(id)], i, mu,
                                _Nux_.data(), _Nuy_.data(), _Nuz_.data());
    const Matrix3d gradU_prev = block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(),
                                                 u_prev_block.data(), N);
    const Matrix3d gradU_next = block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(),
                                                 u_next_block.data(), N);
    fields.gradU = alpha_t * gradU_next + (1-alpha_t) * gradU_prev;
    fields.gradA = stamp_rate(gradU_next, gradU_prev, t_prev, t_next);
  }

  // The scalars are the macro cell's own, as every other checkpoint loader reads them
  if constexpr (Scalars){
    std::array<double, Cell::n_dofs_max> _Np_;
    if (include_pressure){
      cell_basis(cells_[std::size_t(id)], pos.bary, ncoeffs_p, _Np_.data(), "p");
      std::array<double, Cell::n_dofs_max> p_prev_block, p_next_block;
      gather_stamps_nodes<Cell::n_verts, Cell::n_dofs_max, 1>(p_dofs_[std::size_t(id)],
                          p_dofs_.stride(), p_prev_, p_next_,
                          p_prev_block.data(), p_next_block.data());
      const double P_prev = block_scalar(_Np_.data(), p_prev_block.data(), ncoeffs_p);
      const double P_next = block_scalar(_Np_.data(), p_next_block.data(), ncoeffs_p);
      fields.P = alpha_t * P_next + (1-alpha_t) * P_prev;
    }

    if (include_phi){
      // The pressure's basis where phi is in its space
      std::array<double, Cell::n_dofs_max> _Nphi_;
      const double* Nphi = _Np_.data();
      if (!include_pressure || ncoeffs_phi != ncoeffs_p){
        cell_basis(cells_[std::size_t(id)], pos.bary, ncoeffs_phi, _Nphi_.data(), "phi");
        Nphi = _Nphi_.data();
      }
      std::array<double, Cell::n_dofs_max> phi_prev_block, phi_next_block;
      gather_stamps_nodes<Cell::n_verts, Cell::n_dofs_max, 1>(phi_dofs_[std::size_t(id)],
                          phi_dofs_.stride(), phi_prev_, phi_next_,
                          phi_prev_block.data(), phi_next_block.data());
      const double Phi_prev = block_scalar(Nphi, phi_prev_block.data(), ncoeffs_phi);
      const double Phi_next = block_scalar(Nphi, phi_next_block.data(), ncoeffs_phi);
      fields.Phi = alpha_t * Phi_next + (1-alpha_t) * Phi_prev;
    }

    fields.cell_type = cell_type_[std::size_t(id)];
  }
}

template void SplitInterpol<Triangle>::evaluate(const Vector3d&, double, const CellPos&, PointValues&);
template void SplitInterpol<Triangle>::evaluate_motion(const Vector3d&, double, const CellPos&, PointValues&);
template void SplitInterpol<Tet>::evaluate(const Vector3d&, double, const CellPos&, PointValues&);
template void SplitInterpol<Tet>::evaluate_motion(const Vector3d&, double, const CellPos&, PointValues&);
