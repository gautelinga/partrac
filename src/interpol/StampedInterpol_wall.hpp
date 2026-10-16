#ifndef __STAMPEDINTERPOL_WALL_HPP
#define __STAMPEDINTERPOL_WALL_HPP

// StampedInterpol's P2 wall cells, in a unit of their own per format: apart
// from the evaluations, whose inlining they would change

#include <array>
#include <cstdint>
#include "StampedInterpol_eval.hpp"
#include "near_wall.hpp"

// A P2 wall cell: P2 plus the terms, zero for a stamp whose wall vertices are not at rest
template<typename Cell, typename Format>
template<bool Velocity>
__attribute__((noinline))
void StampedInterpol<Cell, Format>::upgrade_at(const int id, const std::array<double, 4>& bary,
                                               const double alpha_t, const double* prev,
                                               const double* next, const near_wall::TermView& c,
                                               PointValues& fields) const
{
  constexpr int nv = Cell::n_verts;
  const Uint n2 = Cell::n_dofs_max;
  const double* cprev = c.prev;
  const double* cnext = c.next;
  const std::uint32_t terms = c.mask;
  std::array<double, Cell::n_dofs_max> _Nu_, _Nux_, _Nuy_, _Nuz_;   // _Nuz_ unused in 2D
  cell_basis(cells_[id], bary, n2, _Nu_.data(), "u");
  const near_wall::TermProducts<Cell> P(bary);

  Vector3d U_prev = block_value<D>(_Nu_.data(), prev, n2);
  Vector3d U_next = block_value<D>(_Nu_.data(), next, n2);
  // The terms
  std::array<double, near_wall::n_terms<Cell>> W;
  near_wall::term_weights<Cell>(bary, P, W.data());
  int j = 0;
  for (std::uint32_t m = terms; m; m &= m - 1, j += D){
    const double Wk = W[std::size_t(__builtin_ctz(m))];
    for (int d = 0; d < D; ++d){
      U_prev[d] += Wk*cprev[j + d];
      U_next[d] += Wk*cnext[j + d];
    }
  }
  fields.U = alpha_t * U_next + (1-alpha_t) * U_prev;
  fields.A = stamp_rate(U_next, U_prev, t_prev, t_next);

  if (!Velocity && wants_gradient()){
    cell_deriv(cells_[id], bary, n2, _Nux_.data(), _Nuy_.data(), _Nuz_.data(), "u");
    Matrix3d gradU_prev = block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(), prev, n2);
    Matrix3d gradU_next = block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(), next, n2);
    // The terms' gradients
    std::array<Vector3d, nv> g;
    for (int q = 0; q < nv; ++q) g[std::size_t(q)] = cells_[id].bary_grad(q);
    int i = 0;
    for (std::uint32_t m = terms; m; m &= m - 1, i += D){
      const Vector3d gk = near_wall::term_gradient<Cell>(__builtin_ctz(m), bary, P, g);
      for (int d = 0; d < D; ++d)
        for (int x = 0; x < D; ++x){
          gradU_prev(d, x) += cprev[i + d]*gk[x];
          gradU_next(d, x) += cnext[i + d]*gk[x];
        }
    }
    fields.gradU = alpha_t * gradU_next + (1-alpha_t) * gradU_prev;
    fields.gradA = stamp_rate(gradU_next, gradU_prev, t_prev, t_next);
  }
}

// The P2 wall cells of one format's two cells
#define STAMPED_WALL_INSTANCES(Format) \
  template void StampedInterpol<Triangle, Format>::upgrade_at<false>(int, const std::array<double, 4>&, double, \
    const double*, const double*, const near_wall::TermView&, PointValues&) const; \
  template void StampedInterpol<Triangle, Format>::upgrade_at<true>(int, const std::array<double, 4>&, double, \
    const double*, const double*, const near_wall::TermView&, PointValues&) const; \
  template void StampedInterpol<Tet, Format>::upgrade_at<false>(int, const std::array<double, 4>&, double, \
    const double*, const double*, const near_wall::TermView&, PointValues&) const; \
  template void StampedInterpol<Tet, Format>::upgrade_at<true>(int, const std::array<double, 4>&, double, \
    const double*, const double*, const near_wall::TermView&, PointValues&) const;

#endif
