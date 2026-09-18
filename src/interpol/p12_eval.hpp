#ifdef USE_DOLFIN
#ifndef __P12_EVAL_HPP
#define __P12_EVAL_HPP

// The pieces of a P1 or P2 evaluation on one cell: the dofs of a field, the
// value of a vector field and its gradient. A vector field's components are
// consecutive blocks of ncoeffs in the gathered dofs.

#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
#include "typedefs.hpp"

template<std::size_t N>
inline void gather_n(const std::uint32_t* dofs, const std::vector<double>& prev,
                     const std::vector<double>& next, double* block_prev, double* block_next){
  for (std::size_t i = 0; i < N; ++i){
    block_prev[i] = prev[dofs[i]];
    block_next[i] = next[dofs[i]];
  }
}

// A field's dofs in this cell, at both stamps; P1 and P2 gather a fixed count
template<std::size_t N1, std::size_t N2>
inline void gather_stamps(const std::uint32_t* dofs, const std::size_t stride,
                          const std::vector<double>& prev, const std::vector<double>& next,
                          double* block_prev, double* block_next){
  if (stride == N1)
    gather_n<N1>(dofs, prev, next, block_prev, block_next);
  else
    gather_n<N2>(dofs, prev, next, block_prev, block_next);
}

inline double block_scalar(const double* N, const double* block, const Uint ncoeffs){
  return std::inner_product(N, N + ncoeffs, block, 0.0);
}

template<int Dim>
inline Vector3d block_value(const double* N, const double* block, const Uint ncoeffs){
  const double Ux = std::inner_product(N, N + ncoeffs, block, 0.0);
  const double Uy = std::inner_product(N, N + ncoeffs, block + ncoeffs, 0.0);
  if constexpr (Dim == 2) return {Ux, Uy, 0.};
  else return {Ux, Uy, std::inner_product(N, N + ncoeffs, block + 2*ncoeffs, 0.0)};
}

// gradU(i, j) is dU_i/dx_j
template<int Dim>
__attribute__((always_inline))
inline Matrix3d block_gradient(const double* dNx, const double* dNy, const double* dNz,
                               const double* block, const Uint ncoeffs){
  const double* dN[3] = {dNx, dNy, dNz};
  const auto d = [&](const int j, const int i){
    return std::inner_product(dN[j], dN[j] + ncoeffs, block + i*ncoeffs, 0.0);
  };
  Matrix3d gradU;
  if constexpr (Dim == 2)
    gradU << d(0, 0), d(1, 0), 0.,
             d(0, 1), d(1, 1), 0.,
             0.,      0.,      0.;
  else
    gradU << d(0, 0), d(1, 0), d(2, 0),
             d(0, 1), d(1, 1), d(2, 1),
             d(0, 2), d(1, 2), d(2, 2);
  return gradU;
}

// P1 or P2 by the number of coefficients; what names it in the error
template<typename Cell>
inline void cell_basis(const Cell& cell, const std::array<double, 4>& bary,
                       const Uint ncoeffs, double* N, const char* what){
  if (ncoeffs == Cell::n_verts){
    if constexpr (Cell::n_verts == 3) cell.linearbasis(bary[0], bary[1], bary[2], N);
    else                              cell.linearbasis(bary[0], bary[1], bary[2], bary[3], N);
  }
  else if (ncoeffs == Cell::n_dofs_max){
    if constexpr (Cell::n_verts == 3) cell.quadbasis(bary[0], bary[1], bary[2], N);
    else                              cell.quadbasis(bary[0], bary[1], bary[2], bary[3], N);
  }
  else {
    std::cout << "Unrecognized ncoeffs_" << what << " = " << ncoeffs << std::endl;
    exit(1);
  }
}

template<typename Cell>
inline void cell_deriv(const Cell& cell, const std::array<double, 4>& bary,
                       const Uint ncoeffs, double* dNx, double* dNy, double* dNz,
                       const char* what){   // dNz unused in 2D
  if (ncoeffs == Cell::n_verts){
    if constexpr (Cell::n_verts == 3) cell.linearderiv(bary[0], bary[1], bary[2], dNx, dNy);
    else                              cell.linearderiv(bary[0], bary[1], bary[2], bary[3], dNx, dNy, dNz);
  }
  else if (ncoeffs == Cell::n_dofs_max){
    if constexpr (Cell::n_verts == 3) cell.quadderiv(bary[0], bary[1], bary[2], dNx, dNy);
    else                              cell.quadderiv(bary[0], bary[1], bary[2], bary[3], dNx, dNy, dNz);
  }
  else {
    std::cout << "Unrecognized ncoeffs_" << what << " = " << ncoeffs << std::endl;
    exit(1);
  }
}

#endif
#endif
