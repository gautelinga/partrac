#include "SimplexFreqInterpol.hpp"
#include "near_wall.hpp"
#include "wall_split.hpp"
#include "p12_eval.hpp"
#include <array>
#include <utility>

// A wall cell's evaluation, in a unit of its own: apart from the evaluations,
// whose inlining it would change

namespace {

// The components' blocks of a cell summed with the weights, N nodes a
// component: the velocity's (w) and its rate's (wt)
template<std::size_t N, int D>
inline void sum_blocks(const std::uint32_t* nodes, const std::vector<std::vector<double>>& u,
                       const FreqWeights& wf, const double* blocks, const std::size_t stride,
                       double* __restrict B, double* __restrict Bt){
  for (std::size_t i = 0; i < D*N; ++i) B[i] = Bt[i] = 0.;
  std::array<double, D*N> g;
  for (std::size_t k = 0; k < u.size(); ++k){
    const double* b = blocks ? blocks + k*stride : g.data();
    if (!blocks) gather_node_values_n<N, D>(nodes, u[k].data(), g.data());
    for (std::size_t i = 0; i < D*N; ++i){
      B[i] += wf.w[k]*b[i];
      Bt[i] += wf.wt[k]*b[i];
    }
  }
}

// A thread's last P2 wall cell, ready at the weights' owner and time last asked
template<typename Cell>
using Memo = wall_split::Memo<Cell, std::pair<std::uint64_t, double>>;
template<typename Cell>
Memo<Cell>& memo(){ static thread_local Memo<Cell> m; return m; }

}  // namespace

template<typename Cell>
void SimplexFreqInterpol<Cell>::wall_at(const int id, const std::array<double, 4>& bary, const FreqWeights& wf,
                                        const double* blocks, const std::size_t stride,
                                        PointValues& fields, const int sub) const
{
  constexpr int nv = Cell::n_verts;
  constexpr std::size_t n2 = Cell::n_dofs_max;
  if (!upgrade_p2_){
    const std::size_t w = std::size_t(wall_index_[id]);
    const unsigned rest = rest_[w];
    const bool gradient = wants_gradient();
    // P1: the summed blocks' P2 ones (the edge rule)
    std::array<double, D*n2> B, Bt;
    sum_blocks<nv, D>(u_dofs_[id], u_nodes_, wf, blocks, stride, B.data(), Bt.data());
    std::array<double, D*n2> B2, Bt2;
    near_wall::raise_block<Cell>(B.data(), B2.data(), wall_cells_[w], rest);
    near_wall::raise_block<Cell>(Bt.data(), Bt2.data(), wall_cells_[w], rest);
    std::array<double, n2> N2, N2x, N2y, N2z;   // N2z unused in 2D
    cell_basis(cells_[id], bary, n2, N2.data(), "u");
    fields.U = block_value<D>(N2.data(), B2.data(), n2);
    fields.A = block_value<D>(N2.data(), Bt2.data(), n2);
    if (gradient){
      cell_deriv(cells_[id], bary, n2, N2x.data(), N2y.data(), N2z.data(), "u");
      fields.gradU = block_gradient<D>(N2x.data(), N2y.data(), N2z.data(), B2.data(), n2);
      fields.gradA = block_gradient<D>(N2x.data(), N2y.data(), N2z.data(), Bt2.data(), n2);
    }
    return;
  }

  upgrade_at(id, bary, wf, blocks, stride, fields, sub);
}

template<typename Cell>
void SimplexFreqInterpol<Cell>::upgrade_at(const int id, const std::array<double, 4>& bary, const FreqWeights& wf,
                                           const double* blocks, const std::size_t stride,
                                           PointValues& fields, const int sub) const
{
  constexpr std::size_t n2 = Cell::n_dofs_max;
  const wall_split::WallCell& c = split_.cell[std::size_t(wall_index_[id])];
  Memo<Cell>& m = memo<Cell>();
  const std::pair<std::uint64_t, double> key{wf.owner, wf.t};
  if (m.key != key || m.id != id){
    std::array<double, D*n2> B, Bt;
    sum_blocks<n2, D>(u_dofs_[id], u_nodes_, wf, blocks, stride, B.data(), Bt.data());
    wall_split::ready<Cell>(cells_[id], split_, c, B.data(), Bt.data(), m.ready);
    m.key = key;
    m.id = id;
  }
  Vector3d U = Vector3d::Zero(), A = Vector3d::Zero();
  if (wants_gradient()){
    Matrix3d gradU = Matrix3d::Zero(), gradA = Matrix3d::Zero();
    wall_split::add_ready<Cell, true>(split_, c, m.ready, bary, U, A, gradU, gradA, sub);
    fields.gradU = gradU;
    fields.gradA = gradA;
  }
  else {
    Matrix3d unused;
    wall_split::add_ready<Cell, false>(split_, c, m.ready, bary, U, A, unused, unused, sub);
  }
  fields.U = U;
  fields.A = A;
}

template void SimplexFreqInterpol<Triangle>::wall_at(const int, const std::array<double, 4>&, const FreqWeights&,
                                                     const double*, const std::size_t, PointValues&, const int) const;
template void SimplexFreqInterpol<Tet>::wall_at(const int, const std::array<double, 4>&, const FreqWeights&,
                                                const double*, const std::size_t, PointValues&, const int) const;
template void SimplexFreqInterpol<Triangle>::upgrade_at(const int, const std::array<double, 4>&, const FreqWeights&,
                                                        const double*, const std::size_t, PointValues&, const int) const;
template void SimplexFreqInterpol<Tet>::upgrade_at(const int, const std::array<double, 4>&, const FreqWeights&,
                                                   const double*, const std::size_t, PointValues&, const int) const;
