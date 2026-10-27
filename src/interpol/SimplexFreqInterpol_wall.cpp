#include "SimplexFreqInterpol.hpp"
#include "near_wall.hpp"
#include "p12_eval.hpp"
#include <array>

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

}  // namespace

template<typename Cell>
void SimplexFreqInterpol<Cell>::wall_at(const int id, const std::array<double, 4>& bary, const FreqWeights& wf,
                                        const double* blocks, const std::size_t stride,
                                        PointValues& fields) const
{
  constexpr int nv = Cell::n_verts;
  constexpr std::size_t n2 = Cell::n_dofs_max;
  const std::size_t w = std::size_t(wall_index_[id]);
  const unsigned rest = rest_[w];
  const bool gradient = wants_gradient();
  std::array<double, D*n2> B, Bt;

  if (!upgrade_p2_){
    // P1: the summed blocks' P2 ones (the edge rule)
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

  sum_blocks<n2, D>(u_dofs_[id], u_nodes_, wf, blocks, stride, B.data(), Bt.data());
  upgrade_at(id, bary, B.data(), Bt.data(), fields);
}

template<typename Cell>
void SimplexFreqInterpol<Cell>::upgrade_at(const int id, const std::array<double, 4>& bary, const double* B,
                                           const double* Bt, PointValues& fields) const
{
  constexpr int nv = Cell::n_verts;
  constexpr std::size_t n2 = Cell::n_dofs_max;
  const std::size_t w = std::size_t(wall_index_[id]);
  const unsigned rest = rest_[w];
  const bool gradient = wants_gradient();
  const near_wall::UpgradeCell& wc = upgrade_.cells[w];
  const std::uint32_t terms = upgrade_.spans[w].mask;
  constexpr int nt = near_wall::n_terms<Cell>;
  std::array<double, D*nt> c, ct;
  near_wall::upgrade_coefficients<Cell>(B, c.data(), wc, terms, upgrade_.normals.data() + wc.normals, rest);
  near_wall::upgrade_coefficients<Cell>(Bt, ct.data(), wc, terms, upgrade_.normals.data() + wc.normals, rest);
  std::array<double, n2> N, Nx, Ny, Nz;   // Nz unused in 2D
  cell_basis(cells_[id], bary, n2, N.data(), "u");
  const near_wall::TermProducts<Cell> P(bary);
  std::array<double, nt> W;
  near_wall::term_weights<Cell>(bary, P, W.data());
  Vector3d U = block_value<D>(N.data(), B, n2);
  Vector3d A = block_value<D>(N.data(), Bt, n2);
  near_wall::add_terms<Cell>(W.data(), terms, c.data(), ct.data(), U, A);
  fields.U = U;
  fields.A = A;
  if (gradient){
    cell_deriv(cells_[id], bary, n2, Nx.data(), Ny.data(), Nz.data(), "u");
    Matrix3d gradU = block_gradient<D>(Nx.data(), Ny.data(), Nz.data(), B, n2);
    Matrix3d gradA = block_gradient<D>(Nx.data(), Ny.data(), Nz.data(), Bt, n2);
    std::array<Vector3d, nv> g;
    for (int q = 0; q < nv; ++q) g[std::size_t(q)] = cells_[id].bary_grad(q);
    near_wall::add_term_gradients<Cell>(bary, P, g, terms, c.data(), ct.data(), gradU, gradA);
    fields.gradU = gradU;
    fields.gradA = gradA;
  }
}

template void SimplexFreqInterpol<Triangle>::wall_at(const int, const std::array<double, 4>&, const FreqWeights&,
                                                     const double*, const std::size_t, PointValues&) const;
template void SimplexFreqInterpol<Tet>::wall_at(const int, const std::array<double, 4>&, const FreqWeights&,
                                                const double*, const std::size_t, PointValues&) const;
template void SimplexFreqInterpol<Triangle>::upgrade_at(const int, const std::array<double, 4>&, const double*,
                                                        const double*, PointValues&) const;
template void SimplexFreqInterpol<Tet>::upgrade_at(const int, const std::array<double, 4>&, const double*,
                                                   const double*, PointValues&) const;
