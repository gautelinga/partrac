#include "SimplexFreqInterpol.hpp"
#include "p12_eval.hpp"
#include <array>
#include <atomic>
#include <cmath>

// SimplexFreqInterpol's evaluation, located and held, in a unit of its own;
// a wall cell's is out of line in another

namespace {

// Times a thread keeps the weights of: an RK4 step asks three distinct ones
constexpr int weight_slots = 4;

// One component's share of the velocity, its rate and their gradients at a point
template<int D>
__attribute__((always_inline))
inline void add_component(const double* N, const double* Nx, const double* Ny, const double* Nz,
                          const double* block, const Uint ncoeffs, const double w, const double wt,
                          Vector3d& U, Vector3d& A, Matrix3d& gradU, Matrix3d& gradA, const bool gradient){
  const Vector3d U_f = block_value<D>(N, block, ncoeffs);
  U += w*U_f;
  A += wt*U_f;
  if (gradient){
    const Matrix3d G_f = block_gradient<D>(Nx, Ny, Nz, block, ncoeffs);
    gradU += w*G_f;
    gradA += wt*G_f;
  }
}

// Serial of each loader built, so kept weights are never read as a later
// loader's at the same address
std::atomic<std::uint64_t> next_id{1};

// The weights last asked, per thread; plain data, read without a guard
thread_local FreqWeights kept[weight_slots];
thread_local int next_kept = 0;

}  // namespace

template<typename Cell>
void SimplexFreqInterpol<Cell>::new_id()
{
  id_ = next_id++;
}

template<typename Cell>
inline const FreqWeights& SimplexFreqInterpol<Cell>::weights(const double t)
{
  for (const FreqWeights& e : kept)
    if (e.owner == id_ && e.t == t)
      return e;
  return fill_weights(t);
}

template<typename Cell>
const FreqWeights& SimplexFreqInterpol<Cell>::fill_weights(const double t)
{
  // The values behind each slot: the cosines, then the rates
  static thread_local std::array<std::vector<double>, weight_slots> store;
  const int i = next_kept;
  next_kept = (next_kept + 1) % weight_slots;
  const std::size_t n = std::size_t(fs.size());
  std::vector<double>& v = store[std::size_t(i)];
  v.resize(2*n);
  double* w = v.data();
  double* wt = v.data() + n;
  for (std::size_t iFreq = 0; iFreq < n; ++iFreq){
    FreqStamp& f = fs.get(int(iFreq));
    double a = f.a;
    double t_shift = f.t;
    if (fs.omega_given()){
      w[iFreq] = a * cos(f.omega * t + t_shift);
      wt[iFreq] = - a * f.omega * sin(f.omega * t + t_shift);
    }
    else {
      w[iFreq] = a * cos(omega0 * (iFreq * t + t_shift));
      wt[iFreq] = - a * omega0 * iFreq * sin(omega0 * (iFreq * t + t_shift));
    }
  }
  kept[i] = {id_, t, w, wt};
  return kept[i];
}

template<typename Cell>
void SimplexFreqInterpol<Cell>::evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  evaluate_impl<true>(x, t, pos, fields);
}

template<typename Cell>
void SimplexFreqInterpol<Cell>::evaluate_motion(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  evaluate_impl<false>(x, t, pos, fields);
}

template<typename Cell>
template<bool Scalars>
void SimplexFreqInterpol<Cell>::evaluate_impl(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  // Assume found in fluid domain
  const int id = pos.id;
  const FreqWeights& wf = weights(t);
  const std::size_t nfreq = u_nodes_.size();

  // Quadratic near walls at rest
  if (wall_upgrade_ && wall_index_[id] >= 0){
    wall_at(id, pos.bary, wf, nullptr, 0, fields);
    if constexpr (Scalars)
      if (include_pressure) fields.P = pressure_at(id, pos.bary, wf);
    return;
  }

  const bool gradient = wants_gradient();
  if constexpr (!Scalars)
    if (!gradient){
      velocity_at(id, pos.bary, wf, fields);
      return;
    }

  // Compute Pk-Pl basis at x
  std::array<double, Cell::n_dofs_max> Nu_, Np_, Nux_, Nuy_, Nuz_;   // Nuz_ unused in 2D

  cell_basis(cells_[id], pos.bary, ncoeffs_u, Nu_.data(), "u");
  if constexpr (Scalars)
    if (include_pressure)
      cell_basis(cells_[id], pos.bary, ncoeffs_p, Np_.data(), "p");

  // Gathered by node: the D components of a node are consecutive
  const std::uint32_t* u_row = u_dofs_[id];
  std::array<double, Cell::n_dofs_max*3> u_block;

  if (gradient)
    cell_deriv(cells_[id], pos.bary, ncoeffs_u, Nux_.data(), Nuy_.data(), Nuz_.data(), "u");

  Vector3d U = Vector3d::Zero();
  Vector3d A = Vector3d::Zero();
  Matrix3d gradU = Matrix3d::Zero();
  Matrix3d gradA = Matrix3d::Zero();
  double P = 0.;
  for (std::size_t iFreq = 0; iFreq < nfreq; ++iFreq){
    gather_cell_nodes<Cell::n_verts, Cell::n_dofs_max, D>(u_row, u_dofs_.stride(),
                      u_nodes_[iFreq].data(), u_block.data());
    add_component<D>(Nu_.data(), Nux_.data(), Nuy_.data(), Nuz_.data(), u_block.data(), ncoeffs_u,
                     wf.w[iFreq], wf.wt[iFreq], U, A, gradU, gradA, gradient);
    if constexpr (Scalars){
      if (include_pressure){
        std::array<double, Cell::n_dofs_max> p_block;
        gather_cell_nodes<Cell::n_verts, Cell::n_dofs_max, 1>(p_dofs_[id], p_dofs_.stride(),
                          p_nodes_[iFreq].data(), p_block.data());
        P += wf.w[iFreq]*block_scalar(Np_.data(), p_block.data(), ncoeffs_p);
      }
    }
  }

  // Update
  fields.U = U;
  fields.A = A;
  if constexpr (Scalars){
    if (include_pressure)
      fields.P = P;
  }
  if (gradient){
    fields.gradU = gradU;
    fields.gradA = gradA;
  }
}

template<typename Cell>
void SimplexFreqInterpol<Cell>::velocity_at(const int id, const std::array<double, 4>& bary,
                                            const FreqWeights& wf, PointValues& fields) const
{
  std::array<double, Cell::n_dofs_max> Nu_;
  std::array<double, Cell::n_dofs_max*3> u_block;
  cell_basis(cells_[id], bary, ncoeffs_u, Nu_.data(), "u");
  Vector3d U = Vector3d::Zero();
  Vector3d A = Vector3d::Zero();
  for (std::size_t iFreq = 0; iFreq < u_nodes_.size(); ++iFreq){
    gather_cell_nodes<Cell::n_verts, Cell::n_dofs_max, D>(u_dofs_[id], u_dofs_.stride(),
                      u_nodes_[iFreq].data(), u_block.data());
    const Vector3d U_f = block_value<D>(Nu_.data(), u_block.data(), ncoeffs_u);
    U += wf.w[iFreq]*U_f;
    A += wf.wt[iFreq]*U_f;
  }
  fields.U = U;
  fields.A = A;
}

template<typename Cell>
inline double SimplexFreqInterpol<Cell>::pressure_at(const int id, const std::array<double, 4>& bary,
                                                     const FreqWeights& wf) const
{
  std::array<double, Cell::n_dofs_max> Np_, p_block;
  cell_basis(cells_[id], bary, ncoeffs_p, Np_.data(), "p");
  double P = 0.;
  for (std::size_t iFreq = 0; iFreq < p_nodes_.size(); ++iFreq){
    gather_cell_nodes<Cell::n_verts, Cell::n_dofs_max, 1>(p_dofs_[id], p_dofs_.stride(),
                      p_nodes_[iFreq].data(), p_block.data());
    P += wf.w[iFreq]*block_scalar(Np_.data(), p_block.data(), ncoeffs_p);
  }
  return P;
}

template<typename Cell>
void SimplexFreqInterpol<Cell>::hold(const Region& R, Held& h) const
{
  const int id = R.id;
  constexpr std::size_t bs = Cell::n_dofs_max*3;
  const std::size_t nfreq = u_nodes_.size();
  h.blocks.resize(nfreq*bs);
  h.wall = wall_upgrade_ && wall_index_[id] >= 0;
  for (std::size_t iFreq = 0; iFreq < nfreq; ++iFreq)
    gather_cell_nodes<Cell::n_verts, Cell::n_dofs_max, D>(u_dofs_[id], u_dofs_.stride(),
                      u_nodes_[iFreq].data(), h.blocks.data() + iFreq*bs);
}

template<typename Cell>
void SimplexFreqInterpol<Cell>::held_motion(const int id, const std::array<double, 4>& lev,
                                            const double t, const Held& h, PointValues& fields)
{
  constexpr std::size_t bs = Cell::n_dofs_max*3;
  const FreqWeights& wf = weights(t);
  if (h.wall){
    wall_at(id, lev, wf, h.blocks.data(), bs, fields);
    return;
  }
  const std::size_t nfreq = u_nodes_.size();
  std::array<double, Cell::n_dofs_max> Nu_, Nux_, Nuy_, Nuz_;   // Nuz_ unused in 2D
  cell_basis(cells_[id], lev, ncoeffs_u, Nu_.data(), "u");
  const bool gradient = wants_gradient();
  if (gradient)
    cell_deriv(cells_[id], lev, ncoeffs_u, Nux_.data(), Nuy_.data(), Nuz_.data(), "u");

  Vector3d U = Vector3d::Zero();
  Vector3d A = Vector3d::Zero();
  Matrix3d gradU = Matrix3d::Zero();
  Matrix3d gradA = Matrix3d::Zero();
  for (std::size_t iFreq = 0; iFreq < nfreq; ++iFreq)
    add_component<D>(Nu_.data(), Nux_.data(), Nuy_.data(), Nuz_.data(), h.blocks.data() + iFreq*bs,
                     ncoeffs_u, wf.w[iFreq], wf.wt[iFreq], U, A, gradU, gradA, gradient);
  fields.U = U;
  fields.A = A;
  if (gradient){
    fields.gradU = gradU;
    fields.gradA = gradA;
  }
}

#define FREQ_EVAL_INSTANCES(Cell) \
  template void SimplexFreqInterpol<Cell>::new_id(); \
  template const FreqWeights& SimplexFreqInterpol<Cell>::fill_weights(const double); \
  template void SimplexFreqInterpol<Cell>::evaluate(const Vector3d&, const double, const CellPos&, PointValues&); \
  template void SimplexFreqInterpol<Cell>::evaluate_motion(const Vector3d&, const double, const CellPos&, PointValues&); \
  template void SimplexFreqInterpol<Cell>::hold(const Region&, Held&) const; \
  template void SimplexFreqInterpol<Cell>::held_motion(const int, const std::array<double, 4>&, const double, \
                                                       const Held&, PointValues&);

FREQ_EVAL_INSTANCES(Triangle)
FREQ_EVAL_INSTANCES(Tet)
