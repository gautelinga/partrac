#include "SimplexFreqInterpol.hpp"
#include "loader_params.hpp"
#include "mesh_tables.hpp"
#include "p12_eval.hpp"
#include "phase_timing.hpp"
#include "simplex_load.hpp"
#include <array>
#include <atomic>
#include <cmath>
#include <iostream>

namespace {

// The components are one function space written many times: a later file that
// holds another element is not summable with the first
void same_element(const simplex_load::Element& a, const simplex_load::Element& b,
                  const std::string& path, const std::string& field){
  if (a.degree != b.degree || a.ncomp != b.ncomp)
    partrac::fail(path, ": '", field, "' is of degree ", b.degree, " with ", b.ncomp,
                  " components, the first component's degree ", a.degree, " with ", a.ncomp);
}

// Times a thread keeps the weights of: an RK4 step asks three distinct ones
constexpr int weight_slots = 4;

// Serial of each loader built, so kept weights are never read as a later
// loader's at the same address
std::atomic<std::uint64_t> next_id{1};

// The weights last asked, per thread; plain data, read without a guard
thread_local FreqWeights kept[weight_slots];
thread_local int next_kept = 0;

}  // namespace

template<typename Cell>
SimplexFreqInterpol<Cell>::SimplexFreqInterpol(const std::string& infilename)
  : MeshCore<Cell>(infilename)
{
  constexpr int nv = Cell::n_verts;
  id_ = next_id++;
  partrac::phase_begin("load");
  dolfin_params = partrac::parse_file_or_exit(simplex_freq_schema(mode), infilename);

  const std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  // Using the FreqStamps class to hold frequency data
  fs.initialize(get_folder() + "/" + dolfin_params.template get<std::string>("freqstamps"));
  if (fs.size() < 1)
    partrac::fail(infilename, ": the frequency stamps list no component");

  // base frequency, for lines that give a harmonic by their number
  if (fs.omega_given() == dolfin_params.has("tau"))
    partrac::fail(infilename, fs.omega_given() ? ": tau is not read when the frequency stamps give omega"
                                               : ": tau is needed when the frequency stamps give no omega");
  if (!fs.omega_given()){
    const double tau = dolfin_params.template get<double>("tau");
    if (tau > 0)
      omega0 = 2 * M_PI / tau;
  }

  read_mesh_params();
  partrac::phase("params");

  simplex_load::Request req;
  req.infilename = infilename;
  // The element comes from the first component
  req.field_file = get_folder() + "/" + fs.get(0).filename;
  req.what = "SimplexFreqInterpol";
  req.nv = nv;
  req.include_pressure = include_pressure;
  req.n_dofs_max = Cell::n_dofs_max;
  req.periodic = periodic;
  req.periodic_tol = periodic_tol;
  simplex_load::request_from_params(req, dolfin_params, get_folder());
  const std::string& u_field = req.u_field;
  const std::string& p_field = req.p_field;
  simplex_load::Tables t;
  simplex_load::build_tables(req, t);

  this->adopt_tables(t);
  ncells_ = t.mesh.ncells;
  nverts_ = t.mesh.nverts;

  // Every component through one dof -> nodes mapping: the components are one
  // function space, so only the first has to be read with its dof table
  const std::size_t nfreq = std::size_t(fs.size());
  u_nodes_.resize(nfreq);
  if (include_pressure)
    p_nodes_.resize(nfreq);
  mesh_tables::DofNodes u_map, p_map;
  for (std::size_t iFreq = 0; iFreq < nfreq; ++iFreq){
    const std::string path = get_folder() + "/" + fs.get(int(iFreq)).filename;
    if (iFreq == 0){
      simplex_load::read_field_by_node(path, u_field, t, t.el_u, t.node_order(t.el_u),
                                       u_nodes_[0], u_map);
      if (include_pressure)
        simplex_load::read_field_by_node(path, p_field, t, t.el_p, t.node_order(t.el_p),
                                         p_nodes_[0], p_map);
    }
    else {
      same_element(t.el_u, simplex_load::read_element(path, u_field, nv), path, u_field);
      simplex_load::read_into(path, u_field, u_map, u_nodes_[0], u_nodes_[iFreq]);
      if (include_pressure){
        same_element(t.el_p, simplex_load::read_element(path, p_field, nv), path, p_field);
        simplex_load::read_into(path, p_field, p_map, p_nodes_[0], p_nodes_[iFreq]);
      }
    }
  }
  partrac::phase("component values");

  topo_ = std::move(t.mesh.topo);
  coords_ = std::move(t.mesh.coords);
  mesh_tables::build_cells_and_tree<Cell>(topo_, coords_, ncells_, nverts_, dim, cells_, tree_);
  partrac::phase_total();

  std::cout << "Setting max threads: " << omp_get_max_threads() << std::endl;
}

template<typename Cell>
void SimplexFreqInterpol<Cell>::update(const double t)
{
  is_initialized = true;
  t_update = t;
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

  // Compute Pk-Pl basis at x
  std::array<double, Cell::n_dofs_max> Nu_, Np_, Nux_, Nuy_, Nuz_;   // Nuz_ unused in 2D

  cell_basis(cells_[id], pos.bary, ncoeffs_u, Nu_.data(), "u");
  if constexpr (Scalars)
    if (include_pressure)
      cell_basis(cells_[id], pos.bary, ncoeffs_p, Np_.data(), "p");

  // Gathered by node: the D components of a node are consecutive
  const std::uint32_t* u_row = u_dofs_[id];
  std::array<double, Cell::n_dofs_max*3> u_block;

  const bool gradient = wants_gradient();
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
    const Vector3d U_f = block_value<D>(Nu_.data(), u_block.data(), ncoeffs_u);
    U += wf.w[iFreq]*U_f;
    A += wf.wt[iFreq]*U_f;
    if (gradient){
      const Matrix3d G_f = block_gradient<D>(Nux_.data(), Nuy_.data(), Nuz_.data(),
                                             u_block.data(), ncoeffs_u);
      gradU += wf.w[iFreq]*G_f;
      gradA += wf.wt[iFreq]*G_f;
    }
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

template class SimplexFreqInterpol<Triangle>;
template class SimplexFreqInterpol<Tet>;
