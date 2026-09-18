#ifdef USE_DOLFIN
#include "SimplexInterpol.hpp"
#include "loader_params.hpp"
#include "dolfin_spaces.hpp"
#include "p12_eval.hpp"
#include "geometry.hpp"
#include "PeriodicBC.hpp"
#include "dolfin_helpers.hpp"
#include <array>
#include <numeric>
#include <cassert>

template<typename Cell>
SimplexInterpol<Cell>::SimplexInterpol(const std::string& infilename)
  : MeshInterpol<Cell>(infilename)
{
  dolfin_params = partrac::parse_file_or_exit(dolfin_h5_schema(mode), infilename);

  std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  ts.initialize(get_folder() + "/" + dolfin_params.template get<std::string>("timestamps"));

  read_mesh_params();

  std::string meshfilename = get_folder() + "/" + dolfin_params.template get<std::string>("mesh");
  dolfin::HDF5File meshfile(MPI_COMM_WORLD, meshfilename, "r");

  dolfin::Mesh mesh_in;
  meshfile.read(mesh_in, "mesh", false);

  mesh = std::make_shared<dolfin::Mesh>(mesh_in);
  init_mesh_geometry();

  auto constrained_domain = std::make_shared<PeriodicBC>(periodic, x_min, x_max, dim);

  std::string u_el = dolfin_params.template get<std::string>("velocity_space");
  std::string p_el = dolfin_params.template get<std::string>("pressure_space");

  taylor_hood_spaces<Cell>(u_el, p_el, include_pressure, mesh, constrained_domain,
                           u_space_, p_space_, ncoeffs_u, ncoeffs_p);

  build_cells(*u_space_->dofmap());

  u_prev_ = std::make_shared<dolfin::Function>(u_space_);
  u_next_ = std::make_shared<dolfin::Function>(u_space_);

  if (include_pressure){
    p_prev_ = std::make_shared<dolfin::Function>(p_space_);
    p_next_ = std::make_shared<dolfin::Function>(p_space_);
  }

  check_dofs_fit(ncoeffs_u, ncoeffs_p, Cell::n_dofs_max, "SimplexInterpol");

  // Precompute dofs of all cells
  u_dofs_.build(*u_space_->dofmap(), dolfin_cells_, "SimplexInterpol");
  u_dofs_.check_stride(D*ncoeffs_u, "SimplexInterpol");
  if (include_pressure){
    p_dofs_.build(*p_space_->dofmap(), dolfin_cells_, "SimplexInterpol");
    p_dofs_.check_stride(ncoeffs_p, "SimplexInterpol");
  }

  std::cout << "Setting max threads: " << omp_get_max_threads() << std::endl;
}

template<typename Cell>
void SimplexInterpol<Cell>::update(const double t)
{
  StampPair sp = ts.get(t);

  // Always load once; keep last bracket past t_max
  if ( !is_initialized || ((t_prev != sp.prev.t || t_next != sp.next.t) && t < ts.get_t_max()) ){
    const std::string u_field = dolfin_params.template get<std::string>("velocity_field");
    const std::string p_field = dolfin_params.template get<std::string>("pressure_field");

    // Swap if possible
    if (is_initialized && t_next == sp.prev.t){
      std::cout << "Prev: Timestep = " << sp.prev.t << ", swapping... " << std::endl;
      u_prev_data_.swap(u_next_data_);
      if (include_pressure)
        p_prev_data_.swap(p_next_data_);
    }
    else {
      std::cout << "Prev: Timestep = " << sp.prev.t << ", filename = " << sp.prev.filename << std::endl;
      dolfin::HDF5File prevfile(MPI_COMM_WORLD, get_folder() + "/" + sp.prev.filename, "r");
      prevfile.read(*u_prev_, u_field);
      u_prev_->vector()->get_local(u_prev_data_);
      if (include_pressure){
        prevfile.read(*p_prev_, p_field);
        p_prev_->vector()->get_local(p_prev_data_);
      }
    }

    std::cout << "Next: Timestep = " << sp.next.t << ", filename = " << sp.next.filename << std::endl;
    // Single stamp: copy prev
    if (sp.next.filename == sp.prev.filename){
      u_next_data_ = u_prev_data_;
      if (include_pressure)
        p_next_data_ = p_prev_data_;
    }
    else {
      dolfin::HDF5File nextfile(MPI_COMM_WORLD, get_folder() + "/" + sp.next.filename, "r");
      nextfile.read(*u_next_, u_field);
      u_next_->vector()->get_local(u_next_data_);
      if (include_pressure){
        nextfile.read(*p_next_, p_field);
        p_next_->vector()->get_local(p_next_data_);
      }
    }

    is_initialized = true;
    t_prev = sp.prev.t;
    t_next = sp.next.t;
  }
  t_update = t;
}

template<typename Cell>
void SimplexInterpol<Cell>::evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  evaluate_impl<true>(x, t, pos, fields);
}

template<typename Cell>
void SimplexInterpol<Cell>::evaluate_motion(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  evaluate_impl<false>(x, t, pos, fields);
}

template<typename Cell>
template<bool Scalars>
void SimplexInterpol<Cell>::evaluate_impl(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields)
{
  // Assuming inside fluid
  assert(t <= t_next && t >= t_prev);
  const double alpha_t = stamp_weight(t, t_prev, t_next);

  // Compute Pk-Pl basis at x
  const int id = pos.id;
  std::array<double, Cell::n_dofs_max> _Nu_, _Np_, _Nux_, _Nuy_, _Nuz_;   // _Nuz_ unused in 2D

  cell_basis(cells_[id], pos.bary, ncoeffs_u, _Nu_.data(), "u");
  if constexpr (Scalars)
    if (include_pressure)
      cell_basis(cells_[id], pos.bary, ncoeffs_p, _Np_.data(), "p");

  // Gathered: restrict() is not thread-safe
  std::array<double, Cell::n_dofs_max*3> u_prev_block, u_next_block;
  gather_stamps<D*Cell::n_verts, D*Cell::n_dofs_max>(u_dofs_[id], u_dofs_.stride(), u_prev_data_, u_next_data_,
                u_prev_block.data(), u_next_block.data());

  // Evaluate
  const Vector3d U_prev = block_value<D>(_Nu_.data(), u_prev_block.data(), ncoeffs_u);
  const Vector3d U_next = block_value<D>(_Nu_.data(), u_next_block.data(), ncoeffs_u);

  // Update
  fields.U = alpha_t * U_next + (1-alpha_t) * U_prev;
  fields.A = stamp_rate(U_next, U_prev, t_prev, t_next);

  if constexpr (Scalars) if (include_pressure){
    std::array<double, Cell::n_dofs_max> p_prev_block, p_next_block;
    gather_stamps<Cell::n_verts, Cell::n_dofs_max>(p_dofs_[id], p_dofs_.stride(), p_prev_data_, p_next_data_,
                  p_prev_block.data(), p_next_block.data());
    const double P_prev = block_scalar(_Np_.data(), p_prev_block.data(), ncoeffs_p);
    const double P_next = block_scalar(_Np_.data(), p_next_block.data(), ncoeffs_p);
    fields.P = alpha_t * P_next + (1-alpha_t) * P_prev;
  }

  if (wants_gradient()){
    cell_deriv(cells_[id], pos.bary, ncoeffs_u, _Nux_.data(), _Nuy_.data(), _Nuz_.data(), "u");
    const Matrix3d gradU_prev = block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(), u_prev_block.data(), ncoeffs_u);
    const Matrix3d gradU_next = block_gradient<D>(_Nux_.data(), _Nuy_.data(), _Nuz_.data(), u_next_block.data(), ncoeffs_u);
    fields.gradU = alpha_t * gradU_next + (1-alpha_t) * gradU_prev;
    fields.gradA = stamp_rate(gradU_next, gradU_prev, t_prev, t_next);
  }
}

template class SimplexInterpol<Triangle>;
template class SimplexInterpol<Tet>;

#endif
