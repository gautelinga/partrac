#ifdef USE_DOLFIN
#ifndef __SIMPLEXINTERPOL_HPP
#define __SIMPLEXINTERPOL_HPP

// Velocity and pressure stamps written as dolfin HDF5 checkpoints, on
// triangles or tets: a P1 or P2 Taylor-Hood pair, evaluated through this
// code's own cell walk and basis rather than dolfin's.

#include "MeshInterpol.hpp"
#include "strings.hpp"
#include "Timestamps.hpp"
#include "Triangle.hpp"
#include "Tet.hpp"
#include "cell_locate.hpp"
#include <omp.h>

template<typename Cell>
class SimplexInterpol final
  : public MeshInterpol<Cell>
{
public:
  SimplexInterpol(const std::string& infilename);
  ~SimplexInterpol() { std::cout << "Destructing SimplexInterpol (" << mode << ")." << std::endl; };
  void update(const double t);
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& );
  double get_t_min() { return ts.get_t_min(); };
  double get_t_max() { return ts.get_t_max(); };
protected:
  static constexpr int D = Cell::n_verts - 1;
  static constexpr const char* mode = D == 2 ? "triangle" : "tet";
  using Base = MeshInterpol<Cell>;
  // Names the base owns
  using Base::dolfin_params; using Base::periodic; using Base::include_pressure;
  using Base::mesh; using Base::dim; using Base::u_space_; using Base::p_space_;
  using Base::ncoeffs_u; using Base::ncoeffs_p;
  using Base::cells_; using Base::dolfin_cells_;
  using Base::u_dofs_; using Base::p_dofs_; using Base::t_prev; using Base::t_next;
  using Base::x_min; using Base::x_max; using Base::is_initialized; using Base::t_update;
  using Base::get_folder; using Base::set_folder; using Base::wants_gradient;
  using Base::read_mesh_params; using Base::init_mesh_geometry; using Base::build_cells;

  Timestamps ts;

  std::shared_ptr<dolfin::Function> u_prev_;
  std::shared_ptr<dolfin::Function> u_next_;
  std::shared_ptr<dolfin::Function> p_prev_;
  std::shared_ptr<dolfin::Function> p_next_;

  std::vector<double> u_prev_data_;
  std::vector<double> u_next_data_;
  std::vector<double> p_prev_data_;
  std::vector<double> p_next_data_;
};

#endif
#endif
