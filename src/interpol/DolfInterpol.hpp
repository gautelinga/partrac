#ifndef __DOLFINTERPOL_HPP
#define __DOLFINTERPOL_HPP

#include <memory>
#include "Interpol.hpp"
#include "Params.hpp"
#include "Timestamps.hpp"

//#include "H5Cpp.h"
//#define hid_t aa_hid_t
//#undef hid_t
#include <dolfin.h>
//#define hid_t bb_hid_t
//#include <dolfin/io/HDF5File.h>
//#include <dolfin/io/XDMFFile.h>
//#undef hid_t
//#define hid_t ambiguous use aa_hid_t or bb_hid_t

#include "dolfin_elements/vP1_2.h"
#include "dolfin_elements/vP2_2.h"
#include "dolfin_elements/vP3_2.h"
#include "dolfin_elements/vP1_3.h"
#include "dolfin_elements/vP2_3.h"
#include "dolfin_elements/vP3_3.h"
#include "dolfin_elements/P1_2.h"
#include "dolfin_elements/P2_2.h"
#include "dolfin_elements/P3_2.h"
#include "dolfin_elements/P1_3.h"
#include "dolfin_elements/P2_3.h"
#include "dolfin_elements/P3_3.h"
#include "PeriodicBC.hpp"
#include "Triangle.hpp"
#include "Tet.hpp"
#include "cell_locate.hpp"
#include "H5Cpp.h"
#include "strings.hpp"


using namespace H5;
//using namespace std;
//using namespace dolfin;


class DolfInterpol final
  : public Interpol {
public:
  DolfInterpol(const std::string& infilename);
  void update(const double t);
  using Interpol::locate;
  using Interpol::evaluate;
  bool locate(const Vector3d &x, const double t, CellPos& pos);
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& ptvals);
  // Walk a move off the walls
  bool reflect(const Vector3d& x, Vector3d& dx, CellPos& pos);
  void enable_reflection();
  double hmin() const { return mesh->hmin(); }
  double get_t_min() { return ts.get_t_min(); };
  double get_t_max() { return ts.get_t_max(); };
  double get_rho() {
    if (dolfin_params.has("rho"))
      return dolfin_params.get<double>("rho");
    else {
      std::cout << "dolfin_params does not contain \"rho\"" << std::endl;
      exit(1);
    }
  };
protected:
  Timestamps ts;
  double t_prev = 0.;
  double t_next = 0.;

  std::vector<bool> periodic = {false, false, false};
  //Vector3d x_min = {0., 0., 0.};
  //Vector3d x_max = {0., 0., 0.};

  Vector3d U = {0., 0., 0.};
  //double Uy = 0.;
  //double Uz = 0.;
  Vector3d A = {0., 0., 0.};
  double P = 0.;
  Matrix3d gradU, gradA;

  bool inside;

  partrac::Params dolfin_params;

  std::shared_ptr<dolfin::Mesh> mesh;
  Uint dim;

  std::shared_ptr<dolfin::FunctionSpace> u_space;
  std::shared_ptr<dolfin::FunctionSpace> p_space;

  std::shared_ptr<dolfin::Function> u_prev_;
  std::shared_ptr<dolfin::Function> u_next_;
  std::shared_ptr<dolfin::Function> p_prev_;
  std::shared_ptr<dolfin::Function> p_next_;
  // Per cell, once: what evaluate used to rebuild on every call
  std::vector<dolfin::Cell> dolfin_cells_;
  std::vector<int> cell_orientations_;   // empty when the mesh carries none
  std::vector<double> coordinate_dofs_;  // ncoords_ of them per cell, flat
  Uint ncoords_ = 0;
  std::vector<Triangle> triangles_;   // one of these two is filled, by dim
  std::vector<Tet> tets_;
  std::vector<CellNeighbours> cell2cells_;
  std::vector<std::int32_t> facet_neigh_;   // reflect_in_cells
  Vector3d period_ = Vector3d::Zero();
  std::vector<std::uint32_t> dolfin2local_;   // empty: dolfin's cell order
  std::shared_ptr<const dolfin::FiniteElement> u_element_, p_element_;
  Uint u_dim_ = 0, p_dim_ = 0;
  // Read out whole at each load: dolfin's vector is not safe to read in parallel
  std::vector<double> u_prev_data_, u_next_data_, p_prev_data_, p_next_data_;
  CellDofs u_dofs_, p_dofs_;
  std::vector<FoundCounts> found_;
  Vector3d _modx(const Vector3d&);


};

#endif
