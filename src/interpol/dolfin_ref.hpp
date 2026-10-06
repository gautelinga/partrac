#ifdef USE_DOLFIN
#ifndef __DOLFIN_REF_HPP
#define __DOLFIN_REF_HPP

// What the reference loader (DolfInterpol, mode fenics) takes from dolfin
// besides its own basis: the periodic subdomain, the Lagrange spaces, the
// per-cell dof table and the cell order a dofmap implies, the facet table off
// the mesh entities, and the bounding box tree. The tests that write their
// fixtures with dolfin use the same pieces. Nothing else in src/ does.

#include <dolfin.h>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include "typedefs.hpp"
#include "PointValues.hpp"
#include "cell_walk.hpp"

// The min face of each periodic direction, and the map from its image
class PeriodicBC : public dolfin::SubDomain {
public:
  PeriodicBC(const std::vector<bool> &periodic,
      const Vector3d &x_min, const Vector3d &x_max, const Uint dim) {
    periodic_x = periodic[0] && dim > 0;
    periodic_y = periodic[1] && dim > 1;
    periodic_z = periodic[2] && dim > 2;
    this->x_min = x_min;
    this->x_max = x_max;
    this->dim = dim;
  };
  bool inside(const dolfin::Array<double>& xx, bool on_boundary) const {
    std::vector<double> x = {0., 0., 0.};
    for (std::size_t d=0; d<dim; ++d){
      x[d] = xx[d];
    }
    return (on_boundary &&
            ((periodic_x && x[0] < x_min[0] + DOLFIN_EPS_LARGE) ||
             (periodic_y && x[1] < x_min[1] + DOLFIN_EPS_LARGE) ||
             (periodic_z && x[2] < x_min[2] + DOLFIN_EPS_LARGE)) &&
            !(
              (periodic_x && periodic_y &&
               x[0] < x_min[0] + DOLFIN_EPS_LARGE &&
               x[1] > x_max[1] - DOLFIN_EPS_LARGE) ||
              (periodic_x && periodic_y &&
               x[0] > x_max[0] - DOLFIN_EPS_LARGE &&
               x[1] < x_min[1] + DOLFIN_EPS_LARGE) ||
              (periodic_x && periodic_z &&
               x[0] < x_min[0] + DOLFIN_EPS_LARGE &&
               x[2] > x_max[2] - DOLFIN_EPS_LARGE) ||
              (periodic_x && periodic_z &&
               x[0] > x_max[0] - DOLFIN_EPS_LARGE &&
               x[2] < x_min[2] + DOLFIN_EPS_LARGE) ||
              (periodic_y && periodic_z &&
               x[1] < x_min[1] + DOLFIN_EPS_LARGE &&
               x[2] > x_max[2] - DOLFIN_EPS_LARGE) ||
              (periodic_y && periodic_z &&
               x[1] > x_max[1] - DOLFIN_EPS_LARGE &&
               x[2] < x_min[2] + DOLFIN_EPS_LARGE)
              )
            );
  }
  void map(const dolfin::Array<double>& xx, dolfin::Array<double>& yy) const {
    std::vector<double> x = {0., 0., 0.};
    std::vector<double> y = {0., 0., 0.};
    for (std::size_t d=0; d<dim; ++d){
      x[d] = xx[d];
    }
    if (periodic_x && periodic_y && periodic_z &&
        x[0] > x_max[0] - DOLFIN_EPS_LARGE &&
        x[1] > x_max[1] - DOLFIN_EPS_LARGE &&
        x[2] > x_max[2] - DOLFIN_EPS_LARGE){
      y[0] = x[0] - (x_max[0]-x_min[0]);
      y[1] = x[1] - (x_max[1]-x_min[1]);
      y[2] = x[2] - (x_max[2]-x_min[2]);
    }
    else if (periodic_x && periodic_y &&
             x[0] > x_max[0] - DOLFIN_EPS_LARGE &&
             x[1] > x_max[1] - DOLFIN_EPS_LARGE){
      y[0] = x[0] - (x_max[0]-x_min[0]);
      y[1] = x[1] - (x_max[1]-x_min[1]);
      y[2] = x[2];
    }
    else if (periodic_x && periodic_z &&
             x[0] > x_max[0] - DOLFIN_EPS_LARGE &&
             x[2] > x_max[2] - DOLFIN_EPS_LARGE){
      y[0] = x[0] - (x_max[0]-x_min[0]);
      y[1] = x[1];
      y[2] = x[2] - (x_max[2]-x_min[2]);
    }
    else if (periodic_y && periodic_z &&
             x[1] > x_max[1] - DOLFIN_EPS_LARGE &&
             x[2] > x_max[2] - DOLFIN_EPS_LARGE){
      y[0] = x[0];
      y[1] = x[1] - (x_max[1]-x_min[1]);
      y[2] = x[2] - (x_max[2]-x_min[2]);
    }
    else if (periodic_x && x[0] > x_max[0] - DOLFIN_EPS_LARGE){
      y[0] = x[0] - (x_max[0]-x_min[0]);
      y[1] = x[1];
      y[2] = x[2];
    }
    else if (periodic_y && x[1] > x_max[1] - DOLFIN_EPS_LARGE){
      y[0] = x[0];
      y[1] = x[1] - (x_max[1]-x_min[1]);
      y[2] = x[2];
    }
    else if (periodic_z && x[2] > x_max[2] - DOLFIN_EPS_LARGE){
      y[0] = x[0];
      y[1] = x[1];
      y[2] = x[2] - (x_max[2]-x_min[2]);
    }
    else {
      y[0] = x[0]-1000;
      y[1] = x[1]-1000;
      y[2] = x[2]-1000;
    }
    for (std::size_t d=0; d<dim; ++d){
      yy[d] = y[d];
    }
  }
private:
  double periodic_x;
  double periodic_y;
  double periodic_z;
  Vector3d x_min;
  Vector3d x_max;
  double dim;
};

// A Lagrange space P1-P3 by name, vector or scalar, of dimension D
template<int D, bool Vector>
std::shared_ptr<dolfin::FunctionSpace> lagrange_space(const std::string& el,
                                                      std::shared_ptr<dolfin::Mesh> mesh,
                                                      std::shared_ptr<const dolfin::SubDomain> cd,
                                                      const char* what);

// Every cell's dofs, sorted, flat, fixed stride
std::vector<int> sorted_dof_table(const dolfin::GenericDofMap& dofmap,
                                  const std::size_t ncells, std::size_t& stride);

// Fraction of consecutive cells sharing a dof; low in a poorly ordered mesh
double dof_sharing(const dolfin::GenericDofMap& dofmap, const std::size_t ncells);

// Cells in the order of their sorted dofs; returns map[old] -> new
std::vector<std::uint32_t> order_cells_by_dofs(const dolfin::GenericDofMap& dofmap,
                                               const std::size_t ncells);

// Cell order: dolfin's, unless consecutive cells rarely share a dof
std::vector<std::uint32_t> cell_order(const dolfin::GenericDofMap& dofmap,
                                      const std::size_t ncells,
                                      const std::string& mode_in,
                                      std::vector<std::uint32_t>& dolfin2local);

// Per cell, the neighbour across the facet facing each vertex (walk_to_cell, reflect_in_cells)
void build_facet_neighbours(std::vector<std::int32_t>& across,
                            std::shared_ptr<dolfin::Mesh> mesh,
                            const std::vector<dolfin::Cell>& dolfin_cells_,
                            const std::vector<std::uint32_t>* dolfin2local,
                            const std::vector<bool>& periodic,
                            const Vector3d& x_min,
                            const Vector3d& x_max,
                            const Uint dim,
                            const double tol);

// The bounding-box tree: every cell, from nothing known
template<typename Cell>
inline bool tree_to_cell(const std::vector<Cell>& cells,
                         const dolfin::Mesh& mesh,
                         const Uint dim,
                         const Vector3d& xx,
                         CellPos& pos,
                         const std::vector<std::uint32_t>* dolfin2local = nullptr,
                         FoundCounts* count = nullptr){
  const dolfin::Point point(dim, xx.data());
  const unsigned int id = mesh.bounding_box_tree()->compute_first_entity_collision(point);
  if (id == std::numeric_limits<unsigned int>::max())
    return false;
  if (count) ++count->tree;
  pos.id = dolfin2local ? int((*dolfin2local)[id]) : int(id);
  // Tree tolerance: may sit just outside
  cells[pos.id].contains(xx, pos.bary);
  return true;
}

#endif
#endif
