#ifndef __DOLFIN_HELPERS_HPP
#define __DOLFIN_HELPERS_HPP
#ifdef USE_DOLFIN


#include "typedefs.hpp"
#include "cell_locate.hpp"
#include <tuple>
#include <dolfin.h>

// 0 bulk, 1 on a wall, 2 next to a cell on a wall; from the facet table
void label_cell_type(std::vector<int>& cell_type_, const std::vector<std::int32_t>& across, const Uint nv);

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

// Box lengths on periodic axes, zero on the others
Vector3d periodic_lengths(const std::vector<bool>& periodic, const Vector3d& x_min,
                          const Vector3d& x_max, const Uint dim);
#endif
#endif