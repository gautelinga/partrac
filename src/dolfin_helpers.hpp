#ifndef __DOLFIN_HELPERS_HPP
#define __DOLFIN_HELPERS_HPP
#ifdef USE_DOLFIN


#include "typedefs.hpp"
#include "cell_locate.hpp"
#include <tuple>
#include <dolfin.h>

void build_neighbor_list( std::vector<CellNeighbours> &cell2cells_
                        , std::shared_ptr<dolfin::Mesh> mesh
                        , std::vector<dolfin::Cell> &dolfin_cells_
                        , const std::vector<std::uint32_t>* dolfin2local = nullptr);

void label_cell_type(std::vector<int>& cell_type_, std::vector<CellNeighbours>& cell2cells_, const Uint dim);

void apply_periodic_boundaries(std::vector<CellNeighbours>& cell2cells_,
                               //std::vector<int>& cell_type_,
                               const std::vector<bool>& periodic,
                               const Vector3d& x_min,
                               const Vector3d& x_max,
                               std::shared_ptr<dolfin::Mesh> mesh,
                               const std::vector<dolfin::Cell> &dolfin_cells_,
                               const Uint dim,
                               const double tol);
#endif
#endif