#ifndef __DOLFIN_HELPERS_HPP
#define __DOLFIN_HELPERS_HPP
#ifdef USE_DOLFIN


#include "typedefs.hpp"
#include <tuple>
#include <dolfin.h>

void build_neighbor_list( std::vector<std::set<Uint>> &cell2cells_
                        , std::shared_ptr<dolfin::Mesh> mesh
                        , std::vector<dolfin::Cell> &dolfin_cells_);

std::tuple<Uint, bool> locate_cell( int &id_prev
                                  , dolfin::Array<double>& x_loc
                                  , std::shared_ptr<dolfin::Mesh> mesh
                                  , std::vector<std::set<Uint>>& cell2cells_);

#endif
#endif