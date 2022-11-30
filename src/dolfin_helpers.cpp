#include "dolfin_helpers.hpp"
#ifdef USE_DOLFIN

void build_neighbor_list( std::vector<std::set<Uint>> &cell2cells_
                        , std::shared_ptr<dolfin::Mesh> mesh
                        , std::vector<dolfin::Cell> &dolfin_cells_)
{
  Uint dim = mesh->geometry().dim();
  for (std::size_t i = 0; i < mesh->num_cells(); ++i){
    //std::cout << "Num cells:  " << dolfin_cells_[i].num_entities(dim) << std::endl;
    //std::cout << "Num facets: " << dolfin_cells_[i].num_entities(dim-1) << std::endl;
    
    auto facets = dolfin_cells_[i].entities(dim-1);
    for ( std::size_t j = 0; j < dolfin_cells_[i].num_entities(dim-1); ++j ){
      //std::cout << "Index: " << facets[j] << std::endl;
      dolfin::Facet dolfin_facet(*mesh, facets[j]);
      auto neighbor_cells = dolfin_facet.entities(dim);
      for (std::size_t k = 0; k < dolfin_facet.num_entities(dim); ++k){
        //std::cout << "Neigh: " << neighbor_cells[k] << std::endl;
        if (i != neighbor_cells[k]){
          cell2cells_[i].insert(neighbor_cells[k]);
        }
      }
    }
  }
  // TODO: Include periodic neighbor cells
}

std::tuple<Uint, bool> locate_cell( int &id_prev
                                  , dolfin::Array<double>& x_loc
                                  , std::shared_ptr<dolfin::Mesh> mesh
                                  , std::vector<std::set<Uint>>& cell2cells_){
  // Index of cell containing point
  Uint dim = mesh->geometry().dim();
  const dolfin::Point point(dim, x_loc.data());

  bool found = false;

  unsigned int id = 0;
  bool inside = false;
  // Search in neighborhood first
  if (id_prev >= 0){
    dolfin::Cell prev_cell(*mesh, id_prev);
    if (prev_cell.contains(point)){
      id = id_prev;
      inside = true;
      found = true;
    }
    else {
      for ( auto neigh_id : cell2cells_[id_prev]){
        dolfin::Cell neigh_cell(*mesh, neigh_id);
        if (neigh_cell.contains(point)){
          inside = true;
          found = true;
          id = neigh_id;
          break;
        }
      }
    }
  }
  if (!found){
    id = mesh->bounding_box_tree()->compute_first_entity_collision(point);
    inside = (id != std::numeric_limits<unsigned int>::max());
    if (inside) found = true;
  }
  if (found){
    id_prev = id;
  }
  return {id, inside};
}

#endif