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

void label_cell_type(std::vector<int>& cell_type_, std::vector<std::set<Uint>>& cell2cells_, const Uint dim){
  //cell_type_.clear(); // set all to zero
  for ( Uint i=0; i < cell2cells_.size(); ++i)
  {
    cell_type_[i] = 0;
  }

  for ( Uint i=0; i < cell2cells_.size(); ++i )
  {
    auto & cells_loc = cell2cells_[i];
    if (cells_loc.size() < dim+1){
      cell_type_[i] = 1;

      for ( auto & cell_loc : cells_loc ){
        if (cell_type_[cell_loc] == 0) cell_type_[cell_loc] = 2;
      }
    }
  }
}

void apply_periodic_boundaries(std::vector<std::set<Uint>>& cell2cells_,
                               //std::vector<int>& cell_type_,
                               const std::vector<bool>& periodic,
                               const Vector3d& x_min,
                               const Vector3d& x_max,
                               std::shared_ptr<dolfin::Mesh> mesh,
                               const std::vector<dolfin::Cell> &dolfin_cells_,
                               const Uint dim,
                               const double tol)
{
  // Needs to be generalized for 3D use
  std::vector<std::vector<std::pair<Vector3d, Uint>>> bdry_l(dim);
  std::vector<std::vector<std::pair<Vector3d, Uint>>> bdry_h(dim);

  for ( Uint i=0; i < mesh->num_cells(); ++i ){
    //if ( cell_type_[i] == 1 ){
      auto facets = dolfin_cells_[i].entities(dim-1);
      for ( std::size_t j = 0; j < dolfin_cells_[i].num_entities(dim-1); ++j ){
        dolfin::Facet dolfin_facet(*mesh, facets[j]);

        if (dolfin_facet.exterior()){
          Vector3d pt(dolfin_facet.midpoint().coordinates());

          for ( Uint k=0; k < dim; ++k)
          {
            if (periodic[k]){

              if (pt[k] < x_min[k] + tol)
              {
                //std::cout << pt << std::endl;
                pt[k] += x_max[k] - x_min[k];
                bdry_l[k].push_back({pt, i});
              }
              else if (pt[k] > x_max[k] - tol)
              {
                bdry_h[k].push_back({pt, i});
              }
            }
          }
        }
      }
    //}
  }

  for ( Uint k=0; k < dim; ++k ){
    std::cout << k << " " << bdry_l[k].size() << " " << bdry_h[k].size() << std::endl;
  }

  for ( Uint k=0; k < dim; ++k ){
    for ( auto & item1 : bdry_l[k] ){
      auto x1 = item1.first;
      auto id1 = item1.second;

      for ( auto & item2 : bdry_h[k] ){
        auto x2 = item2.first;
        auto id2 = item2.second;

        double dx = (x1-x2).norm();

        if ( dx < tol ){
          //std::cout << dx << " " << id1 << " " << id2 << std::endl;
          cell2cells_[id1].insert(id2);
          cell2cells_[id2].insert(id1);
        }
      }
    }
  }
}

#endif