#include "dolfin_helpers.hpp"
#ifdef USE_DOLFIN
#include <algorithm>
#include <cmath>

void build_neighbor_list( std::vector<CellNeighbours> &cell2cells_
                        , std::shared_ptr<dolfin::Mesh> mesh
                        , std::vector<dolfin::Cell> &dolfin_cells_
                        , const std::vector<std::uint32_t>* dolfin2local)
{
  // Cell ids must fit int
  if (mesh->num_cells() > std::size_t(std::numeric_limits<int>::max())){
    std::cout << "Mesh has " << mesh->num_cells() << " cells, more than a cell id can hold" << std::endl;
    exit(1);
  }
  Uint dim = mesh->geometry().dim();
  for (std::size_t i = 0; i < mesh->num_cells(); ++i){
    //std::cout << "Num cells:  " << dolfin_cells_[i].num_entities(dim) << std::endl;
    //std::cout << "Num facets: " << dolfin_cells_[i].num_entities(dim-1) << std::endl;
    
    // dolfin's own index
    const std::size_t self = dolfin_cells_[i].index();
    auto facets = dolfin_cells_[i].entities(dim-1);
    for ( std::size_t j = 0; j < dolfin_cells_[i].num_entities(dim-1); ++j ){
      //std::cout << "Index: " << facets[j] << std::endl;
      dolfin::Facet dolfin_facet(*mesh, facets[j]);
      auto neighbor_cells = dolfin_facet.entities(dim);
      for (std::size_t k = 0; k < dolfin_facet.num_entities(dim); ++k){
        //std::cout << "Neigh: " << neighbor_cells[k] << std::endl;
        if (self != neighbor_cells[k]){
          cell2cells_[i].insert(dolfin2local ? (*dolfin2local)[neighbor_cells[k]] : neighbor_cells[k]);
        }
      }
    }
  }
  // TODO: Include periodic neighbor cells
}


void label_cell_type(std::vector<int>& cell_type_, std::vector<CellNeighbours>& cell2cells_, const Uint dim){
  //cell_type_.clear(); // set all to zero
  // Cell types:
  // 0: bulk cell
  // 1: boundary cell
  // 2: next to boundary cell

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

void apply_periodic_boundaries(std::vector<CellNeighbours>& cell2cells_,
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

  // Match periodic facets within a sorted window
  for ( Uint k=0; k < dim; ++k ){
    const Uint a = (k + 1) % dim;
    auto & hi = bdry_h[k];
    std::sort(hi.begin(), hi.end(), [a](const auto& p, const auto& q){ return p.first[a] < q.first[a]; });
    for ( auto & item1 : bdry_l[k] ){
      auto x1 = item1.first;
      auto id1 = item1.second;
      const double w = 2*tol + 4*std::numeric_limits<double>::epsilon()*(std::abs(x1[a]) + 1.);
      auto it = std::lower_bound(hi.begin(), hi.end(), x1[a] - w,
                                 [a](const auto& p, const double v){ return p.first[a] < v; });
      for ( ; it != hi.end() && it->first[a] <= x1[a] + w; ++it ){
        auto x2 = it->first;
        auto id2 = it->second;

        double dx = (x1-x2).norm();

        if ( dx < tol ){
          cell2cells_[id1].insert(id2);
          cell2cells_[id2].insert(id1);
        }
      }
    }
  }
}


void build_facet_neighbours(std::vector<std::int32_t>& across,
                            std::shared_ptr<dolfin::Mesh> mesh,
                            const std::vector<dolfin::Cell>& dolfin_cells_,
                            const std::vector<std::uint32_t>* dolfin2local,
                            const std::vector<bool>& periodic,
                            const Vector3d& x_min,
                            const Vector3d& x_max,
                            const Uint dim,
                            const double tol)
{
  const std::size_t nv = dim + 1;
  const std::size_t ncells = dolfin_cells_.size();
  across.assign(ncells*nv, facet_wall);
  // Periodic facets by axis: (midpoint, slot), the low side shifted onto the high
  std::vector<std::vector<std::pair<Vector3d, std::size_t>>> lo(dim), hi(dim);
  for (std::size_t i = 0; i < ncells; ++i){
    const dolfin::Cell& cell = dolfin_cells_[i];
    const auto verts = cell.entities(0);
    const auto facets = cell.entities(dim-1);
    for (std::size_t j = 0; j < cell.num_entities(dim-1); ++j){
      dolfin::Facet f(*mesh, facets[j]);
      // Slot: the cell vertex off the facet
      const auto fv = f.entities(0);
      const auto fv_end = fv + f.num_entities(0);
      std::size_t k = 0;
      while (k < nv && std::find(fv, fv_end, verts[k]) != fv_end) ++k;
      assert(k < nv);
      const std::size_t slot = i*nv + k;
      if (!f.exterior()){
        const auto nc = f.entities(dim);
        const std::size_t other = (nc[0] == cell.index()) ? nc[1] : nc[0];
        across[slot] = std::int32_t(dolfin2local ? (*dolfin2local)[other] : other);
        continue;
      }
      const Vector3d pt(f.midpoint().coordinates());
      for (Uint a = 0; a < dim; ++a){
        if (!periodic[a]) continue;
        if (pt[a] < x_min[a] + tol){
          Vector3d q = pt;
          q[a] += x_max[a] - x_min[a];
          lo[a].push_back({q, slot});
          break;
        }
        if (pt[a] > x_max[a] - tol){
          hi[a].push_back({pt, slot});
          break;
        }
      }
    }
  }

  // Match periodic facets within a sorted window
  std::size_t unmatched = 0;
  for (Uint a = 0; a < dim; ++a){
    const Uint b = (a + 1) % dim;
    auto& h = hi[a];
    std::sort(h.begin(), h.end(), [b](const auto& p, const auto& q){ return p.first[b] < q.first[b]; });
    std::vector<bool> h_matched(h.size(), false);
    for (const auto& l : lo[a]){
      const double w = 2*tol + 4*std::numeric_limits<double>::epsilon()*(std::abs(l.first[b]) + 1.);
      auto it = std::lower_bound(h.begin(), h.end(), l.first[b] - w,
                                 [b](const auto& p, const double v){ return p.first[b] < v; });
      bool matched = false;
      for ( ; it != h.end() && it->first[b] <= l.first[b] + w; ++it){
        if ((l.first - it->first).norm() < tol){
          across[l.second] = facet_periodic(std::int32_t(it->second / nv));
          across[it->second] = facet_periodic(std::int32_t(l.second / nv));
          h_matched[it - h.begin()] = true;
          matched = true;
          break;
        }
      }
      if (!matched) ++unmatched;
    }
    unmatched += std::count(h_matched.begin(), h_matched.end(), false);
  }
  if (unmatched)
    std::cout << unmatched << " periodic facets have no partner and reflect as walls" << std::endl;
}

Vector3d periodic_lengths(const std::vector<bool>& periodic, const Vector3d& x_min,
                          const Vector3d& x_max, const Uint dim)
{
  Vector3d period = Vector3d::Zero();
  for (Uint a = 0; a < dim; ++a)
    if (periodic[a]) period[a] = x_max[a] - x_min[a];
  return period;
}

#endif