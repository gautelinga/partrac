#include "dolfin_helpers.hpp"
#ifdef USE_DOLFIN
#include <algorithm>
#include <cmath>

// Cell types: 0 bulk, 1 on a wall, 2 across a facet from one on a wall
void label_cell_type(std::vector<int>& cell_type_, const std::vector<std::int32_t>& across, const Uint nv){
  const std::size_t ncells = across.size() / nv;
  cell_type_.assign(ncells, 0);
  for (std::size_t i = 0; i < ncells; ++i)
    for (Uint k = 0; k < nv; ++k)
      if (across[i*nv + k] == facet_wall) cell_type_[i] = 1;
  for (std::size_t i = 0; i < ncells; ++i){
    if (cell_type_[i] != 1) continue;
    for (Uint k = 0; k < nv; ++k){
      const std::int32_t a = across[i*nv + k];
      if (a == facet_wall) continue;
      const std::int32_t j = a >= 0 ? a : facet_periodic(a);
      if (cell_type_[j] == 0) cell_type_[j] = 2;
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
  // Cell ids must fit int
  if (mesh->num_cells() > std::size_t(std::numeric_limits<int>::max())){
    std::cout << "Mesh has " << mesh->num_cells() << " cells, more than a cell id can hold" << std::endl;
    exit(1);
  }
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