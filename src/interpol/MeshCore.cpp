#include "MeshCore.hpp"
#include "Triangle.hpp"
#include "Tet.hpp"

template<typename Cell>
void MeshCore<Cell>::read_mesh_params(){
  if (dolfin_params.get<bool>("periodic_x")){
    periodic[0] = true;
  }
  if (dolfin_params.get<bool>("periodic_y")){
    periodic[1] = true;
  }
  if (dolfin_params.get<bool>("periodic_z")){
    periodic[2] = true;
  }
  if (dolfin_params.get<bool>("ignore_pressure")){
    include_pressure = false;
  }
}

template<typename Cell>
void MeshCore<Cell>::set_period(){
  period_ = Vector3d::Zero();
  for (Uint d = 0; d < dim; ++d)
    if (periodic[d]) period_[d] = x_max[d] - x_min[d];
}

template<typename Cell>
bool MeshCore<Cell>::locate_tree(const Vector3d& xx, CellPos& pos)
{
  if (!tree_)
    return false;
  const int id = tree_->locate(xx);
  if (id < 0)
    return false;
  pos.id = id;
  // The exact test decided the cell; the barycentrics are the cell's own, as
  // every other path in this code computes them
  cells_[id].contains(xx, pos.bary);
  return true;
}

template<typename Cell>
void MeshCore<Cell>::enable_reflection()
{
  can_reflect = true;
}

template<typename Cell>
Vector3d MeshCore<Cell>::get_boundary_normal(const Vector3d &x, int& cell_id)
{
  if (cell_id < 0 || std::size_t(cell_id) >= cells_.size())
    return Vector3d::Zero();
  // The gradient of barycentric k points into the cell from the facet facing vertex k
  Vector3d n = Vector3d::Zero();
  for (int k = 0; k < Cell::n_verts; ++k)
    if (facet_neigh_[std::size_t(cell_id)*Cell::n_verts + k] == facet_wall)
      n -= cells_[cell_id].bary_grad(k).normalized();
  const double len = n.norm();
  return len > 0. ? Vector3d(n/len) : n;
}

template<typename Cell>
bool MeshCore<Cell>::reflect(const Vector3d& x, Vector3d& dx, CellPos& pos)
{
  return reflect_in_cells(cells_, facet_neigh_, Cell::n_verts, period_, x, dx, pos,
                          [this](const Vector3d& p){ return _modx(p); });
}

template class MeshCore<Triangle>;
template class MeshCore<Tet>;
