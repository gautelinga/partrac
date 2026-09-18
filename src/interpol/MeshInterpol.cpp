#ifdef USE_DOLFIN
#include "MeshInterpol.hpp"
#include "Triangle.hpp"
#include "Tet.hpp"
#include "PeriodicBC.hpp"
#include "dolfin_helpers.hpp"
#include "geometry.hpp"
#include <cassert>

template<typename Cell>
void MeshInterpol<Cell>::read_mesh_params(){
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
void MeshInterpol<Cell>::init_mesh_geometry(){
  dim = mesh->geometry().dim();
  mesh->init();
  mesh->bounding_box_tree();

  std::vector<double> xx = mesh->coordinates();

  for (Uint i=0; i<dim; ++i){
    x_min[i] = xx[i];
    x_max[i] = xx[i];
  }

  for (Uint i=0; i<xx.size(); ++i){
    Uint i_loc = i % dim;
    x_min[i_loc] = std::min(x_min[i_loc], xx[i]);
    x_max[i_loc] = std::max(x_max[i_loc], xx[i]);
  }
}

template<typename Cell>
void MeshInterpol<Cell>::build_cells(const dolfin::GenericDofMap& dofmap){
  const std::size_t ncells = mesh->num_cells();
  cells_.resize(ncells);
  dolfin_cells_.resize(ncells);

  const std::vector<std::uint32_t> order =
    cell_order(dofmap, ncells, dolfin_params.get<std::string>("renumber_cells"), dolfin2local_);
  for (std::size_t l = 0; l < ncells; ++l)
  {
    dolfin::Cell dolfin_cell(*mesh, order[l]);
    cells_[l] = Cell(dolfin_cell);
    dolfin_cells_[l] = dolfin_cell;
  }
  build_facet_table();
}

template<typename Cell>
bool MeshInterpol<Cell>::locate_tree(const Vector3d& xx, CellPos& pos)
{
  return tree_to_cell(cells_, *mesh, dim, xx, pos,
                      dolfin2local_.empty() ? nullptr : &dolfin2local_);
}

template<typename Cell>
void MeshInterpol<Cell>::build_facet_table()
{
  build_facet_neighbours(facet_neigh_, mesh, dolfin_cells_,
                         dolfin2local_.empty() ? nullptr : &dolfin2local_,
                         periodic, x_min, x_max, dim, periodic_tol);
  period_ = periodic_lengths(periodic, x_min, x_max, dim);
}

template<typename Cell>
void MeshInterpol<Cell>::enable_reflection()
{
  can_reflect = true;
}

template<typename Cell>
Vector3d MeshInterpol<Cell>::get_boundary_normal(const Vector3d &x, int& cell_id)
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
bool MeshInterpol<Cell>::reflect(const Vector3d& x, Vector3d& dx, CellPos& pos)
{
  return reflect_in_cells(cells_, facet_neigh_, Cell::n_verts, period_, x, dx, pos,
                          [this](const Vector3d& p){ return _modx(p); });
}

template class MeshInterpol<Triangle>;
template class MeshInterpol<Tet>;

#endif
