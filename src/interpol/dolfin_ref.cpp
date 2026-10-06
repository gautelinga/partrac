#ifdef USE_DOLFIN
#include "Error.hpp"
#include "dolfin_ref.hpp"
#include "mesh_tables.hpp"
#include "Triangle.hpp"
#include "Tet.hpp"
#include "dolfin_elements/vP1_2.h"
#include "dolfin_elements/vP2_2.h"
#include "dolfin_elements/vP3_2.h"
#include "dolfin_elements/vP1_3.h"
#include "dolfin_elements/vP2_3.h"
#include "dolfin_elements/vP3_3.h"
#include "dolfin_elements/P1_2.h"
#include "dolfin_elements/P2_2.h"
#include "dolfin_elements/P3_2.h"
#include "dolfin_elements/P1_3.h"
#include "dolfin_elements/P2_3.h"
#include "dolfin_elements/P3_3.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>

template<int D, bool Vector>
std::shared_ptr<dolfin::FunctionSpace> lagrange_space(const std::string& el,
                                                      std::shared_ptr<dolfin::Mesh> mesh,
                                                      std::shared_ptr<const dolfin::SubDomain> cd,
                                                      const char* what){
  if constexpr (D == 2 && Vector){
    if (el == "P1") return std::make_shared<vP1_2::FunctionSpace>(mesh, cd);
    if (el == "P2") return std::make_shared<vP2_2::FunctionSpace>(mesh, cd);
    if (el == "P3") return std::make_shared<vP3_2::FunctionSpace>(mesh, cd);
  }
  else if constexpr (D == 2){
    if (el == "P1") return std::make_shared<P1_2::FunctionSpace>(mesh, cd);
    if (el == "P2") return std::make_shared<P2_2::FunctionSpace>(mesh, cd);
    if (el == "P3") return std::make_shared<P3_2::FunctionSpace>(mesh, cd);
  }
  else if constexpr (Vector){
    if (el == "P1") return std::make_shared<vP1_3::FunctionSpace>(mesh, cd);
    if (el == "P2") return std::make_shared<vP2_3::FunctionSpace>(mesh, cd);
    if (el == "P3") return std::make_shared<vP3_3::FunctionSpace>(mesh, cd);
  }
  else {
    if (el == "P1") return std::make_shared<P1_3::FunctionSpace>(mesh, cd);
    if (el == "P2") return std::make_shared<P2_3::FunctionSpace>(mesh, cd);
    if (el == "P3") return std::make_shared<P3_3::FunctionSpace>(mesh, cd);
  }
  partrac::fail("unrecognized ", what, " element: ", el);
}

template std::shared_ptr<dolfin::FunctionSpace> lagrange_space<2, true>(
    const std::string&, std::shared_ptr<dolfin::Mesh>,
    std::shared_ptr<const dolfin::SubDomain>, const char*);
template std::shared_ptr<dolfin::FunctionSpace> lagrange_space<2, false>(
    const std::string&, std::shared_ptr<dolfin::Mesh>,
    std::shared_ptr<const dolfin::SubDomain>, const char*);
template std::shared_ptr<dolfin::FunctionSpace> lagrange_space<3, true>(
    const std::string&, std::shared_ptr<dolfin::Mesh>,
    std::shared_ptr<const dolfin::SubDomain>, const char*);
template std::shared_ptr<dolfin::FunctionSpace> lagrange_space<3, false>(
    const std::string&, std::shared_ptr<dolfin::Mesh>,
    std::shared_ptr<const dolfin::SubDomain>, const char*);

void CellDofs::build(const dolfin::GenericDofMap& dofmap,
                     const std::vector<dolfin::Cell>& cells, const char* what){
  stride_ = cells.empty() ? 0 : dofmap.cell_dofs(cells[0].index()).size();
  dofs_.resize(cells.size()*stride_);
  for (std::size_t id = 0; id < cells.size(); ++id){
    const auto dofs = dofmap.cell_dofs(cells[id].index());
    if (std::size_t(dofs.size()) != stride_){
      partrac::fail(what, ": cell ", id, " has ", dofs.size(), " dofs, the first cell ", stride_);
    }
    for (std::size_t i = 0; i < stride_; ++i){
      if (dofs[i] < 0 || std::uint64_t(dofs[i]) > std::numeric_limits<std::uint32_t>::max()){
        partrac::fail(what, ": dof index ", dofs[i], " does not fit 32 bits");
      }
      dofs_[id*stride_ + i] = std::uint32_t(dofs[i]);
    }
  }
}

std::vector<int> sorted_dof_table(const dolfin::GenericDofMap& dofmap,
                                  const std::size_t ncells, std::size_t& stride){
  stride = ncells ? dofmap.cell_dofs(0).size() : 0;
  std::vector<int> table(ncells * stride);
  for (std::size_t i = 0; i < ncells; ++i){
    const auto d = dofmap.cell_dofs(i);
    if (std::size_t(d.size()) != stride){
      partrac::fail("cell ", i, " has ", d.size(), " dofs, the first cell ", stride);
    }
    int* row = table.data() + i*stride;
    std::copy(d.data(), d.data() + stride, row);
    std::sort(row, row + stride);
  }
  return table;
}

double dof_sharing(const dolfin::GenericDofMap& dofmap, const std::size_t ncells){
  if (ncells < 2) return 1.;
  std::size_t stride = 0;
  const std::vector<int> table = sorted_dof_table(dofmap, ncells, stride);
  std::size_t shared = 0;
  for (std::size_t i = 1; i < ncells; ++i){
    const int* prev = table.data() + (i-1)*stride;
    const int* cur = table.data() + i*stride;
    std::size_t a = 0, b = 0;
    while (a < stride && b < stride){
      if (prev[a] == cur[b]){ ++shared; break; }
      (prev[a] < cur[b]) ? ++a : ++b;
    }
  }
  return double(shared) / (ncells - 1);
}

std::vector<std::uint32_t> order_cells_by_dofs(const dolfin::GenericDofMap& dofmap,
                                               const std::size_t ncells){
  std::size_t stride = 0;
  const std::vector<int> table = sorted_dof_table(dofmap, ncells, stride);
  std::vector<std::uint32_t> by_key(ncells);
  std::iota(by_key.begin(), by_key.end(), 0);
  std::stable_sort(by_key.begin(), by_key.end(),
                   [&](const std::uint32_t a, const std::uint32_t b){
                     const int* ra = table.data() + a*stride;
                     const int* rb = table.data() + b*stride;
                     return std::lexicographical_compare(ra, ra + stride, rb, rb + stride);
                   });
  std::vector<std::uint32_t> map(ncells);
  for (std::size_t l = 0; l < ncells; ++l) map[by_key[l]] = l;
  return map;
}

std::vector<std::uint32_t> cell_order(const dolfin::GenericDofMap& dofmap,
                                      const std::size_t ncells,
                                      const std::string& mode_in,
                                      std::vector<std::uint32_t>& dolfin2local){
  std::vector<std::uint32_t> order(ncells);
  std::iota(order.begin(), order.end(), 0);
  const std::string mode = mode_in.empty() ? "auto" : mode_in;
  if (mode != "auto" && mode != "never" && mode != "always"){
    partrac::fail("renumber_cells must be auto, never or always, not ", mode);
  }
  bool renumber = mode == "always";
  if (mode == "auto"){
    const double sharing = dof_sharing(dofmap, ncells);
    renumber = sharing < 0.5;
    std::cout << "Consecutive cells sharing a dof: " << sharing << std::endl;
  }
  if (renumber){
    std::cout << "Cell order: renumbering cells by dofs" << std::endl;
    dolfin2local = order_cells_by_dofs(dofmap, ncells);
    for (std::size_t i = 0; i < ncells; ++i) order[dolfin2local[i]] = i;
  }
  return order;
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
    partrac::fail("mesh has ", mesh->num_cells(), " cells, more than a cell id can hold");
  }
  const std::size_t nv = dim + 1;
  const std::size_t ncells = dolfin_cells_.size();
  across.assign(ncells*nv, facet_wall);
  // The exterior facets, (midpoint, slot), for the periodic match
  std::vector<std::pair<Vector3d, std::size_t>> exterior;
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
      exterior.push_back({Vector3d(f.midpoint().coordinates()), slot});
    }
  }
  mesh_tables::match_periodic_facets(across, exterior, periodic, x_min, x_max, dim, nv, tol);
}

#endif
