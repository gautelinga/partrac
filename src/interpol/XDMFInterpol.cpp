#include "Error.hpp"
#include "XDMFInterpol.hpp"
#include "StampedInterpol_impl.hpp"
#include "loader_params.hpp"
#include "mesh_tables.hpp"
#include "phase_timing.hpp"
#include "simplex_load.hpp"
#include "xdmf_helpers.hpp"
#include <iostream>

partrac::Schema XDMFFormat::schema(const int D)
{
  return xdmf_schema(D == 2 ? "xdmftriangle" : "xdmftet");
}

namespace {

// One stamp of one field by vertex, the first ncols of each stored row
void read_columns(const std::vector<std::string>& path, std::vector<double>& data,
                  const int ncols, const std::size_t nverts)
{
  read_dataset_columns(path[0], path[1], data, ncols);
  if (data.size() != nverts*std::size_t(ncols)){
    partrac::fail(path[0], ": '", path[1], "' holds ", data.size()/std::size_t(ncols),
                  " rows, the mesh has ", nverts, " vertices");
  }
}

}  // namespace

template<typename I>
void XDMFFormat::load(I& intp, const std::string&)
{
  using Cell = std::decay_t<decltype(intp.cells_[0])>;
  constexpr int nv = Cell::n_verts;
  constexpr int ne = nv*(nv-1)/2;
  auto& prm = intp.dolfin_params;
  const std::string folder = intp.get_folder();

  if constexpr (nv == 4)
    intp.periodic_tol = 1e-4;   // the XDMF tet meshes need a looser match

  // One xdmf per field; the velocity's carries the mesh
  const std::string xdmffilename_u = folder + "/" + prm.template get<std::string>("u");
  const std::string xdmffilename_p = folder + "/"
    + (prm.has("p") ? prm.template get<std::string>("p") : "");
  const std::string xdmffilename_phi = folder + "/"
    + (prm.has("phi") ? prm.template get<std::string>("phi") : "");

  std::string h5filename_u, topology_path, geometry_path;
  ts.initialize(parse_xdmf(xdmffilename_u, h5filename_u, topology_path, geometry_path));

  std::cout << "mesh: " << h5filename_u << ": " << topology_path << " " << geometry_path << std::endl;

  if (intp.include_pressure)
    ts.add("p", parse_xdmf(xdmffilename_p));
  if (intp.include_phi)
    ts.add("phi", parse_xdmf(xdmffilename_phi));
  partrac::phase("params");

  // The mesh the xdmf names, its cells along the Morton curve
  simplex_load::MeshData m;
  std::vector<std::uint32_t> cell_perm;
  simplex_load::read_mesh_arrays(h5filename_u, topology_path, geometry_path, nv, m, cell_perm);
  cell_perm.clear();
  cell_perm.shrink_to_fit();
  intp.dim = m.gdim;
  intp.x_min = m.x_min;
  intp.x_max = m.x_max;
  intp.ncells_ = m.ncells;
  intp.nverts_ = m.nverts;
  intp.set_period();

  // A periodic space gave a vertex and its images one dof, so they read one value
  intp.vclass_ = mesh_tables::match_periodic_vertices(m.coords, m.nverts, intp.dim, intp.periodic,
                                                      intp.x_min, intp.x_max, intp.periodic_tol);
  partrac::phase("periodic vertices");

  // The dofs are the vertices: no dofmap, and one node table for every field
  intp.ncoeffs_u = Uint(nv);
  intp.ncoeffs_p = Uint(nv);
  intp.ncoeffs_phi = Uint(nv);
  check_dofs_fit(intp.ncoeffs_u, intp.ncoeffs_p, Cell::n_dofs_max, "XDMFInterpol");
  intp.u_dofs_.fill(m.topo, std::vector<std::uint32_t>(), m.ncells, nv, ne, false, m.nverts,
                    intp.vclass_.data());
  intp.u_dofs_.check_stride(intp.ncoeffs_u, "XDMFInterpol");
  partrac::phase("cell dofs");

  mesh_tables::build_facet_neighbours<nv>(m.topo, m.ncells, m.coords, intp.dim, intp.periodic,
                                          intp.x_min, intp.x_max, intp.periodic_tol,
                                          intp.facet_neigh_);
  partrac::phase("facet table");

  intp.hmin_ = simplex_load::shortest_edge(m, nv);
  intp.topo_ = std::move(m.topo);
  intp.coords_ = std::move(m.coords);
}

template<typename I>
void XDMFFormat::read(I& intp, const Key& it, typename I::Stamp& s)
{
  read_columns(ts.get_path("u", it), s.u, int(intp.dim), intp.nverts_);
  if (intp.include_pressure)
    read_columns(ts.get_path("p", it), s.p, 1, intp.nverts_);
  if (intp.include_phi)
    read_columns(ts.get_path("phi", it), s.phi, 1, intp.nverts_);
}

template void XDMFFormat::load(XDMFInterpol<Triangle>&, const std::string&);
template void XDMFFormat::load(XDMFInterpol<Tet>&, const std::string&);
template void XDMFFormat::read(XDMFInterpol<Triangle>&, const Key&, XDMFInterpol<Triangle>::Stamp&);
template void XDMFFormat::read(XDMFInterpol<Tet>&, const Key&, XDMFInterpol<Tet>::Stamp&);

template class StampedInterpol<Triangle, XDMFFormat>;
template class StampedInterpol<Tet, XDMFFormat>;
