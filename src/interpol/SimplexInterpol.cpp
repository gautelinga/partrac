#include "SimplexInterpol.hpp"
#include "StampedInterpol_impl.hpp"
#include "loader_params.hpp"
#include "mesh_tables.hpp"
#include "simplex_load.hpp"
#include "phase_timing.hpp"
#include <algorithm>
#include <cstdio>
#include <iostream>

partrac::Schema DolfinH5Format::schema(const int D)
{
  return dolfin_h5_schema(D == 2 ? "triangle" : "tet");
}

template<typename I>
void DolfinH5Format::load(I& intp, const std::string& infilename)
{
  using Cell = std::decay_t<decltype(intp.cells_[0])>;
  constexpr int nv = Cell::n_verts;
  const char* mode = nv == 3 ? "triangle" : "tet";
  auto& prm = intp.dolfin_params;
  // The split field is a class of its own, which the factory picks; this loader
  // would read the same file and evaluate it as plain P2
  if (prm.template get<bool>("divfree"))
    partrac::fail(infilename, ": divfree = true is read by SplitInterpol, which the factory "
                  "picks from this key; reached here, the key was not seen as a boolean");
  const std::string folder = intp.get_folder();
  const bool include_pressure = intp.include_pressure;
  const bool include_phi = intp.include_phi;

  ts.initialize(folder + "/" + prm.template get<std::string>("timestamps"));
  partrac::phase("params");

  // The first stamp, by the path update builds for every other one
  const std::string first = folder + "/" + ts.get(ts.get_t_min()).prev.filename;
  simplex_load::Request req;
  req.infilename = infilename;
  req.field_file = first;
  req.what = "SimplexInterpol";
  req.nv = nv;
  req.include_pressure = include_pressure;
  req.include_phi = include_phi;
  req.n_dofs_max = Cell::n_dofs_max;
  req.periodic = intp.periodic;
  req.periodic_tol = intp.periodic_tol;
  // The declared elements are checked whether the tables are rebuilt or read
  // from the cache
  simplex_load::request_from_params(req, prm, folder);
  u_field = req.u_field;
  p_field = req.p_field;
  phi_field = req.phi_field;
  const std::string& mesh_file = req.mesh_file;

  // The native cache, opt-in: everything below is deterministic from the two
  // input files, so a later run reads the tables back instead of rebuilding
  const bool use_cache = prm.template get<bool>("mesh_cache");
  const std::string cache_file = mesh_file.substr(0, mesh_file.find_last_of('.'))
                               + "_partrac_" + mode + ".h5";
  std::string key;
  if (use_cache){
    key = simplex_load::cache_key({mesh_file, first});
    if (!key.empty()){
      // The facet table and the node tables are built from the periodicity too
      char tol[32];
      std::snprintf(tol, sizeof tol, "%.17g", intp.periodic_tol);
      key += "|u=" + u_field + "|p=" + (include_pressure ? p_field : std::string())
           + "|phi=" + (include_phi ? phi_field : std::string())
           + "|periodic=" + (intp.periodic[0] ? "1" : "0") + (intp.periodic[1] ? "1" : "0")
                          + (intp.periodic[2] ? "1" : "0")
           + "|periodic_tol=" + tol;
    }
  }
  auto& a = intp.stamps_.a();
  simplex_load::CacheTables cache;
  const bool cached = use_cache && simplex_load::cache_read(cache_file, key, cache)
                   && (!include_pressure || cache.ncoeffs_p > 0)
                   && (!include_phi || cache.ncoeffs_phi > 0);
  if (cached){
    std::cout << "Mesh cache: read from " << cache_file << std::endl;
    intp.dim = cache.gdim;
    intp.x_min = cache.x_min;
    intp.x_max = cache.x_max;
    intp.hmin_ = cache.hmin;
    intp.ncoeffs_u = cache.ncoeffs_u;
    intp.ncoeffs_p = include_pressure ? cache.ncoeffs_p : 0;
    intp.ncoeffs_phi = include_phi ? cache.ncoeffs_phi : 0;
    check_dofs_fit(intp.ncoeffs_u, std::max(intp.ncoeffs_p, intp.ncoeffs_phi), Cell::n_dofs_max, "SimplexInterpol");
    intp.topo_ = std::move(cache.topo);
    intp.coords_ = std::move(cache.coords);
    intp.facet_neigh_ = std::move(cache.facets);
    intp.u_dofs_.adopt(std::move(cache.u_nodes), intp.ncoeffs_u);
    a.u = std::move(cache.u_values);
    u_map = std::move(cache.u_map);
    if (include_pressure){
      intp.p_dofs_.adopt(std::move(cache.p_nodes), intp.ncoeffs_p);
      a.p = std::move(cache.p_values);
      p_map = std::move(cache.p_map);
    }
    if (include_phi){
      intp.phi_dofs_.adopt(std::move(cache.phi_nodes), intp.ncoeffs_phi);
      a.phi = std::move(cache.phi_values);
      phi_map = std::move(cache.phi_map);
    }
    intp.stamps_.hold_a(cache.stamp);
    intp.ncells_ = intp.topo_.size()/nv;
    intp.nverts_ = intp.coords_.size()/intp.dim;
    intp.set_period();
    partrac::phase("mesh cache");
    return;
  }

  simplex_load::Tables t;
  simplex_load::build_tables(req, t);

  intp.adopt_tables(t);
  intp.ncoeffs_phi = t.ncoeffs_phi;
  intp.phi_dofs_ = std::move(t.phi_dofs);
  intp.ncells_ = t.mesh.ncells;
  intp.nverts_ = t.mesh.nverts;

  simplex_load::read_field_by_node(first, u_field, t, t.el_u, t.node_order(t.el_u), a.u, u_map);
  if (include_pressure)
    simplex_load::read_field_by_node(first, p_field, t, t.el_p, t.node_order(t.el_p), a.p, p_map);
  if (include_phi)
    simplex_load::read_field_by_node(first, phi_field, t, t.el_phi, t.node_order(t.el_phi), a.phi, phi_map);
  intp.stamps_.hold_a(first);

  intp.topo_ = std::move(t.mesh.topo);
  intp.coords_ = std::move(t.mesh.coords);

  if (use_cache){
    simplex_load::CacheTables out;
    out.topo = intp.topo_;
    out.coords = intp.coords_;
    out.facets = intp.facet_neigh_;
    out.u_nodes = intp.u_dofs_.table();
    out.u_values = a.u;
    out.u_map = u_map;
    out.p_nodes = intp.p_dofs_.table();
    out.p_values = a.p;
    out.p_map = p_map;
    out.phi_nodes = intp.phi_dofs_.table();
    out.phi_values = a.phi;
    out.phi_map = phi_map;
    out.ncells = intp.ncells_;
    out.nverts = intp.nverts_;
    out.gdim = intp.dim;
    out.ncoeffs_u = intp.ncoeffs_u;
    out.ncoeffs_p = intp.ncoeffs_p;
    out.ncoeffs_phi = intp.ncoeffs_phi;
    out.ncomp_u = std::size_t(nv - 1);
    out.hmin = intp.hmin_;
    out.x_min = intp.x_min;
    out.x_max = intp.x_max;
    out.stamp = intp.stamps_.key_a();
    simplex_load::cache_write(cache_file, key, out);
    partrac::phase("mesh cache written");
  }
}

template<typename I>
void DolfinH5Format::read(I& intp, const Key& file, typename I::Stamp& s)
{
  // The node counts every stamp shares
  const auto& a = intp.stamps_.a();
  simplex_load::read_into(file, u_field, u_map, a.u, s.u);
  if (intp.include_pressure) simplex_load::read_into(file, p_field, p_map, a.p, s.p);
  if (intp.include_phi)      simplex_load::read_into(file, phi_field, phi_map, a.phi, s.phi);
}

template void DolfinH5Format::load(SimplexInterpol<Triangle>&, const std::string&);
template void DolfinH5Format::load(SimplexInterpol<Tet>&, const std::string&);
template void DolfinH5Format::read(SimplexInterpol<Triangle>&, const Key&, SimplexInterpol<Triangle>::Stamp&);
template void DolfinH5Format::read(SimplexInterpol<Tet>&, const Key&, SimplexInterpol<Tet>::Stamp&);

template class StampedInterpol<Triangle, DolfinH5Format>;
template class StampedInterpol<Tet, DolfinH5Format>;
