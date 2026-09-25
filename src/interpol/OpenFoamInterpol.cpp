#include "OpenFoamInterpol.hpp"
#include "StampedInterpol_impl.hpp"
#include "Error.hpp"
#include "loader_params.hpp"
#include "mesh_tables.hpp"
#include "openfoam_split.hpp"
#include "phase_timing.hpp"
#include "simplex_load.hpp"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <sstream>

namespace {
const char* const phase_note = "; phase_field names a volScalarField such as alpha.water "
                               "(OpenFOAM's phi is the volume flux through the faces)";

// A scalar field's conditions OpenFOAM cannot build, into the log
void report_unknown(const openfoam_load::FieldData& f, const std::vector<openfoam_load::Patch>& patches){
  for (std::size_t i = 0; i < f.patches.size(); ++i)
    if (f.patches[i].unknown)
      std::cout << f.name << ": patch " << patches[i].name << " is " << f.patches[i].condition
                << ", which OpenFOAM cannot build here (its library is not loaded); its faces are not data"
                << std::endl;
}

// A volume's error relative to its target; n/a for a target of zero
std::string relative(const double v, const double target){
  if (target == 0.) return "n/a";
  std::ostringstream s;
  s << (v - target)/target;
  return s.str();
}
}

partrac::Schema OpenFoamFormat::schema(const int)
{
  return openfoam_schema();
}

void OpenFoamFormat::report_boundary(const std::vector<openfoam_load::Patch>& patches,
                                     const std::vector<std::int32_t>& facet_patch) const
{
  using openfoam_nodes::PatchClass;
  std::vector<std::size_t> per_patch(patches.size(), 0);
  for (const std::int32_t p : facet_patch) if (p >= 0) ++per_patch[std::size_t(p)];
  std::cout << "Boundary facets:";
  for (const PatchClass k : {PatchClass::wall, PatchClass::moving_wall, PatchClass::cyclic, PatchClass::other}){
    std::size_t n = 0;
    std::string names;
    for (std::size_t p = 0; p < patches.size(); ++p){
      if (patch_class[p] != k || !per_patch[p]) continue;
      n += per_patch[p];
      names += (names.empty() ? "" : " ") + patches[p].name;
    }
    std::cout << (k == PatchClass::wall ? " " : ", ") << n << " " << openfoam_nodes::class_name(k)
              << (names.empty() ? "" : " (" + names + ")");
  }
  std::cout << std::endl;
}

template<typename I>
void OpenFoamFormat::load(I& intp, const std::string& infilename)
{
  using Cell = std::decay_t<decltype(intp.cells_[0])>;
  constexpr int nv = Cell::n_verts;
  auto& prm = intp.dolfin_params;
  const std::string folder = intp.get_folder();
  u_field = prm.template get<std::string>("velocity_field");
  p_field = prm.template get<std::string>("pressure_field");
  phi_field = prm.has("phase_field") ? prm.template get<std::string>("phase_field") : std::string();
  intp.include_phi = !phi_field.empty();
  const int tets_per_hex = prm.template get<std::string>("split") == "6" ? 6 : 12;

  openfoam_load::CaseData c = openfoam_load::read_case(folder, u_field);
  partrac::phase("case read");
  if ((c.empty_axis >= 0) != (nv == 3))
    partrac::fail(folder, ": a ", c.empty_axis >= 0 ? "2D" : "3D", " case read as ",
                  nv == 3 ? "triangles" : "tets");

  // The split
  openfoam_split::SplitData s = openfoam_split::split(c, tets_per_hex);
  {
    std::cout << "Split " << (nv == 4 ? "tets" : "triangles") << ": " << s.nsimplices() << " for "
              << c.ncells << " cells, " << s.nnodes() << " nodes";
    // The counts that are not zero
    const char* sep = "; ";
    const auto count = [&](const std::size_t n, const std::string& what){
      if (!n) return;
      std::cout << sep << n << " " << what;
      sep = ", ";
    };
    count(s.fan_cells, "cells fanned");
    count(s.five_tet_hexes, "hexes in 5");
    count(s.fallback, "split the other way");
    std::ostringstream invalid;
    invalid << "cells with " << s.invalid_simplices << " simplices under " << openfoam_split::valid_fraction
            << " of their cell";
    count(s.invalid_cells, invalid.str());
    count(s.base_failures, "faces without a valid base");
    std::cout << std::endl;
  }
  partrac::phase("split");

  // The node values' operators, from the first stamp's conditions
  first_u = openfoam_load::read_field(c, c.times.front(), u_field);
  if (first_u.ncomp != 3)
    partrac::fail(folder, "/", c.times.front(), "/", u_field, " has ", first_u.ncomp, " components, not a vector");
  if (intp.include_pressure){
    const std::string p0 = folder + "/" + c.times.front() + "/" + p_field;
    if (!std::filesystem::exists(p0) && !std::filesystem::exists(p0 + ".gz"))
      partrac::fail(p0, ": no such file; set pressure_field, or ignore_pressure=true");
    first_p = openfoam_load::read_field(c, c.times.front(), p_field);
    if (first_p.ncomp != 1)
      partrac::fail(folder, "/", c.times.front(), "/", p_field, " has ", first_p.ncomp, " components, not a scalar");
    report_unknown(first_p, c.patches);
  }
  if (intp.include_phi){
    for (const std::string& t : c.times){
      const std::string f = folder + "/" + t + "/" + phi_field;
      if (!std::filesystem::exists(f) && !std::filesystem::exists(f + ".gz"))
        partrac::fail(f, ": no such file; the phase field must be in every time directory");
    }
    first_phi = openfoam_load::read_field(c, c.times.front(), phi_field);
    if (first_phi.ncomp != 1)
      partrac::fail(folder, "/", c.times.front(), "/", phi_field, " has ", first_phi.ncomp,
                    " components, not a scalar", phase_note);
    report_unknown(first_phi, c.patches);
  }
  patch_class = openfoam_nodes::classify_patches(c, first_u);
  openfoam_nodes::check_walled(c, patch_class, tets_per_hex, std::cerr);
  partrac::phase("first stamp read");
  least_squares = prm.template get<std::string>("nodes") == "least_squares";
  {
    const openfoam_nodes::Geometry g(c);
    if (least_squares){
      w_u = openfoam_nodes::LeastSquares(g, s, first_u);
      if (intp.include_pressure) w_p = openfoam_nodes::LeastSquares(g, s, first_p);
      for (const auto* w : {&w_u, &w_p}){
        if (!w->nnodes()) continue;
        std::cout << "Nodes by least squares, " << (w == &w_u ? u_field : p_field) << ": "
                  << w->nrows() << " rows, " << w->nonzeros() << " nonzeros, " << double(w->bytes())/1e6 << " MB; "
                  << w->fixed << " fixed, "
                  << w->second_ring << " on the second ring (" << w->deficient << " still deficient), "
                  << w->zero_won << " where zero wins, " << w->mirrored << " mirrored" << std::endl;
      }
    }
    if (!least_squares || intp.include_phi){
      idw = openfoam_nodes::InverseDistance(g, s);
      std::cout << "Nodes by inverse distance" << (least_squares ? ", " + phi_field : std::string()) << ", "
                << double(idw.bytes())/1e6 << " MB" << std::endl;
    }
    if (intp.include_phi){
      g_phi = openfoam_nodes::Gradient(g, s);
      std::cout << "Gradient of " << phi_field << " by least squares: " << g_phi.nrows() << " rows, "
                << g_phi.nonzeros() << " nonzeros, " << double(g_phi.bytes())/1e6 << " MB; "
                << g_phi.second_ring << " on the second ring (" << g_phi.deficient << " still deficient), "
                << g_phi.mirrored << " mirrored" << std::endl;
    }
  }
  comps.clear();
  if (nv == 3) comps = s.inplane;
  else comps = {0, 1, 2};
  partrac::phase("node weights");
  report_boundary(c.patches, s.facet_patch);
  std::vector<std::int32_t>().swap(s.facet_patch);

  // The time directories, and what a field's read needs of the case; its mesh arrays are done with
  std::vector<std::pair<double, std::string>> times;
  for (std::size_t i = 0; i < c.times.size(); ++i) times.push_back({c.time_values[i], c.times[i]});
  layout.dir = c.dir;
  layout.ncells = c.ncells;
  layout.n_internal = c.n_internal;
  layout.patches = c.patches;
  const int empty_axis = c.empty_axis;
  c = openfoam_load::CaseData();

  // The periodic axes and the box from the cyclic pairs
  std::vector<bool> periodic = {false, false, false};
  std::vector<double> period = {0., 0., 0.};
  for (const auto& p : layout.patches){
    if (p.type != "cyclic" || !p.owner) continue;
    int axis = -1;
    const double len = std::sqrt(p.separation[0]*p.separation[0] + p.separation[1]*p.separation[1]
                                 + p.separation[2]*p.separation[2]);
    for (int a = 0; a < 3; ++a)
      if (std::abs(p.separation[std::size_t(a)]) > (1. - 1e-12)*len) axis = a;
    if (axis < 0 || axis == empty_axis)
      partrac::fail(folder, ": cyclic ", p.name, " is separated by (", p.separation[0], ", ",
                    p.separation[1], ", ", p.separation[2], "), not along an in-plane axis");
    const int a = nv == 3 ? int(std::find(s.inplane.begin(), s.inplane.end(), axis) - s.inplane.begin()) : axis;
    periodic[std::size_t(a)] = true;
    period[std::size_t(a)] = len;
  }

  // The split's arrays, cells and then vertices along the Morton curve
  simplex_load::Tables t;
  simplex_load::MeshData& m = t.mesh;
  m.gdim = Uint(nv - 1);
  m.ncells = s.nsimplices();
  m.nverts = s.nnodes();
  m.topo = std::move(s.cells);
  m.coords = std::move(s.node_x);
  std::vector<std::uint32_t> perm;
  simplex_load::morton_cells(m, nv, perm);
  std::vector<std::uint32_t> node_map;
  if (simplex_load::node_span(m, nv) > simplex_load::node_span_max){
    node_map = simplex_load::morton_node_order(m, std::vector<std::uint32_t>(), 0, nv);
    simplex_load::renumber_vertices(m, nv, node_map);
  }
  slot.resize(m.nverts);
  for (std::size_t n = 0; n < m.nverts; ++n) slot[n] = node_map.empty() ? std::uint32_t(n) : node_map[n];
  if (intp.include_phi){
    // Each OpenFOAM cell's rows and centre, as the tables number them
    std::vector<std::int32_t> of(perm.size());
    for (std::size_t i = 0; i < perm.size(); ++i) of[i] = s.cell_of[perm[i]];
    std::vector<std::int32_t> centre(std::size_t(layout.ncells), -1);
    for (std::size_t n = 0; n < s.nnodes(); ++n)
      if (s.node_kind[n]) centre[std::size_t(s.node_cell[n])] = std::int32_t(slot[n]);
    phase_cells = openfoam_phase::Cells(std::size_t(layout.ncells), of, centre);
  }
  for (int a = 0; a < nv - 1; ++a)
    if (periodic[std::size_t(a)] && std::abs(m.x_max[a] - m.x_min[a] - period[std::size_t(a)]) > 1e-9*period[std::size_t(a)])
      partrac::fail(folder, ": the cyclic pair along ", "xyz"[a], " is ", period[std::size_t(a)],
                    " apart, but the mesh spans ", m.x_max[a] - m.x_min[a], "; the periodic box is the mesh's");

  // A vertex and its cyclic images read one value, the lowest image's
  t.np.master.resize(m.nverts);
  for (std::size_t n = 0; n < m.nverts; ++n) t.np.master[slot[n]] = slot[s.node_master[n]];
  t.nv = nv;
  t.el_u = {1, std::size_t(nv - 1)};
  t.el_p = {1, 1};
  t.ncoeffs_u = Uint(nv);
  t.ncoeffs_p = intp.include_pressure ? Uint(nv) : 0;
  intp.ncoeffs_phi = intp.include_phi ? Uint(nv) : 0;
  simplex_load::Request r;
  r.infilename = infilename;
  r.what = "OpenFoamInterpol";
  r.nv = nv;
  r.include_pressure = false;   // p reads the velocity's table
  r.n_dofs_max = Cell::n_dofs_max;
  check_dofs_fit(t.ncoeffs_u, t.ncoeffs_p, Cell::n_dofs_max, r.what);
  simplex_load::tables_from_mesh(r, t);

  // The cyclic facets paired by the case's addressing
  std::vector<std::uint32_t> row_of(perm.size());
  for (std::size_t i = 0; i < perm.size(); ++i) row_of[perm[i]] = std::uint32_t(i);
  for (std::size_t q = 0; q < s.facet_partner.size(); ++q){
    const std::int64_t o = s.facet_partner[q];
    if (o < 0) continue;
    const std::size_t here = std::size_t(row_of[q/nv])*nv + q % nv;
    t.facets[here] = facet_periodic(std::int32_t(row_of[std::size_t(o)/nv]));
  }
  partrac::phase("cyclic facets");

  intp.periodic = periodic;
  intp.vclass_ = t.np.master;
  intp.adopt_tables(t);
  intp.ncells_ = m.ncells;
  intp.nverts_ = m.nverts;
  intp.topo_ = std::move(m.topo);
  intp.coords_ = std::move(m.coords);

  ts.initialize(times);
  auto& a = intp.stamps_.a();
  intp.phase_gradient_ = intp.include_phi;
  read(intp, times.front().second, a);
  intp.stamps_.hold_a(times.front().second);
  partrac::phase("first stamp");
}

template<typename I>
void OpenFoamFormat::read(I& intp, const Key& time, typename I::Stamp& s)
{
  using Cell = std::decay_t<decltype(intp.cells_[0])>;
  constexpr int nv = Cell::n_verts;
  // Node values, k a node, by vertex
  const auto place_values = [&](const std::vector<double>& v, const std::size_t k, std::vector<double>& out){
    out.assign(v.size(), 0.);
    for (std::size_t n = 0; n < slot.size(); ++n)
      for (std::size_t j = 0; j < k; ++j) out[std::size_t(slot[n])*k + j] = v[n*k + j];
  };
  const auto place = [&](const openfoam_load::FieldData& f, const openfoam_nodes::LeastSquares& w,
                         const std::vector<int>& cs, std::vector<double>& out){
    std::vector<double> v;
    if (least_squares) w.apply(f, cs, v);
    else idw.apply(f, cs, v);
    place_values(v, cs.size(), out);
  };
  // The first stamp's fields once, as read at load
  const auto field = [&](const std::string& name, openfoam_load::FieldData& first){
    openfoam_load::FieldData f;
    if (!first.name.empty() && first.time == time) f = std::move(first);
    else f = openfoam_load::read_field(layout, time, name);
    first = openfoam_load::FieldData();
    return f;
  };
  const openfoam_load::FieldData u = field(u_field, first_u);
  if (u.ncomp != 3)
    partrac::fail(layout.dir, "/", time, "/", u_field, " has ", u.ncomp, " components, not a vector");
  place(u, w_u, comps, s.u);
  if (intp.include_pressure){
    const openfoam_load::FieldData p = field(p_field, first_p);
    if (p.ncomp != 1)
      partrac::fail(layout.dir, "/", time, "/", p_field, " has ", p.ncomp, " components, not a scalar");
    place(p, w_p, {0}, s.p);
  }
  if (intp.include_phi){
    const openfoam_load::FieldData phi = field(phi_field, first_phi);
    if (phi.ncomp != 1)
      partrac::fail(layout.dir, "/", time, "/", phi_field, " has ", phi.ncomp, " components, not a scalar", phase_note);
    std::vector<double> v;
    idw.apply(phi, {0}, v);
    place_values(v, 1, s.phi);
    const openfoam_phase::Report r = openfoam_phase::conserve(phase_cells, intp.topo_, intp.coords_, nv,
                                                              phi.internal, s.phi);
    std::cout << "Phase " << phi_field << " at " << time << ": volume above 1/2 " << r.after << " for "
              << r.target << ", relative error " << relative(r.before, r.target) << " at the nodes, "
              << relative(r.after, r.target) << " corrected; worst cell " << r.worst_before << " of its volume, then "
              << r.worst_after << " (cell " << r.worst_cell << "); cells:";
    // The counts that are not zero
    const char* sep = " ";
    const auto count = [&](const std::size_t n, const char* what){
      if (!n) return;
      std::cout << sep << n << " " << what;
      sep = ", ";
    };
    count(r.pure, "pure");
    count(r.corrected, "corrected");
    count(r.clipped, "clipped to [0, 1]");
    count(r.unreachable, "unreachable");
    count(r.unfanned, "off without a centre");
    std::cout << std::endl;
    // Its gradient after its values
    std::vector<double> g;
    g_phi.apply(phi, v);
    place_values(v, std::size_t(nv - 1), g);
    openfoam_phase::centre_means(phase_cells, intp.topo_, nv, nv - 1, g);
    s.phi.insert(s.phi.end(), g.begin(), g.end());
  }
}

template void OpenFoamFormat::load(OpenFoamInterpol<Triangle>&, const std::string&);
template void OpenFoamFormat::load(OpenFoamInterpol<Tet>&, const std::string&);
template void OpenFoamFormat::read(OpenFoamInterpol<Triangle>&, const Key&, OpenFoamInterpol<Triangle>::Stamp&);
template void OpenFoamFormat::read(OpenFoamInterpol<Tet>&, const Key&, OpenFoamInterpol<Tet>::Stamp&);

template class StampedInterpol<Triangle, OpenFoamFormat>;
template class StampedInterpol<Tet, OpenFoamFormat>;
