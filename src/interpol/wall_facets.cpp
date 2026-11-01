#include "wall_facets.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>
#include <parallel/algorithm>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Error.hpp"
#include "Params.hpp"
#include "cell_walk.hpp"
#include "h5direct.hpp"

namespace wall_facets {

namespace {

namespace pt = boost::property_tree;
using partrac::trim;

// A point matches a vertex within this much of the shortest boundary edge there
constexpr double rel_tol = 1e-3;

// An integer, the whole of s
std::int64_t integer(const std::string& s, const std::string& list){
  std::size_t end = 0;
  long long v = 0;
  try { v = std::stoll(s, &end); }
  catch (const std::exception&){ end = 0; }
  if (s.empty() || end != s.size())
    partrac::fail("wall_markers = ", list, ": '", s, "' is not a marker; give values and ranges as 3,5,7:106");
  return std::int64_t(v);
}

// The HDF5 file and dataset of a DataItem, as "file.h5:/path" beside the XDMF
std::pair<std::string, std::string> data_item(const pt::ptree& node, const std::string& what,
                                              const std::string& xdmf){
  const auto item = node.get_child_optional("DataItem");
  if (!item) partrac::fail(xdmf, ": its ", what, " has no DataItem");
  const std::string format = item->get<std::string>("<xmlattr>.Format", "HDF");
  if (format != "HDF")
    partrac::fail(xdmf, ": its ", what, " is of Format ", format, "; only HDF data is read");
  const std::string text = trim(item->data());
  const std::size_t colon = text.find(':');
  if (colon == std::string::npos)
    partrac::fail(xdmf, ": its ", what, " names no dataset, as file.h5:/data");
  std::filesystem::path file = text.substr(0, colon);
  if (file.is_relative()) file = std::filesystem::path(xdmf).parent_path() / file;
  return {file.string(), text.substr(colon + 1)};
}

// The first n (2 or 3) of k ascending
void sort_key(std::array<std::uint32_t, 3>& k, const int n){
  if (k[0] > k[1]) std::swap(k[0], k[1]);
  if (n == 3){
    if (k[1] > k[2]) std::swap(k[1], k[2]);
    if (k[0] > k[1]) std::swap(k[0], k[1]);
  }
}

}  // namespace

Markers::Markers(const std::string& list)
{
  if (trim(list).empty())
    partrac::fail("wall_markers is empty; give the markers of the no-slip walls, as 3,5,7:106");
  std::size_t b = 0;
  while (b <= list.size()){
    const std::size_t e = std::min(list.find(',', b), list.size());
    const std::string tok = trim(list.substr(b, e - b));
    const std::size_t colon = tok.find(':');
    std::int64_t lo, hi;
    if (colon == std::string::npos) lo = hi = integer(tok, list);
    else {
      lo = integer(trim(tok.substr(0, colon)), list);
      hi = integer(trim(tok.substr(colon + 1)), list);
      if (hi < lo) partrac::fail("wall_markers = ", list, ": the range ", tok, " is empty");
    }
    ranges_.push_back({lo, hi});
    b = e + 1;
  }
}

bool Markers::selects(const std::int64_t m) const
{
  for (const auto& r : ranges_)
    if (m >= r.first && m <= r.second) return true;
  return false;
}

File read(const std::string& xdmf, const Uint gdim)
{
  if (!std::filesystem::exists(xdmf)) partrac::fail("wall_facets: ", xdmf, " does not exist");
  pt::ptree tree;
  try { pt::read_xml(xdmf, tree); }
  catch (const pt::xml_parser_error& e){ partrac::fail(xdmf, ": not an XDMF file (", e.what(), ")"); }
  const auto top = tree.get_child_optional("Xdmf.Domain.Grid");
  if (!top) partrac::fail(xdmf, ": no Grid in its Domain");
  // A meshio file's grid, else the first of a collection
  const pt::ptree* g = &*top;
  if (!g->get_child_optional("Topology"))
    if (const auto inner = g->get_child_optional("Grid")) g = &*inner;
  const auto topology = g->get_child_optional("Topology");
  const auto geometry = g->get_child_optional("Geometry");
  if (!topology || !geometry) partrac::fail(xdmf, ": its grid has no Topology or no Geometry");
  // The markers: the attribute named so, else the grid's one cell attribute
  const pt::ptree* marker = nullptr;
  std::size_t cell_attributes = 0;
  for (const auto& c : *g){
    if (c.first != "Attribute") continue;
    if (c.second.get<std::string>("<xmlattr>.Center", "Cell") != "Cell") continue;
    ++cell_attributes;
    if (c.second.get<std::string>("<xmlattr>.Name", "") == "markers") marker = &c.second;
    else if (!marker && cell_attributes == 1) marker = &c.second;
  }
  if (cell_attributes > 1 && marker && marker->get<std::string>("<xmlattr>.Name", "") != "markers")
    partrac::fail(xdmf, ": ", cell_attributes, " cell attributes, none named markers");

  File f;
  const auto [geo_file, geo_set] = data_item(*geometry, "Geometry", xdmf);
  const auto [top_file, top_set] = data_item(*topology, "Topology", xdmf);
  {
    const partrac::H5Id h = partrac::h5_open_read(geo_file);
    const partrac::H5DatasetInfo info = partrac::h5_dataset_info(h, geo_set);
    if (info.rank() != 2 || info.cols() < gdim)
      partrac::fail(xdmf, ": its points (", geo_file, ":", geo_set, ") have ", info.cols(),
                    " coordinates, the mesh ", gdim);
    partrac::h5_read(h, geo_set, f.points, gdim);
  }
  {
    const partrac::H5Id h = partrac::h5_open_read(top_file);
    const partrac::H5DatasetInfo info = partrac::h5_dataset_info(h, top_set);
    if (info.rank() != 2 || info.cols() != gdim)
      partrac::fail(xdmf, ": its facets (", top_file, ":", top_set, ") have ", info.cols(),
                    " points; a ", gdim, "D mesh's have ", gdim);
    partrac::h5_read(h, top_set, f.facets);
  }
  if (marker){
    const auto [m_file, m_set] = data_item(*marker, "markers", xdmf);
    const partrac::H5Id h = partrac::h5_open_read(m_file);
    const partrac::H5DatasetInfo info = partrac::h5_dataset_info(h, m_set);
    if (info.cols() != 1 || info.rows() != f.facets.size()/gdim)
      partrac::fail(xdmf, ": ", info.count(), " markers (", m_file, ":", m_set, ") for ",
                    f.facets.size()/gdim, " facets");
    partrac::h5_read(h, m_set, f.markers);
  }
  return f;
}

std::vector<std::uint8_t> match(const File& f, const Markers& sel, const std::vector<std::uint32_t>& topo,
                                const std::vector<double>& coords, const std::size_t ncells,
                                const Uint gdim, const std::vector<std::int32_t>& facet_neigh,
                                const std::string& name, const bool verbose)
{
  const int nv = int(gdim) + 1, D = int(gdim);
  const std::size_t nverts = coords.size()/gdim, npoints = f.points.size()/gdim;
  const std::size_t nfacets = f.facets.size()/gdim;
  if (f.markers.empty())
    partrac::fail(name, " holds no markers, which wall_markers selects the walls by");
  if (f.markers.size() != nfacets)
    partrac::fail(name, ": ", f.markers.size(), " markers for ", nfacets, " facets");
  const auto facet_name = [&](const std::size_t j){
    return "facet " + std::to_string(j) + " (marker " + std::to_string(f.markers[j]) + ") of " + name;
  };

  // The selected facets, and the points they need
  std::vector<std::size_t> chosen;
  std::vector<std::int64_t> vertex_of(npoints, -1);
  constexpr std::int64_t needed = -2, unmatched = -3;
  for (std::size_t j = 0; j < nfacets; ++j){
    if (!sel.selects(f.markers[j])) continue;
    chosen.push_back(j);
    for (int a = 0; a < D; ++a){
      const std::int64_t p = f.facets[j*gdim + std::size_t(a)];
      if (p < 0 || std::size_t(p) >= npoints)
        partrac::fail(facet_name(j), ": point ", p, " of ", npoints);
      vertex_of[std::size_t(p)] = needed;
    }
  }
  if (chosen.empty())
    partrac::fail(name, ": wall_markers selects none of its ", nfacets, " facets");

  // The boundary vertices and the shortest boundary edge at each
  std::vector<double> h(nverts, std::numeric_limits<double>::infinity());
  const auto dist = [&](const double* x, const double* y){
    double s = 0.;
    for (int d = 0; d < D; ++d) s += (x[d] - y[d])*(x[d] - y[d]);
    return std::sqrt(s);
  };
  for (std::size_t s = 0; s < ncells*std::size_t(nv); ++s){
    if (facet_neigh[s] >= 0) continue;
    const std::uint32_t* row = topo.data() + (s/std::size_t(nv))*std::size_t(nv);
    std::uint32_t w[3];
    for (int j = 0, m = 0; j < nv; ++j)
      if (j != int(s % std::size_t(nv))) w[m++] = row[j];
    for (int a = 0; a < D; ++a)
      for (int b = a + 1; b < D; ++b){
        const double l = dist(&coords[w[a]*gdim], &coords[w[b]*gdim]);
        h[w[a]] = std::min(h[w[a]], l);
        h[w[b]] = std::min(h[w[b]], l);
      }
  }

  // The boundary vertices by the bits of their coordinates, -0 as 0, for the exact matches
  using Bits = std::array<std::uint64_t, 3>;
  const auto bits_of = [&](const double* x){
    Bits b{0, 0, 0};
    for (int d = 0; d < D; ++d){
      const double c = x[d] == 0. ? 0. : x[d];
      std::memcpy(&b[std::size_t(d)], &c, sizeof c);
    }
    return b;
  };
  struct Exact { Bits bits; std::uint32_t v; };
  std::vector<Exact> exact;
  for (std::size_t v = 0; v < nverts; ++v)
    if (std::isfinite(h[v])) exact.push_back({bits_of(&coords[v*gdim]), std::uint32_t(v)});
  if (exact.empty()) partrac::fail(name, ": the mesh has no boundary edge to match its points by");
  const auto by_bits = [](const Exact& a, const Exact& b){ return a.bits < b.bits || (a.bits == b.bits && a.v < b.v); };
  __gnu_parallel::sort(exact.begin(), exact.end(), by_bits);

  // Each needed point's vertex where its coordinates are one's bit for bit
  std::size_t inexact = 0;
#pragma omp parallel for schedule(dynamic, 1024) reduction(+:inexact)
  for (std::size_t p = 0; p < npoints; ++p){
    if (vertex_of[p] != needed) continue;
    const Exact e{bits_of(&f.points[p*gdim]), 0};
    const auto at = std::lower_bound(exact.begin(), exact.end(), e, by_bits);
    if (at != exact.end() && at->bits == e.bits) vertex_of[p] = at->v;
    else ++inexact;
  }

  // The others: the nearest within its tolerance, the lowest of equals
  if (inexact){
    // The candidates by level, a vertex's tolerance below half its level's
    // box, then by box: a point's within reach in two boxes an axis
    using Key = std::array<std::int64_t, 3>;
    struct Boxed { int level; Key key; std::uint32_t v; };
    static constexpr double key_max = 4.6e18;
    const auto box_of = [](const double c){ return std::int64_t(std::clamp(std::floor(c), -key_max, key_max)); };
    std::vector<Boxed> grid;
    for (std::size_t v = 0; v < nverts; ++v){
      if (!std::isfinite(h[v]) || !(rel_tol*h[v] > 0.)) continue;
      Boxed b{std::ilogb(rel_tol*h[v]) + 2, {0, 0, 0}, std::uint32_t(v)};
      for (int d = 0; d < D; ++d) b.key[std::size_t(d)] = box_of(std::ldexp(coords[v*gdim + std::size_t(d)], -b.level));
      grid.push_back(b);
    }
    const auto by_box = [](const Boxed& a, const Boxed& b){
      return a.level < b.level || (a.level == b.level && (a.key < b.key || (a.key == b.key && a.v < b.v)));
    };
    __gnu_parallel::sort(grid.begin(), grid.end(), by_box);
    std::vector<std::size_t> level_at;   // each level's first entry, then the end
    for (std::size_t i = 0; i < grid.size(); ++i)
      if (i == 0 || grid[i].level != grid[i - 1].level) level_at.push_back(i);
    level_at.push_back(grid.size());

#pragma omp parallel for schedule(dynamic, 1024)
    for (std::size_t p = 0; p < npoints; ++p){
      if (vertex_of[p] != needed) continue;
      const double* x = &f.points[p*gdim];
      std::int64_t best = unmatched;
      double best_d = std::numeric_limits<double>::infinity();
      for (std::size_t l = 0; l + 1 < level_at.size(); ++l){
        const auto first = grid.begin() + std::ptrdiff_t(level_at[l]), last = grid.begin() + std::ptrdiff_t(level_at[l + 1]);
        const int level = first->level;
        Key lo{0, 0, 0}, hi{0, 0, 0};
        for (int d = 0; d < D; ++d){
          const double c = std::ldexp(x[d], -level);
          lo[std::size_t(d)] = box_of(c - 0.5);
          hi[std::size_t(d)] = box_of(c + 0.5);
        }
        Key q;
        for (q[0] = lo[0]; q[0] <= hi[0]; ++q[0])
          for (q[1] = lo[1]; q[1] <= hi[1]; ++q[1])
            for (q[2] = lo[2]; q[2] <= hi[2]; ++q[2]){
              auto it = std::lower_bound(first, last, Boxed{level, q, 0}, by_box);
              for (; it != last && it->key == q; ++it){
                const double d = dist(x, &coords[std::size_t(it->v)*gdim]);
                if (d <= rel_tol*h[it->v] && (d < best_d || (d == best_d && it->v < best))){
                  best_d = d;
                  best = it->v;
                }
              }
            }
      }
      vertex_of[p] = best;
    }
  }
  for (const std::size_t j : chosen)
    for (int a = 0; a < D; ++a){
      const std::size_t p = std::size_t(f.facets[j*gdim + std::size_t(a)]);
      if (vertex_of[p] != unmatched) continue;
      std::string at;
      for (int d = 0; d < D; ++d) at += (d ? ", " : "") + std::to_string(f.points[p*gdim + std::size_t(d)]);
      partrac::fail("point ", p, " (", at, ") of ", facet_name(j), " matches no vertex on the mesh's boundary");
    }

  // The exterior facets by their sorted vertices, the periodic sides' among them
  struct Facet { std::array<std::uint32_t, 3> key; std::size_t slot; };
  std::vector<Facet> ext;
  for (std::size_t s = 0; s < ncells*std::size_t(nv); ++s){
    if (facet_neigh[s] >= 0) continue;
    const std::uint32_t* row = topo.data() + (s/std::size_t(nv))*std::size_t(nv);
    Facet e{{0, 0, std::numeric_limits<std::uint32_t>::max()}, s};
    for (int j = 0, m = 0; j < nv; ++j)
      if (j != int(s % std::size_t(nv))) e.key[std::size_t(m++)] = row[j];
    sort_key(e.key, D);
    ext.push_back(e);
  }
  const auto by_facet = [](const Facet& a, const Facet& b){ return a.key < b.key; };
  __gnu_parallel::sort(ext.begin(), ext.end(), by_facet);

  std::vector<std::uint8_t> walls(ncells, 0);
  std::size_t n_walls = 0;
  for (const std::size_t j : chosen){
    Facet q{{0, 0, std::numeric_limits<std::uint32_t>::max()}, 0};
    for (int a = 0; a < D; ++a)
      q.key[std::size_t(a)] = std::uint32_t(vertex_of[std::size_t(f.facets[j*gdim + std::size_t(a)])]);
    sort_key(q.key, D);
    const auto it = std::lower_bound(ext.begin(), ext.end(), q, by_facet);
    if (it == ext.end() || it->key != q.key)
      partrac::fail(facet_name(j), " matches no exterior facet of the mesh");
    if (facet_neigh[it->slot] != facet_wall)
      partrac::fail(facet_name(j), " lies on a periodic side of the mesh, which is no wall");
    std::uint8_t& c = walls[it->slot/std::size_t(nv)];
    const std::uint8_t bit = std::uint8_t(1u << (it->slot % std::size_t(nv)));
    if (!(c & bit)) ++n_walls;
    c |= bit;
  }
  if (verbose){
    const std::size_t n_ext = std::size_t(std::count(facet_neigh.begin(), facet_neigh.end(), facet_wall));
    std::cout << "Wall facets: " << n_walls << " of " << n_ext << " exterior facets off the periodic sides, from "
              << name << std::endl;
  }
  return walls;
}

std::vector<std::uint8_t> load(const std::string& folder, const std::string& file, const std::string& markers,
                               const std::vector<std::uint32_t>& topo, const std::vector<double>& coords,
                               const std::size_t ncells, const Uint gdim,
                               const std::vector<std::int32_t>& facet_neigh, const bool verbose)
{
  const std::string path = std::filesystem::path(file).is_absolute() ? file : folder + "/" + file;
  const Markers sel(markers);
  return match(read(path, gdim), sel, topo, coords, ncells, gdim, facet_neigh, file, verbose);
}

}  // namespace wall_facets
