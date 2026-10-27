#include "SimplexFreqInterpol.hpp"
#include "loader_params.hpp"
#include "mesh_tables.hpp"
#include "p12_eval.hpp"
#include "phase_timing.hpp"
#include "simplex_load.hpp"
#include <array>
#include <cmath>
#include <iostream>

namespace {

// The components are one function space written many times: a later file that
// holds another element is not summable with the first
void same_element(const simplex_load::Element& a, const simplex_load::Element& b,
                  const std::string& path, const std::string& field){
  if (a.degree != b.degree || a.ncomp != b.ncomp)
    partrac::fail(path, ": '", field, "' is of degree ", b.degree, " with ", b.ncomp,
                  " components, the first component's degree ", a.degree, " with ", a.ncomp);
}

}  // namespace

template<typename Cell>
SimplexFreqInterpol<Cell>::SimplexFreqInterpol(const std::string& infilename)
  : MeshCore<Cell>(infilename)
{
  constexpr int nv = Cell::n_verts;
  new_id();
  partrac::phase_begin("load");
  dolfin_params = partrac::parse_file_or_exit(simplex_freq_schema(mode), infilename);

  const std::size_t botDirPos = infilename.find_last_of("/");
  set_folder(infilename.substr(0, botDirPos));

  // Using the FreqStamps class to hold frequency data
  fs.initialize(get_folder() + "/" + dolfin_params.template get<std::string>("freqstamps"));
  if (fs.size() < 1)
    partrac::fail(infilename, ": the frequency stamps list no component");

  // base frequency, for lines that give a harmonic by their number
  if (fs.omega_given() == dolfin_params.has("tau"))
    partrac::fail(infilename, fs.omega_given() ? ": tau is not read when the frequency stamps give omega"
                                               : ": tau is needed when the frequency stamps give no omega");
  if (!fs.omega_given()){
    const double tau = dolfin_params.template get<double>("tau");
    if (tau > 0)
      omega0 = 2 * M_PI / tau;
  }

  read_mesh_params();
  partrac::phase("params");

  simplex_load::Request req;
  req.infilename = infilename;
  // The element comes from the first component
  req.field_file = get_folder() + "/" + fs.get(0).filename;
  req.what = "SimplexFreqInterpol";
  req.nv = nv;
  req.include_pressure = include_pressure;
  req.n_dofs_max = Cell::n_dofs_max;
  req.periodic = periodic;
  req.periodic_tol = periodic_tol;
  simplex_load::request_from_params(req, dolfin_params, get_folder());
  const std::string& u_field = req.u_field;
  const std::string& p_field = req.p_field;
  simplex_load::Tables t;
  simplex_load::build_tables(req, t);

  this->adopt_tables(t);
  ncells_ = t.mesh.ncells;
  nverts_ = t.mesh.nverts;

  // Every component through one dof -> nodes mapping: the components are one
  // function space, so only the first has to be read with its dof table
  const std::size_t nfreq = std::size_t(fs.size());
  u_nodes_.resize(nfreq);
  if (include_pressure)
    p_nodes_.resize(nfreq);
  mesh_tables::DofNodes u_map, p_map;
  for (std::size_t iFreq = 0; iFreq < nfreq; ++iFreq){
    const std::string path = get_folder() + "/" + fs.get(int(iFreq)).filename;
    if (iFreq == 0){
      simplex_load::read_field_by_node(path, u_field, t, t.el_u, t.node_order(t.el_u),
                                       u_nodes_[0], u_map);
      if (include_pressure)
        simplex_load::read_field_by_node(path, p_field, t, t.el_p, t.node_order(t.el_p),
                                         p_nodes_[0], p_map);
    }
    else {
      same_element(t.el_u, simplex_load::read_element(path, u_field, nv), path, u_field);
      simplex_load::read_into(path, u_field, u_map, u_nodes_[0], u_nodes_[iFreq]);
      if (include_pressure){
        same_element(t.el_p, simplex_load::read_element(path, p_field, nv), path, p_field);
        simplex_load::read_into(path, p_field, p_map, p_nodes_[0], p_nodes_[iFreq]);
      }
    }
  }
  partrac::phase("component values");

  topo_ = std::move(t.mesh.topo);
  coords_ = std::move(t.mesh.coords);

  // Near walls: the velocity's order raised or not
  wall_upgrade_ = dolfin_params.template get<bool>("wall_upgrade");
  if (wall_upgrade_){
    upgrade_p2_ = near_wall::velocity_p2<Cell>(ncoeffs_u);
    near_wall::build_tables<Cell>(topo_, coords_, ncells_, nverts_, dim, periodic, x_min, x_max,
                                  periodic_tol, facet_neigh_, upgrade_p2_, vclass_, wall_index_,
                                  wall_cells_, upgrade_, verbose);
    rest_masks();
    partrac::phase("wall edges");
  }

  mesh_tables::build_cells_and_tree<Cell>(topo_, coords_, ncells_, nverts_, dim, cells_, tree_);
  partrac::phase_total();

  std::cout << "Setting max threads: " << omp_get_max_threads() << std::endl;
}

template<typename Cell>
void SimplexFreqInterpol<Cell>::update(const double t)
{
  is_initialized = true;
  t_update = t;
}

template<typename Cell>
void SimplexFreqInterpol<Cell>::rest_masks()
{
  // One tolerance for the sum: round-off of the largest value it can reach
  double tol = 0.;
  for (std::size_t k = 0; k < u_nodes_.size(); ++k)
    tol += std::abs(fs.get(int(k)).a)*near_wall::rest_tolerance(u_nodes_[k]);
  rest_.assign(upgrade_p2_ ? upgrade_.cells.size() : wall_cells_.size(), 0);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < ncells_; ++i){
    const std::int32_t w = wall_index_[i];
    if (w < 0) continue;
    unsigned rest = upgrade_p2_ ? upgrade_.cells[std::size_t(w)].wall : wall_cells_[std::size_t(w)].wall;
    for (std::size_t k = 0; k < u_nodes_.size() && rest; ++k){
      std::array<double, Cell::n_dofs_max*3> u;
      gather_cell_nodes<Cell::n_verts, Cell::n_dofs_max, D>(u_dofs_[i], u_dofs_.stride(), u_nodes_[k].data(),
                                                            u.data());
      rest = near_wall::rest_vertices<Cell>(u.data(), ncoeffs_u, std::uint8_t(rest), tol);
    }
    rest_[std::size_t(w)] = std::uint8_t(rest);
    // None at rest: read as without the rule
    if (!rest) wall_index_[i] = -1;
  }
}

template class SimplexFreqInterpol<Triangle>;
template class SimplexFreqInterpol<Tet>;
