#include <cmath>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>

#include "RunLoop.hpp"
#include "TimeScheme.hpp"

#include "filaments_schema.hpp"

// Reinject whole edges that have a stuck node
inline void reinject_edges(Run& run, Topology& mesh, ParticleSet& ps, const std::vector<Uint>& nodes){
  std::vector<Uint> edge_ids;
  for (const Uint i : nodes)
    edge_ids.insert(edge_ids.end(), mesh.node2edges[i].begin(), mesh.node2edges[i].end());
  std::sort(edge_ids.begin(), edge_ids.end());
  edge_ids.erase(std::unique(edge_ids.begin(), edge_ids.end()), edge_ids.end());
  const auto key = split_string(run.prm.get<std::string>("init_mode"), "_");
  const std::string dirs = key.size() > 2 ? key[2] : key[1];
  const bool rx = contains(dirs, "x"), ry = contains(dirs, "y"), rz = contains(dirs, "z");
  const Vector3d Dx_max = 0.5*(run.intp->get_x_max() - run.intp->get_x_min());
  std::uniform_real_distribution<> ux(-Dx_max[0], Dx_max[0]), uy(-Dx_max[1], Dx_max[1]), uz(-Dx_max[2], Dx_max[2]);
  for (const Uint e : edge_ids){
    const Uint a = mesh.edges[e].first[0];
    const Uint b = mesh.edges[e].first[1];
    const Vector3d x0 = 0.5*(ps.x(a) + ps.x(b));
    const Vector3d dx = ps.x(a) - ps.x(b);
    Vector3d Dx = {0., 0., 0.};
    bool outside = true;
    while (outside){
      if (rx) Dx[0] = ux(run.gens[0]);
      if (ry) Dx[1] = uy(run.gens[0]);
      if (rz) Dx[2] = uz(run.gens[0]);
      const bool inside_a = run.intp->locate(x0 + Dx + 0.5*dx);
      const bool inside_b = run.intp->locate(x0 + Dx - 0.5*dx);
      outside = !(inside_a && inside_b);
    }
    ps.set_x(a, ps.x(a) + Dx);
    ps.set_x(b, ps.x(b) + Dx);
  }
}

int main(int argc, char* argv[])
{

    std::cout << "Initialized FILAMENTS." << std::endl;

  // Input parameters
  if (argc < 2) {
    std::cout << "Specify an input file." << std::endl;
    return 1;
  }
  partrac::Params prm = partrac::parse_or_exit(filaments_schema(), argc, argv);

  Run run = start_run(prm, "Filaments");
  TimeScheme scheme(prm, run.gens);

  ParticleSet ps(run.intp, prm.get<Uint>("Nrw_max"));
  Topology mesh(ps, prm);
  // Doublings
  const bool doublings = prm.get<std::string>("resize") == "doublings";
  mesh.records_doublings = doublings;

  load_or_initialize(run, mesh);

  std::map<std::string, bool> output_fields;
  output_fields["u"] = !prm.get<bool>("minimal_output");
  output_fields["c"] = !prm.get<bool>("minimal_output");
  output_fields["p"] = !prm.get<bool>("minimal_output") && prm.get<bool>("output_all_props");
  output_fields["rho"] = !prm.get<bool>("minimal_output") && prm.get<bool>("output_all_props");
  // H and n need the curvature
  output_fields["H"] = !prm.get<bool>("minimal_output") && mesh.dim() > 0 && mesh.computes_curvature();
  output_fields["n"] = !prm.get<bool>("minimal_output") && mesh.dim() > 1 && mesh.computes_curvature();

  const double dt = prm.get<double>("dt");
  const double resize_intv = prm.get<double>("resize_intv");
  const double resize_to = prm.get<std::string>("resize_target") == "ds_init"
                         ? prm.get<double>("ds_init") : prm.get<double>("ds_max");
  const std::string outside = prm.get<std::string>("outside");
  const bool verbose = prm.get<bool>("verbose");

  RunHooks hooks;

  // Resizing
  hooks.reshape = [&](const int it, const double){
    if (!at_interval(it, resize_intv, dt))
      return;
    const bool resized = doublings ? mesh.resize_doublings(resize_to) : mesh.resize(resize_to);
    if (resized && verbose)
      std::cout << "Resized edges." << std::endl;
  };

  // Pair statistics with doublings
  if (doublings)
    hooks.statistics = [&](const double t, Integrator& counters){
      return pair_stats_columns(t, ps, mesh.edges, mesh.edge_doublings(), counters.get_declined());
    };

  hooks.after_step = [&](const int, const double, const std::vector<Uint>& outside_nodes){
    if (outside_nodes.size() == 0)
      return;
    std::cout << "Some nodes are outside.\n";
    if (outside == "reinject")
      reinject_edges(run, mesh, ps, outside_nodes);
  };

  run_loop(run, ps, mesh, scheme, output_fields, dt, hooks);

  return 0;
}
