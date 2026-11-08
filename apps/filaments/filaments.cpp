#include <iostream>
#include <string>

#include "Error.hpp"
#include "wait_policy.hpp"
#include "RunLoop.hpp"
#include "TimeScheme.hpp"

#include "filaments_schema.hpp"

static int run(int argc, char* argv[])
{

    std::cout << "Initialized FILAMENTS." << std::endl;

  partrac::Params prm = partrac::parse_or_exit(filaments_schema(), argc, argv);

  Run run = start_run(prm, "Filaments");
  TimeScheme scheme(prm, run.gens);

  ParticleSet ps(run.intp, prm.get<Uint>("Nrw_max"));
  scheme.prepare(*run.intp);
  Topology mesh(ps, mesh_options(prm));
  // Doublings
  const bool doublings = prm.get<std::string>("resize") == "doublings";
  mesh.records_doublings = doublings;

  load_or_initialize(run, mesh);
  if (check_only(run, ps, mesh))
    return 0;

  const OutputFields output_fields = mesh_output_fields(prm, mesh);

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

  hooks.outside = [&](const std::vector<Uint>& outside_nodes, const double t){
    handle_outside(mesh, ps, outside, outside_nodes, t, verbose);
  };
  OutsideTotals totals{outside, *run.intp, [&]{ return scheme.counters().get_declined(); }};
  totals.attach(hooks, verbose);

  run_loop(run, ps, mesh, scheme, output_fields, dt, hooks);
  scheme.report(std::cout);
  totals.finish(std::cout);

  return 0;
}

int main(int argc, char* argv[])
{
  partrac::wait_passively(argv);
  return partrac::report_errors([&]{ return run(argc, argv); });
}
