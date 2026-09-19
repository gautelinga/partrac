#include <cmath>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>

#include "Error.hpp"
#include "param_print.hpp"
#include "RunLoop.hpp"
#include "TimeScheme.hpp"

#include "partrac_schema.hpp"

static int run(int argc, char* argv[])
{

    std::cout << "Initialized Partrac." << std::endl;

  // Input parameters
  if (argc < 2) {
    std::cout << "Specify an input file." << std::endl;
    return 1;
  }
  partrac::Params prm = partrac::parse_or_exit(partrac_schema(), argc, argv);

  // Dry run
  const bool dry_run = prm.check_only();

  Run run = start_run(prm, "RandomWalkers", dry_run ? DryRun : DefaultLayout);
  TimeScheme scheme(prm, run.gens);

  ParticleSet ps(run.intp, prm.get<Uint>("Nrw_max"));
  if (prm.get<bool>("output_phi"))
    ps.record_phi();
  Topology mesh(ps, prm);

  if (prm.get<bool>("inject"))
    std::cout << "Injection activated!" << std::endl;

  const bool restarting = load_or_initialize(run, mesh);

  // Parameters checked
  if (dry_run){
      std::cout << "Check OK: " << ps.N() << " particles, dim = " << mesh.dim() << std::endl;
    return 0;
  }

  const bool refine = prm.get<bool>("refine");
  const bool coarsen = prm.get<bool>("coarsen");
  const bool filter = prm.get<bool>("filter");
  const bool inject = prm.get<bool>("inject");
  const bool verbose = prm.get<bool>("verbose");

  // Initial refinement, not on a restart
  if (refine && !restarting && !inject && mesh.dim() > 0){
    std::cout << "Initial refinement" << std::endl;
    Uint n_add = mesh.refine();
    if (verbose)
      std::cout << "Added " << n_add << " edges." << std::endl;
  }
  if (coarsen && !restarting && !inject && mesh.dim() > 0){
    std::cout << "Initial coarsening" << std::endl;
    Uint n_rem = mesh.coarsen(true);
    if (verbose)
      std::cout << "Removed " << n_rem << " edges." << std::endl;
  }

  mesh.compute_interior();

  const double dt = prm.get<double>("dt");
  if (verbose){
    print_param("sqrt(2*Dm*dt)", sqrt(2*prm.get<double>("Dm")*dt));
    print_param("U*dt         ", prm.get<double>("U")*dt);
  }

  std::map<std::string, bool> output_fields;
  output_fields["u"] = !prm.get<bool>("minimal_output");
  output_fields["c"] = !prm.get<bool>("minimal_output");
  output_fields["p"] = !prm.get<bool>("minimal_output") && prm.get<bool>("output_all_props");
  output_fields["rho"] = !prm.get<bool>("minimal_output"); // && prm.output_all_props;
  // H and n need the curvature
  output_fields["H"] = !prm.get<bool>("minimal_output") && mesh.dim() > 0 && mesh.computes_curvature();
  output_fields["n"] = !prm.get<bool>("minimal_output") && mesh.dim() > 1 && mesh.computes_curvature();
  output_fields["tau"] = prm.get<bool>("integrate_tau");
  output_fields["phi"] = prm.get<bool>("output_phi");

  // Coarsen at refine_intv when coarsening is off
  const double coarsen_intv = coarsen ? prm.get<double>("coarsen_intv")
                                      : prm.get<double>("refine_intv");

  // Hook constants
  const bool integrate_tau = prm.get<bool>("integrate_tau");
  const double inject_intv = prm.get<double>("inject_intv");
  const double T_inject = prm.get<double>("T_inject");
  const double refine_intv = prm.get<double>("refine_intv");
  const double filter_intv = prm.get<double>("filter_intv");
  const double dump_intv = prm.get<double>("dump_intv");
  const double tau_intv = prm.get<double>("tau_intv");
  const double tau_max = prm.get<double>("tau_max");
  const double Ln = prm.get<double>("Ln");
  const std::string outside = prm.get<std::string>("outside");

  bool any_exit_plane = (prm.get<std::string>("exit_plane") == "x" || prm.get<std::string>("exit_plane") == "y" || prm.get<std::string>("exit_plane") == "z") && Ln > 0;
  int exit_dim = prm.get<std::string>("exit_plane") == "x" ? 0 : (prm.get<std::string>("exit_plane") == "y" ? 1 : 2);

  RunHooks hooks;

  // Reshaping
  hooks.reshape = [&](const int it, const double t){
    // Injection
    if (inject && it > 0 && at_interval(it, inject_intv, dt) && t <= T_inject){
      mesh.inject();
    }
    // Curvature computation
    if ((refine && at_interval(it, refine_intv, dt)) || at_interval(it, coarsen_intv, dt) || at_interval(it, dump_intv, dt)){
      mesh.compute_interior();
    }
    // Refinement
    if (refine && at_interval(it, refine_intv, dt) && it > 0){
      Uint n_add = mesh.refine();
      if (verbose)
        std::cout << "Added " << n_add << " edges." << std::endl;
    }
    // Coarsening
    if (at_interval(it, coarsen_intv, dt)){
      Uint n_rem = mesh.coarsen(coarsen);
      if (verbose)
        std::cout << "Removed " << n_rem << " edges." << std::endl;
    }
    // Filtering
    if (filter && at_interval(it, filter_intv, dt)){
      bool filtered = mesh.filter();
      if (verbose && filtered)
        std::cout << "Filtered edges." << std::endl;
    }
    // Removal
    if (any_exit_plane && at_interval(it, filter_intv, dt)){
      Uint n_rem = mesh.remove_beyond(exit_dim, Ln);
      if (verbose)
        std::cout << "Removed " << n_rem << " nodes that were beyond." << std::endl;
    }
  };

  hooks.after_step = [&](const int it, const double t, const std::vector<Uint>& outside_nodes){
    // Tau integration
    if (integrate_tau && at_interval(it + 1, tau_intv, dt)){
      mesh.integrate_tau(dt * steps_per(tau_intv, dt), tau_max);
    }
    // Outside nodes
    handle_outside(run, mesh, ps, outside, outside_nodes, t, verbose);
  };

  run_loop(run, ps, mesh, scheme, output_fields, dt, hooks);

  return 0;
}

int main(int argc, char* argv[])
{
  return partrac::report_errors([&]{ return run(argc, argv); });
}
