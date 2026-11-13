#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <set>
#include <string>

#include "Error.hpp"
#include "wait_policy.hpp"
#include "param_print.hpp"
#include "RunLoop.hpp"
#include "TimeScheme.hpp"

#include "partrac_schema.hpp"
#include "crossings.hpp"

static int run(int argc, char* argv[])
{

    std::cout << "Initialized Partrac." << std::endl;

  partrac::Params prm = partrac::parse_or_exit(partrac_schema(), argc, argv);

  Run run = start_run(prm, "RandomWalkers");
  TimeScheme scheme(prm, run.gens);

  ParticleSet ps(run.intp, prm.get<Uint>("Nrw_max"));
  scheme.prepare(*run.intp);
  if (prm.get<bool>("output_phi"))
    ps.record_phi();
  Topology mesh(ps, mesh_options(prm));

  if (prm.get<bool>("inject"))
    std::cout << "Injection activated!" << std::endl;

  const bool restarting = load_or_initialize(run, mesh);

  const bool filter = prm.get<bool>("filter");
  if (filter && !mesh.faces.empty())
    partrac::fail("filter=true filters a strip; a sheet is not supported");
  if (prm.get<double>("cross_dn") > 0. && !mesh.faces.empty())
    partrac::fail("cross_dn records strips and particles; a sheet is not supported");

  if (check_only(run, ps, mesh))
    return 0;

  const bool inject = prm.get<bool>("inject");
  const bool verbose = prm.get<bool>("verbose");

  const double dt = prm.get<double>("dt");
  Remesher remesher(mesh, run, dt);
  remesher.initial(restarting, inject);

  if (verbose){
    print_param("sqrt(2*Dm*dt)", sqrt(2*prm.get<double>("Dm")*dt));
    print_param("U*dt         ", prm.get<double>("U")*dt);
  }

  OutputFields output_fields = mesh_output_fields(prm, mesh);
  // rho whatever output_all_props says
  output_fields.rho = !prm.get<bool>("minimal_output");
  output_fields.tau = prm.get<bool>("integrate_tau");
  output_fields.phi = prm.get<bool>("output_phi");

  // Hook constants
  const bool integrate_tau = prm.get<bool>("integrate_tau");
  const double inject_intv = prm.get<double>("inject_intv");
  const double T_inject = prm.get<double>("T_inject");
  const double filter_intv = prm.get<double>("filter_intv");
  const double tau_intv = prm.get<double>("tau_intv");
  const double tau_max = prm.get<double>("tau_max");
  const double Ln = prm.get<double>("Ln");
  const std::string outside = prm.get<std::string>("outside");

  const std::string exit_plane = prm.get<std::string>("exit_plane");
  const bool any_exit_plane = exit_plane != "none";
  // One normal for the exit plane and the crossing planes: an axis, or cross_n
  const std::string cross_n = prm.get<std::string>("cross_n");
  const Vector3d plane_n = parse_normal(cross_n.empty() ? exit_plane : cross_n);
  std::unique_ptr<Crossings> crossings;
  if (prm.get<double>("cross_dn") > 0.)
    crossings = std::make_unique<Crossings>(ps, !mesh.edges.empty(), plane_n, prm.get<double>("cross_dn"),
                                            run.out.run + "/" + prm.get<std::string>("cross_file"), restarting);

  RunHooks hooks;

  // Reshaping
  hooks.reshape = [&](const int it, const double t){
    // Injection
    if (inject && it > 0 && at_interval(it, inject_intv, dt) && t <= T_inject){
      mesh.inject();
    }
    remesher.reshape(it);
    // Filtering
    if (filter && at_interval(it, filter_intv, dt)){
      const Uint before = mesh.edges.size();
      bool filtered = mesh.filter();
      if (crossings)
        crossings->filtered(t, before, mesh.edges.size());
      if (verbose && filtered)
        std::cout << "Filtered edges." << std::endl;
    }
    // Removal
    if (any_exit_plane && at_interval(it, filter_intv, dt)){
      Uint n_rem = mesh.remove_beyond(plane_n, Ln);
      if (verbose && n_rem > 0)
        std::cout << "Removed " << n_rem << " nodes that were beyond." << std::endl;
    }
  };

  // Outside nodes
  hooks.outside = [&](const std::vector<Uint>& outside_nodes, const double t){
    handle_outside(mesh, ps, outside, outside_nodes, t, verbose);
  };
  // Tau integration
  hooks.after_step = [&](const int it, const double t){
    if (integrate_tau && at_interval(it + 1, tau_intv, dt)){
      mesh.integrate_tau(dt * steps_per(tau_intv, dt), tau_max);
    }
    if (crossings)
      crossings->record(ps, mesh.edges, mesh.node2edges, t, dt);
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
