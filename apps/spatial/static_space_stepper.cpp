#include <iostream>
#include <set>
#include <string>

#include "Error.hpp"
#include "wait_policy.hpp"
#include "RunLoop.hpp"
#include "SpatialIntegrator.hpp"

#include "static_space_stepper_schema.hpp"

static int run(int argc, char* argv[])
{

    std::cout << "Initialized spatial stepper." << std::endl;

  partrac::Params prm = partrac::parse_or_exit(spatial_schema(), argc, argv);

  // March from xn0 to Ln, fields frozen at t0
  Run run = start_run(prm, "StaticSpaceStepper", DefaultLayout, true);
  SpatialIntegrator integrator(prm.get<int>("int_order"), prm.get<double>("u_eps"),
                               prm.get<double>("dx_max"), prm.get<double>("T"));

  ParticleSet ps(run.intp, prm.get<Uint>("Nrw_max"));
  Topology mesh(ps, mesh_options(prm));
  // Checkpoint t_loc
  mesh.records_t_loc = true;

  const bool restarting = load_or_initialize(run, mesh);
  if (check_only(run, ps, mesh))
    return 0;

  const bool verbose = prm.get<bool>("verbose");
  const double dxn = prm.get<double>("dxn");
  Remesher remesher(mesh, prm, dxn);
  remesher.initial(restarting, prm.get<bool>("inject"));

  OutputFields output_fields = mesh_output_fields(prm, mesh);
  output_fields.c = true;
  output_fields.t_loc = true;
  output_fields.tau = true;

  // Hook constants
  const double T_final = prm.get<double>("T");
  const std::string outside = prm.get<std::string>("outside");

  RunHooks hooks;

  // Reshaping
  hooks.reshape = [&](const int it, const double){ remesher.reshape(it); };

  hooks.outside = [&](const std::vector<Uint>& nodes, const double xn){
    // Finished and trapped nodes
    if (verbose && nodes.size() > 0){
      Vector3d x_trapped = {0., 0., 0.};
      Uint n_done = 0, n_trapped = 0;
      for (const Uint i : nodes){
        if (ps.t_loc(i) >= T_final){
          ++n_done;
        }
        else {
          x_trapped += ps.x(i);
          ++n_trapped;
        }
      }
      std::cout << "At xn = " << xn << ": " << n_done
                << " nodes finished their integration time";
      if (n_trapped > 0){
        x_trapped /= n_trapped;
        std::cout << ", " << n_trapped << " could not move, centred on ("
                  << x_trapped[0] << ", " << x_trapped[1] << ", "
                  << x_trapped[2] << ")";
      }
      std::cout << std::endl;
    }
    handle_outside(run, mesh, ps, outside, nodes, xn, false);
  };

  // Stop when no node is left
  hooks.keep_going = [&]{ return ps.N() > 0; };

  // Stepper
  struct Stepper {
    SpatialIntegrator& integrator;
    double t_fields;
    Integrator& counters(){ return integrator; }
    std::vector<Uint> step(Interpol& intp, ParticleSet& ps, const double, const double ds){
      return spatial_step<TransportElement::Point>(integrator, intp, ps, t_fields, ds);
    }
  } stepper{integrator, run.t_fields};

  OutsideTotals totals{outside, *run.intp, [&]{ return integrator.get_left(); }};
  totals.attach(hooks, verbose);

  run_loop(run, ps, mesh, stepper, output_fields, dxn, hooks);
  totals.finish(std::cout);

  return 0;
}

int main(int argc, char* argv[])
{
  partrac::wait_passively(argv);
  return partrac::report_errors([&]{ return run(argc, argv); });
}
