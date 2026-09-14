#include <iostream>
#include <map>
#include <set>
#include <string>

#include "RunLoop.hpp"
#include "SpatialIntegrator.hpp"

#include "static_space_stepper_schema.hpp"

int main(int argc, char* argv[])
{

    std::cout << "Initialized spatial stepper." << std::endl;

  // Input parameters
  if (argc < 2) {
    std::cout << "Specify an input file." << std::endl;
    return 1;
  }
  partrac::Params prm = partrac::parse_or_exit(spatial_schema(), argc, argv);

  // March from xn0 to Ln, fields frozen at t0
  Run run = start_run(prm, "StaticSpaceStepper", DefaultLayout, true);
  SpatialIntegrator integrator(prm.get<int>("int_order"), prm.get<double>("u_eps"),
                               prm.get<double>("dx_max"), prm.get<double>("T"));

  ParticleSet ps(run.intp, prm.get<Uint>("Nrw_max"));
  Topology mesh(ps, prm);
  // Checkpoint t_loc
  mesh.records_t_loc = true;

  load_or_initialize(run, mesh);

  const bool refine = prm.get<bool>("refine");
  const bool coarsen = prm.get<bool>("coarsen");
  const bool verbose = prm.get<bool>("verbose");

  // Initial refinement and coarsening
  if (refine && !prm.get<bool>("inject") && mesh.dim() > 0){
    std::cout << "Initial refinement" << std::endl;
    Uint n_add = mesh.refine();
    if (verbose)
      std::cout << "Added " << n_add << " edges." << std::endl;
  }
  if (coarsen && !prm.get<bool>("inject") && mesh.dim() > 0){
    std::cout << "Initial coarsening" << std::endl;
    Uint n_rem = mesh.coarsen(true);
    if (verbose)
      std::cout << "Removed " << n_rem << " edges." << std::endl;
  }

  mesh.compute_interior();

  std::map<std::string, bool> output_fields;
  output_fields["u"] = !prm.get<bool>("minimal_output");
  output_fields["c"] = true;
  output_fields["p"] = !prm.get<bool>("minimal_output") && prm.get<bool>("output_all_props");
  output_fields["rho"] = !prm.get<bool>("minimal_output") && prm.get<bool>("output_all_props");
  // H and n need the curvature
  output_fields["H"] = !prm.get<bool>("minimal_output") && mesh.dim() > 0 && mesh.computes_curvature();
  output_fields["n"] = !prm.get<bool>("minimal_output") && mesh.dim() > 1 && mesh.computes_curvature();
  output_fields["t_loc"] = true;
  output_fields["tau"] = true;

  // Coarsen at refine_intv when coarsening is off
  const double coarsen_intv = coarsen ? prm.get<double>("coarsen_intv")
                                      : prm.get<double>("refine_intv");

  // Hook constants
  const double refine_intv = prm.get<double>("refine_intv");
  const double dump_intv = prm.get<double>("dump_intv");
  const double dxn = prm.get<double>("dxn");
  const double T_final = prm.get<double>("T");
  const std::string outside = prm.get<std::string>("outside");

  RunHooks hooks;

  // Reshaping
  hooks.reshape = [&](const int it, const double){
    // Curvature computation
    if ((refine && at_interval(it, refine_intv, dxn)) || at_interval(it, coarsen_intv, dxn) || at_interval(it, dump_intv, dxn)){
      mesh.compute_interior();
    }
    // Refinement
    if (refine && at_interval(it, refine_intv, dxn) && it > 0){
      Uint n_add = mesh.refine();
      if (verbose)
        std::cout << "Added " << n_add << " edges." << std::endl;
    }
    // Coarsening
    if (at_interval(it, coarsen_intv, dxn)){
      Uint n_rem = mesh.coarsen(coarsen);
      if (verbose)
        std::cout << "Removed " << n_rem << " edges." << std::endl;
    }
  };

  hooks.after_step = [&](const int, const double xn, const std::vector<Uint>& nodes){
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
      return with_concrete(intp, [&](auto& ip){ return integrator.step(ip, ps, t_fields, ds); });
    }
  } stepper{integrator, run.t_fields};

  run_loop(run, ps, mesh, stepper, output_fields, dxn, hooks);

  return 0;
}
