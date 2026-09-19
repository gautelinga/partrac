#ifndef __TRACERAPP_HPP
#define __TRACERAPP_HPP

// Tracer apps: points, line elements or tensors on the shared run loop

#include <iostream>
#include <map>
#include <set>
#include <string>

#include "RunLoop.hpp"
#include "TimeScheme.hpp"
#include "SpatialIntegrator.hpp"
#include "TransportElement.hpp"
#include "stats.hpp"

// Tracer fields
template<TransportElement E>
void record_tracer_fields(Run& run, ParticleSet& ps, const bool phase_statistics){
  const partrac::Params& prm = run.prm;
  ps.carry(E);
  // Carried elements need the velocity gradient
  if (E != TransportElement::Point)
    run.intp->set_needs_gradient(true);
  if (prm.get<bool>("output_J"))
    ps.record_J();
  // Vector statistics split on phi (0 without a phase field)
  if (phase_statistics || prm.get<bool>("output_phi"))
    ps.record_phi();
  if (prm.get<bool>("output_cell_type"))
    ps.record_cell_type();
  // Dump rhohat as n
  if (E == TransportElement::Vector)
    ps.dump_as("rhohat", "n");
}

// Load checkpoint or initialize; returns whether restarted
template<TransportElement E>
bool start_tracers(Run& run, ParticleSet& ps, Topology& mesh){
  const bool restarting = load_or_initialize(run, mesh);
  // Random slot order on a fresh run (the points initializer sorts by x)
  if (!restarting)
    mesh.shuffle(run.gens[0]);
  // Random initial direction of the line element
  if (E == TransportElement::Vector && !restarting)
    ps.spin_rhohat(split_string(run.prm.get<std::string>("init_mode"), "_").back(), run.gens[0]);
  return restarting;
}

inline std::map<std::string, bool> tracer_output_fields(const partrac::Params& prm){
  std::map<std::string, bool> output_fields;
  output_fields["u"] = !prm.get<bool>("minimal_output");
  output_fields["c"] = !prm.get<bool>("minimal_output");
  output_fields["p"] = !prm.get<bool>("minimal_output") && prm.get<bool>("output_all_props");
  output_fields["rho"] = !prm.get<bool>("minimal_output") && prm.get<bool>("output_all_props");
  output_fields["J"] = prm.get<bool>("output_J");
  output_fields["phi"] = prm.get<bool>("output_phi");
  output_fields["cell_type"] = prm.get<bool>("output_cell_type");
  return output_fields;
}

template<TransportElement E>
int run_tracers(partrac::Params& prm, const std::string& folder){
  Run run = start_run(prm, folder);
  TimeScheme scheme(prm, run.gens);

  ParticleSet ps(run.intp, prm.get<Uint>("Nrw_max"));
  record_tracer_fields<E>(run, ps, E == TransportElement::Vector);
  Topology mesh(ps, prm);
  start_tracers<E>(run, ps, mesh);

  std::map<std::string, bool> output_fields = tracer_output_fields(prm);

  const double dt = prm.get<double>("dt");
  const std::string outside = prm.get<std::string>("outside");
  const bool verbose = prm.get<bool>("verbose");

  // Stepper
  struct Stepper {
    TimeScheme& scheme;
    Integrator& counters(){ return scheme.counters(); }
    std::vector<Uint> step(Interpol& intp, ParticleSet& ps, const double t, const double dt){
      return scheme.template step<E>(intp, ps, t, dt);
    }
  } stepper{scheme};

  RunHooks hooks;
  hooks.statistics = [&](const double t, Integrator& counters){
    return E == TransportElement::Vector
         ? vector_stats_columns(t, ps, true, counters.get_declined())
         : cloud_stats_columns(t, ps, counters.get_declined());
  };
  hooks.after_step = [&](const int, const double t, const std::vector<Uint>& outside_nodes){
    handle_outside(run, mesh, ps, outside, outside_nodes, t, verbose);
  };

  run_loop(run, ps, mesh, stepper, output_fields, dt, hooks);
  return 0;
}

// March in path length, fields frozen
template<TransportElement E>
int run_spatial_tracers(partrac::Params& prm, const std::string& folder){
  Run run = start_run(prm, folder, DefaultLayout, true);
  SpatialIntegrator integrator(prm.get<int>("int_order"), prm.get<double>("u_eps"),
                               prm.get<double>("dx_max"), prm.get<double>("T"));

  ParticleSet ps(run.intp, prm.get<Uint>("Nrw_max"));
  record_tracer_fields<E>(run, ps, false);
  ps.dump_as("t_loc", "tau");
  Topology mesh(ps, prm);
  // Checkpoint t_loc
  mesh.records_t_loc = true;
  start_tracers<E>(run, ps, mesh);

  std::map<std::string, bool> output_fields = tracer_output_fields(prm);
  output_fields["t_loc"] = true;

  const double dxn = prm.get<double>("dxn");
  const std::string outside = prm.get<std::string>("outside");
  const bool verbose = prm.get<bool>("verbose");

  // Stepper
  struct Stepper {
    SpatialIntegrator& integrator;
    double t_fields;
    Integrator& counters(){ return integrator; }
    std::vector<Uint> step(Interpol& intp, ParticleSet& ps, const double, const double ds){
      return spatial_step<E>(integrator, intp, ps, t_fields, ds);
    }
  } stepper{integrator, run.t_fields};

  RunHooks hooks;
  hooks.statistics = [&](const double xn, Integrator& counters){
    return cloud_stats_columns(xn, ps, counters.get_declined());
  };
  hooks.after_step = [&](const int, const double xn, const std::vector<Uint>& nodes){
    handle_outside(run, mesh, ps, outside, nodes, xn, verbose);
  };
  // Stop when no particle is left
  hooks.keep_going = [&]{ return ps.N() > 0; };

  run_loop(run, ps, mesh, stepper, output_fields, dxn, hooks);
  return 0;
}

#endif
