#ifndef __PARTRAC_SCHEMA_HPP
#define __PARTRAC_SCHEMA_HPP

#include "typedefs.hpp"
#include "Params.hpp"
#include "interpol_factory.hpp"
#include "Initializer.hpp"
#include "TracerParams.hpp"

// Parameters accepted by this app
inline partrac::Schema partrac_schema(){
  partrac::Schema s("partrac");
  add_initializer_params(s);
  s.require<std::string>("mode", "interpolator type");
  s.opt<std::string>("outside", "ignore", "a particle that cannot take its step: ignore (it stays), reinject at a random offset, or mark (c = 2)");
  s.opt<int>("sort_every", 5, "reorder particles by cell every this many steps, 0 = never");
  s.opt<bool>("output_phi", false, "dump the phase field at each particle, and split the vector statistics on it");
  s.require<double>("Dm", "molecular diffusivity");
  s.require<double>("dt", "timestep");
  s.require<double>("T", "final time");
  s.require<int>("int_order", "integration order");
  // exit_plane cuts at Ln
  s.require_if<double>("Ln", 0.0,
                       [](const partrac::Params& p){
                         return p.get<std::string>("exit_plane") != "none";
                       },
                       "exit_plane is set", "exit plane position");
  s.require_if<double>("filter_intv", 0.0,
                       [](const partrac::Params& p){
                         return p.get<std::string>("exit_plane") != "none" ||
                                p.get<bool>("filter");
                       },
                       "exit_plane or filter is set", "filter interval");
  s.require_if<double>("inject_intv", 0.0,
                       [](const partrac::Params& p){ return p.get<bool>("inject"); },
                       "inject is set", "injection interval");
  s.require_if<double>("tau_intv", 0.0,
                       [](const partrac::Params& p){ return p.get<bool>("integrate_tau"); },
                       "integrate_tau is set", "tau interval");
  add_run_params(s, "current time");
  s.opt<double>("T_inject", 1e10, "time to stop injecting at");
  s.opt<double>("tau_max", 0.0, "max tau");
  s.optional<double>("t_frozen", "time to freeze the fields at, clamped to the fields' times; default t0");
  add_intervals(s, "dt", {"coarsen_intv", "filter_intv", "inject_intv", "refine_intv", "tau_intv"});
  s.opt<double>("refine_intv", 100.0, "refinement interval");
  s.opt<double>("coarsen_intv", 1000.0, "coarsening interval");
  s.opt<double>("curv_refine_factor", 0.0, "curvature refinement factor");
  s.opt<int>("filter_target", 0, "filter target");
  s.opt<bool>("refine", false, "refine the mesh");
  s.opt<bool>("coarsen", false, "coarsen the mesh");
  s.opt<bool>("filter", false, "filter the mesh");
  s.opt<bool>("inject_edges", true, "inject edges too");
  s.opt<bool>("frozen_fields", false, "freeze the velocity field");
  s.opt<bool>("cut_if_stuck", true, "cut edges that get stuck");
  s.opt<bool>("integrate_tau", false, "integrate the eigentime");
  s.opt<bool>("output_all_props", true, "dump all properties");
  add_scheme(s, "explicit");
  s.opt<std::string>("exit_plane", "none", "plane to remove particles beyond");
  s.choices("mode", interpol_modes());
  s.choices("outside", {"ignore", "reinject", "mark"});
  s.choices("exit_plane", {"none", "x", "y", "z"});
  s.check([](const partrac::Params& p){ return p.get<int>("int_order") <= 2; },
          "int_order must be 1 or 2");
  s.check([](const partrac::Params& p){
            return !(p.get<bool>("inject") && p.get<bool>("filter"));
          },
          "cannot inject and filter at the same time");
  return s;
}

#endif
