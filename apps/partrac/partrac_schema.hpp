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
  add_outside(s, "ignore", "a particle that cannot take its step", false);
  s.opt<int>("sort_every", 5, "reorder particles by cell every this many steps, 0 = never");
  s.opt<bool>("output_phi", false, "dump the phase field at each particle, and split the vector statistics on it");
  s.require<double>("Dm", "molecular diffusivity");
  s.require<double>("dt", "timestep");
  s.require<double>("T", "final time");
  add_int_order(s, "integration order");
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
  add_curv_limit(s);
  s.opt<Uint>("filter_target", 0, "filter target");
  s.opt<bool>("refine", false, "refine the mesh");
  s.opt<bool>("coarsen", false, "coarsen the mesh");
  s.opt<bool>("filter", false, "filter a strip's edges");
  s.opt<bool>("inject_edges", true, "inject edges too");
  s.opt<bool>("frozen_fields", false, "freeze the velocity field");
  s.opt<bool>("cut_if_stuck", true, "cut edges that get stuck");
  s.opt<bool>("integrate_tau", false, "integrate the eigentime");
  s.opt<bool>("output_all_props", true, "dump all properties");
  add_scheme(s, "explicit");
  s.opt<std::string>("exit_plane", "none", "remove nodes with x.n > Ln: n the x, y or z axis, or cross_n with exit_plane=n");
  add_crossing_params(s, false);
  s.choices("mode", interpol_modes());
  s.choices("exit_plane", {"none", "x", "y", "z", "n"});
  s.check([](const partrac::Params& p){
            return !(p.get<bool>("inject") && p.get<bool>("filter"));
          },
          "cannot inject and filter at the same time");
  // One normal: exit_plane=x|y|z, or cross_n
  s.check([](const partrac::Params& p){
            const std::string e = p.get<std::string>("exit_plane");
            return e == "none" || e == "n" || p.get<std::string>("cross_n").empty();
          },
          "exit_plane=x|y|z names the normal; cross_n goes with exit_plane=n");
  s.check([](const partrac::Params& p){
            const std::string e = p.get<std::string>("exit_plane");
            const bool needs = e == "n" || (e == "none" && crossings_on(p));
            return !needs || !p.get<std::string>("cross_n").empty();
          },
          "exit_plane=n, or cross_intv without exit_plane=x|y|z, needs cross_n");
  return s;
}

#endif
