#ifndef __FILAMENTS_SCHEMA_HPP
#define __FILAMENTS_SCHEMA_HPP

#include "typedefs.hpp"
#include "Params.hpp"
#include "interpol_factory.hpp"
#include "strings.hpp"
#include "TracerParams.hpp"

// Parameters accepted by this app
inline partrac::Schema filaments_schema(){
  partrac::Schema s("filaments");
  s.require<std::string>("mode", "interpolator type");
  s.require<std::string>("init_mode", "initial distribution");
  s.require<double>("Dm", "molecular diffusivity");
  s.require<double>("dt", "timestep");
  s.require<double>("T", "final time");
  s.require<Uint>("Nrw", "number of particles");
  // Particles placed and current
  s.runtime<Uint>("Nrw_init", 0, "particles the initializer placed");
  s.runtime<Uint>("Nrw_current", 0, "particles in the set when this was written");
  s.require<Uint>("Nrw_max", "max number of particles");
  s.require<int>("int_order", "integration order");
  // Required unless resize_target is ds_init
  s.require_if<double>("ds_max", 0.0,
                       [](const partrac::Params& p){ return p.get<std::string>("resize_target") == "ds_max"; },
                       "resize_target is ds_max", "max edge length");
  s.opt<double>("ds_min", 0.0, "min edge length");
  s.require<double>("ds_init", "initial edge length");
  s.opt<double>("t0", 0.0, "start time");
  add_run_params(s, "current time");
  s.opt<double>("x0", 0.0, "initial position");
  s.opt<double>("y0", 0.0, "initial position");
  s.opt<double>("z0", 0.0, "initial position");
  add_intervals(s, "dt", {"resize_intv"});
  s.opt<double>("resize_intv", 0.0, "resize interval");
  s.opt<std::string>("resize", "rescale", "rescale an edge to the target, reference with it (ds/ds0 kept); or doublings: halve it until it fits, reference kept, halvings counted and dumped");
  s.opt<std::string>("resize_target", "ds_max", "the length a resize brings an edge back to: ds_max or ds_init");
  s.opt<std::string>("outside", "ignore", "an edge with a node that cannot take its step: ignore, or reinject the whole edge at a random offset");
  s.opt<int>("sort_every", 0, "reorder particles by cell every this many steps, 0 = never");
  s.opt<double>("t_frozen", 0.0, "time to freeze the fields at");
  s.opt<double>("curv_refine_factor", 0.0, "curvature refinement factor");
  s.opt<int>("filter_target", 0, "filter target");
  s.opt<bool>("filter", false, "filter the filament");
  s.opt<bool>("inject", false, "inject new particles");
  s.opt<bool>("inject_edges", true, "inject edges too");
  s.opt<bool>("frozen_fields", false, "freeze the velocity field");
  s.opt<bool>("local_dt", false, "use a local timestep");
  s.opt<bool>("cut_if_stuck", true, "cut edges that get stuck");
  s.opt<bool>("clear_initial_edges", false, "drop the initial edges");
  s.opt<bool>("output_all_props", true, "dump all properties");
  add_scheme(s, "explicit");
  s.choices("mode", interpol_modes());
  s.choices("resize", {"rescale", "doublings"});
  s.choices("resize_target", {"ds_max", "ds_init"});
  s.choices("outside", {"ignore", "reinject"});
  s.check([](const partrac::Params& p){ return p.get<int>("int_order") <= 2; },
          "int_order must be 1 or 2");
  // Pairs only
  s.check([](const partrac::Params& p){
            const auto key = split_string(p.get<std::string>("init_mode"), "_");
            return !key.empty() && (key[0] == "pair" || key[0] == "pairs");
          },
          "init_mode must be pair_* or pairs_*");
  // No injection
  s.check([](const partrac::Params& p){ return !p.get<bool>("inject"); },
          "filaments does not inject");
  // init_mode needs a direction
  s.check([](const partrac::Params& p){
            return split_string(p.get<std::string>("init_mode"), "_").size() >= 2;
          },
          "init_mode is missing a direction, as in pairs_xyz");
  // No reinjection with diffusion
  s.check([](const partrac::Params& p){
            return !(p.get<std::string>("outside") == "reinject"
                     && p.get<std::string>("scheme") == "explicit"
                     && p.get<double>("Dm") > 0.);
          },
          "outside=reinject is for edges stuck in an underresolved field, not for diffusion: it cannot be combined with scheme=explicit and Dm > 0");
  return s;
}

#endif
