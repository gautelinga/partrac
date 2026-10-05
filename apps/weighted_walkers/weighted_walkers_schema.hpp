#ifndef __WEIGHTED_WALKERS_SCHEMA_HPP
#define __WEIGHTED_WALKERS_SCHEMA_HPP

#include <algorithm>
#include <string>
#include <vector>
#include "typedefs.hpp"
#include "Params.hpp"
#include "interpol_factory.hpp"
#include "AppParams.hpp"
#include "strings.hpp"

// Parameters accepted by this app
inline partrac::Schema weighted_walkers_schema(){
  partrac::Schema s("weighted_walkers");
  s.opt<std::string>("mode", "analytic", "interpolator type");
  s.require<std::string>("init_mode", "initial distribution: strip_<along>_<spread> or circle_<normal>_<spread>, as strip_y_x");
  s.require<Uint>("Nrw", "number of particles");
  s.require<Uint>("Nrw_max", "max number of particles");
  s.require<double>("Dm", "molecular diffusivity");
  s.require<double>("dt", "timestep");
  s.require<double>("T", "final time");
  s.require<int>("int_order", "integration order");
  s.require<double>("La", "principal extent");
  s.require<double>("Lb", "gaussian width");
  s.require<double>("ds_max", "distance from the exit plane's axis within which walkers are written as separation data");
  s.opt<double>("Ln", 0.0, "exit plane position");
  s.opt<double>("Lt", 0.0, "tangential extent of the exit plane");
  s.opt<double>("refine_intv", 100.0, "resampling interval");
  s.opt<std::string>("exit_plane", "none", "plane to remove particles beyond");
  s.opt<double>("t0", 0.0, "start time");
  s.opt<bool>("frozen_fields", false, "freeze the velocity field");
  s.opt<double>("t_frozen", 0.0, "time to freeze the fields at");
  add_run_params(s, "current time");
  s.opt<double>("x0", 0.0, "initial position");
  s.opt<double>("y0", 0.0, "initial position");
  s.opt<double>("z0", 0.0, "initial position");
  add_intervals(s, "dt", {"refine_intv"});
  s.opt<bool>("output_all_props", true, "accepted as before; a walker dumps no more for it");
  s.opt<bool>("inject", false, "accepted as before; walkers are never injected");
  s.opt<bool>("clear_initial_edges", false, "accepted as before; walkers have no edges");
  s.runtime<Uint>("Nrw_init", 0, "particles the initializer placed");
  s.runtime<Uint>("Nrw_current", 0, "particles in the set when this was written");
  // Cloud: no edges
  s.runtime<double>("ds_min", 0.0, "a cloud: no edges");
  // Written by older runs; a cloud's Topology takes none of them
  s.retired({"inject_edges", "cut_if_stuck", "curv_refine_factor", "filter_target"});

  s.choices("mode", interpol_modes());
  s.choices("exit_plane", {"none", "x", "y", "z"});
  // Strip or circle, two direction tokens
  s.check([](const partrac::Params& p){
            const std::vector<std::string> key = split_string(p.get<std::string>("init_mode"), "_");
            if (key.size() != 3) return false;
            if (!contains(key[0], "strip") && !contains(key[0], "circle")) return false;
            for (Uint i = 1; i < 3; ++i){
              if (key[i].empty()) return false;
              for (const char c : key[i])
                if (c != 'x' && c != 'y' && c != 'z') return false;
            }
            return true;
          },
          "init_mode must be a strip or a circle with two direction tokens, as strip_y_x or circle_x_yz");
  s.check([](const partrac::Params& p){ return p.get<int>("int_order") <= 2; },
          "int_order must be 1 or 2");
  s.check([](const partrac::Params& p){ return !p.get<bool>("inject"); },
          "weighted_walkers does not inject");
  return s;
}

#endif
