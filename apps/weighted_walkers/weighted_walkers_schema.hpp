#ifndef __WEIGHTED_WALKERS_SCHEMA_HPP
#define __WEIGHTED_WALKERS_SCHEMA_HPP

#include "typedefs.hpp"
#include "Params.hpp"
#include "experimental/initializer.hpp"

// Parameters accepted by this app
inline partrac::Schema weighted_walkers_schema(){
  partrac::Schema s("weighted_walkers");
  add_common_app_params(s);
  add_experimental_initializer_params(s);
  add_restart_params(s);
  s.require<int>("int_order", "integration order");
  s.opt<int>("num_threads", 0, "OpenMP threads, 0 = leave alone");
  s.require<double>("ds_max", "max edge length");
  s.opt<double>("Lt", 0.0, "tangential extent of the exit plane");
  s.require<double>("La", "principal extent");
  s.require<double>("Lb", "gaussian width");
  s.opt<double>("Ln", 0.0, "exit plane position");
  s.opt<double>("refine_intv", 100.0, "refinement interval");
  s.opt<std::string>("exit_plane", "none", "plane to remove particles beyond");
  s.choices("exit_plane", {"none", "x", "y", "z"});
  // an interval of 0 turns that output off; a negative one is a typo
  s.check([](const partrac::Params& p){
            for (const auto& key : {"refine_intv"})
              if (p.get<double>(key) < 0.) return false;
            return true;
          },
          "intervals cannot be negative");
  return s;
}

#endif
