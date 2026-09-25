#ifndef __WEIGHTED_WALKERS_SCHEMA_HPP
#define __WEIGHTED_WALKERS_SCHEMA_HPP

#include <algorithm>
#include <string>
#include <vector>
#include "typedefs.hpp"
#include "Params.hpp"
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
  s.opt<double>("U", 1.0, "velocity scale");
  s.opt<double>("x0", 0.0, "initial position");
  s.opt<double>("y0", 0.0, "initial position");
  s.opt<double>("z0", 0.0, "initial position");
  s.opt<double>("dump_intv", 100.0, "dump interval");
  s.opt<double>("stat_intv", 100.0, "statistics interval");
  s.opt<double>("checkpoint_intv", 1000.0, "checkpoint interval");
  s.opt<int>("seed", 0, "random seed");
  s.opt<bool>("random", true, "draw the seed randomly");
  s.opt<int>("num_threads", 0, "OpenMP threads, 0 = leave alone");
  s.opt<int>("dump_chunk_size", 0, "particles per dump chunk");
  s.opt<bool>("minimal_output", false, "dump less");
  s.opt<bool>("output_all_props", true, "accepted as before; a walker dumps no more for it");
  s.opt<bool>("verbose", false, "print the parameters");
  s.opt<std::string>("tag", "", "appended to the folder name");
  s.opt<std::string>("restart_folder", "", "folder to restart from");
  s.opt<bool>("inject", false, "accepted as before; walkers are never injected");
  s.opt<bool>("clear_initial_edges", false, "accepted as before; walkers have no edges");
  s.runtime<std::string>("folder", "", "output folder");
  s.runtime<double>("t", 0.0, "current time");
  s.runtime<Uint>("it", 0, "current step");
  s.runtime<double>("Lx", 0.0, "domain size, from the interpolator");
  s.runtime<double>("Ly", 0.0, "domain size, from the interpolator");
  s.runtime<double>("Lz", 0.0, "domain size, from the interpolator");
  s.runtime<Uint>("Nrw_init", 0, "particles the initializer placed");
  s.runtime<Uint>("Nrw_current", 0, "particles in the set when this was written");
  // Cloud: no edges
  s.runtime<double>("ds_min", 0.0, "a cloud: no edges");
  s.runtime<bool>("inject_edges", true, "a cloud: no injection");
  s.runtime<bool>("cut_if_stuck", true, "a cloud: no edges");
  s.runtime<double>("curv_refine_factor", 0.0, "a cloud: no refinement");
  s.runtime<int>("filter_target", 0, "a cloud: no filtering");

  s.choices("mode", {"analytic", "structured", "lbm", "felbm", "fenics",
                     "tet", "triangle", "trianglefreq", "tetfreq", "xdmftriangle", "xdmftet", "openfoam"});
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
  // Intervals: 0 is off, negative is an error
  s.check([](const partrac::Params& p){
            for (const auto& key : {"checkpoint_intv", "dump_intv", "stat_intv", "refine_intv"})
              if (p.get<double>(key) < 0.) return false;
            return true;
          },
          "intervals cannot be negative");
  // Floor output intervals at one step
  s.finalize([](partrac::Params& p){
    const double dt = p.get<double>("dt");
    if (p.get<double>("dump_intv") > 0.)
      p.set<double>("dump_intv", std::max(p.get<double>("dump_intv"), dt));
    if (p.get<double>("stat_intv") > 0.)
      p.set<double>("stat_intv", std::max(p.get<double>("stat_intv"), dt));
    p.set<Uint>("Nrw_max", std::max(p.get<Uint>("Nrw_max"), p.get<Uint>("Nrw")));
  });
  return s;
}

#endif
