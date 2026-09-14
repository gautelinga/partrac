#ifndef __FILAMENTS_SCHEMA_HPP
#define __FILAMENTS_SCHEMA_HPP

#include "typedefs.hpp"
#include "Params.hpp"
#include "utils.hpp"

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
  // Sampling weight for points
  s.require_if<std::string>("init_weight",
                            [](const partrac::Params& p){
                              return p.get<std::string>("init_mode").rfind("points", 0) == 0;
                            },
                            "init_mode is a points distribution",
                            "sampling weight");
  s.opt<double>("t0", 0.0, "start time");
  s.opt<double>("U", 1.0, "velocity scale");
  s.opt<double>("x0", 0.0, "initial position");
  s.opt<double>("y0", 0.0, "initial position");
  s.opt<double>("z0", 0.0, "initial position");
  s.opt<double>("dump_intv", 100.0, "dump interval");
  s.opt<double>("stat_intv", 100.0, "statistics interval");
  s.opt<double>("checkpoint_intv", 1000.0, "checkpoint interval");
  s.opt<double>("resize_intv", 0.0, "resize interval");
  s.opt<std::string>("resize", "rescale", "rescale an edge to the target, reference with it (ds/ds0 kept); or doublings: halve it until it fits, reference kept, halvings counted and dumped");
  s.opt<std::string>("resize_target", "ds_max", "the length a resize brings an edge back to: ds_max or ds_init");
  s.opt<std::string>("outside", "ignore", "an edge with a node that cannot take its step: ignore, or reinject the whole edge at a random offset");
  s.opt<int>("sort_every", 0, "reorder particles by cell every this many steps, 0 = never");
  s.opt<double>("t_frozen", 0.0, "time to freeze the fields at");
  s.opt<double>("curv_refine_factor", 0.0, "curvature refinement factor");
  s.opt<int>("seed", 0, "random seed");
  s.opt<int>("num_threads", 0, "OpenMP threads, 0 = leave alone");
  s.opt<int>("dump_chunk_size", 0, "particles per dump chunk");
  s.opt<int>("filter_target", 0, "filter target");
  s.opt<bool>("verbose", false, "print the parameters");
  s.opt<bool>("random", true, "draw the seed randomly");
  s.opt<bool>("filter", false, "filter the filament");
  s.opt<bool>("inject", false, "inject new particles");
  s.opt<bool>("inject_edges", true, "inject edges too");
  s.opt<bool>("frozen_fields", false, "freeze the velocity field");
  s.opt<bool>("local_dt", false, "use a local timestep");
  s.opt<bool>("cut_if_stuck", true, "cut edges that get stuck");
  s.opt<bool>("clear_initial_edges", false, "drop the initial edges");
  s.opt<bool>("output_all_props", true, "dump all properties");
  s.opt<bool>("minimal_output", false, "dump less");
  s.opt<std::string>("scheme", "explicit", "ODE integration scheme");
  s.opt<std::string>("tag", "", "appended to the folder name");
  s.opt<std::string>("restart_folder", "", "folder to restart from");
  s.runtime<std::string>("folder", "", "output folder");
  s.runtime<double>("t", 0.0, "current time");
  s.runtime<Uint>("it", 0, "current step");
  s.runtime<double>("Lx", 0.0, "domain size, from the interpolator");
  s.runtime<double>("Ly", 0.0, "domain size, from the interpolator");
  s.runtime<double>("Lz", 0.0, "domain size, from the interpolator");
  s.choices("mode", {"analytic", "structured", "lbm", "felbm", "fenics",
                     "tet", "triangle", "trianglefreq", "xdmftriangle", "xdmftet"});
  s.choices("scheme", {"explicit", "RK4"});
  s.choices("resize", {"rescale", "doublings"});
  s.choices("resize_target", {"ds_max", "ds_init"});
  s.choices("outside", {"ignore", "reinject"});
  s.check([](const partrac::Params& p){ return p.get<int>("int_order") <= 2; },
          "int_order must be 1 or 2");
  // Pairs or points only
  s.check([](const partrac::Params& p){
            const auto key = split_string(p.get<std::string>("init_mode"), "_");
            return !key.empty() && (key[0] == "pair" || key[0] == "pairs" || key[0] == "points");
          },
          "init_mode must be pair_*, pairs_* or points_*");
  // No injection
  s.check([](const partrac::Params& p){ return !p.get<bool>("inject"); },
          "filaments does not inject");
  // init_mode needs a direction
  s.check([](const partrac::Params& p){
            return split_string(p.get<std::string>("init_mode"), "_").size() >= 2;
          },
          "init_mode is missing a direction, as in pairs_xyz or points_xy");
  // RK4 ignores Dm
  s.warn([](const partrac::Params& p){
           return p.get<std::string>("scheme") == "RK4" && p.get<double>("Dm") != 0.0;
         },
         "scheme=RK4 ignores Dm");
  // No reinjection with diffusion
  s.check([](const partrac::Params& p){
            return !(p.get<std::string>("outside") == "reinject"
                     && p.get<std::string>("scheme") == "explicit"
                     && p.get<double>("Dm") > 0.);
          },
          "outside=reinject is for edges stuck in an underresolved field, not for diffusion: it cannot be combined with scheme=explicit and Dm > 0");
  // Floor output intervals at one step; 0 is off, negative an error
  s.check([](const partrac::Params& p){
            for (const auto& key : {"checkpoint_intv", "dump_intv", "resize_intv", "stat_intv"})
              if (p.get<double>(key) < 0.) return false;
            return true;
          },
          "intervals cannot be negative");
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
