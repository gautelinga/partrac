#ifndef __PARTRAC_SCHEMA_HPP
#define __PARTRAC_SCHEMA_HPP

#include "typedefs.hpp"
#include "Params.hpp"
#include "Initializer.hpp"

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
  s.opt<double>("U", 1.0, "velocity scale");
  s.opt<double>("T_inject", 1e10, "time to stop injecting at");
  s.opt<double>("tau_max", 0.0, "max tau");
  s.opt<double>("t_frozen", 0.0, "time to freeze the fields at");
  s.opt<double>("dump_intv", 100.0, "dump interval");
  s.opt<double>("stat_intv", 100.0, "statistics interval");
  s.opt<double>("checkpoint_intv", 1000.0, "checkpoint interval");
  s.opt<double>("refine_intv", 100.0, "refinement interval");
  s.opt<double>("coarsen_intv", 1000.0, "coarsening interval");
  s.opt<double>("curv_refine_factor", 0.0, "curvature refinement factor");
  s.opt<int>("seed", 0, "random seed");
  s.opt<int>("num_threads", 0, "OpenMP threads, 0 = leave alone");
  s.opt<int>("dump_chunk_size", 0, "particles per dump chunk");
  s.opt<int>("filter_target", 0, "filter target");
  s.opt<bool>("verbose", false, "print the parameters");
  s.opt<bool>("random", true, "draw the seed randomly");
  s.opt<bool>("refine", false, "refine the mesh");
  s.opt<bool>("coarsen", false, "coarsen the mesh");
  s.opt<bool>("filter", false, "filter the mesh");
  s.opt<bool>("inject_edges", true, "inject edges too");
  s.opt<bool>("frozen_fields", false, "freeze the velocity field");
  s.opt<bool>("cut_if_stuck", true, "cut edges that get stuck");
  s.opt<bool>("integrate_tau", false, "integrate the eigentime");
  s.opt<bool>("output_all_props", true, "dump all properties");
  s.opt<bool>("minimal_output", false, "dump less");
  s.opt<std::string>("scheme", "explicit", "ODE integration scheme");
  s.opt<std::string>("exit_plane", "none", "plane to remove particles beyond");
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
  s.choices("outside", {"ignore", "reinject", "mark"});
  s.choices("exit_plane", {"none", "x", "y", "z"});
  s.check([](const partrac::Params& p){ return p.get<int>("int_order") <= 2; },
          "int_order must be 1 or 2");
  s.check([](const partrac::Params& p){
            return !(p.get<bool>("inject") && p.get<bool>("filter"));
          },
          "cannot inject and filter at the same time");
  // RK4 ignores Dm
  s.warn([](const partrac::Params& p){
           return p.get<std::string>("scheme") == "RK4" && p.get<double>("Dm") != 0.0;
         },
         "scheme=RK4 ignores Dm");
  // Floor output intervals at one step; 0 is off, negative an error
  s.check([](const partrac::Params& p){
            for (const auto& key : {"checkpoint_intv", "coarsen_intv", "dump_intv", "filter_intv", "inject_intv", "refine_intv", "stat_intv", "tau_intv"})
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
