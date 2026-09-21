#ifndef __TRACERPARAMS_HPP
#define __TRACERPARAMS_HPP

#include <algorithm>
#include <string>
#include "typedefs.hpp"
#include "Params.hpp"

// Tracer app parameters
struct TracerDefaults {
  std::string outside;     // default for outside
  int int_order;           // fixed order, or 0 if required
  bool output_phi;         // what it dumped by default
  bool output_cell_type;
  bool output_J;
};

namespace tracer_params_detail {

inline void add_common(partrac::Schema& s, const TracerDefaults& d){
  s.require<std::string>("mode", "interpolator type");
  s.require<std::string>("init_mode", "initial distribution: points_<directions>, as points_xy");
  s.require<Uint>("Nrw", "number of particles");
  s.require<Uint>("Nrw_max", "max number of particles");
  if (d.int_order > 0)
    s.opt<int>("int_order", d.int_order, "integration order");
  else
    s.require<int>("int_order", "integration order");
  s.opt<std::string>("init_weight", "uniform", "what the initial points are sampled by: uniform, u, ux, uy or uz");
  s.opt<int>("sort_every", 0, "reorder particles by cell every this many steps, 0 = never");
  s.opt<bool>("output_J", d.output_J, "dump the velocity gradient at each particle");
  s.opt<bool>("output_phi", d.output_phi, "dump the phase field at each particle");
  s.opt<bool>("output_cell_type", d.output_cell_type, "dump the cell marker at each particle");
  s.opt<double>("t0", 0.0, "start time");
  s.opt<double>("x0", 0.0, "initial position");
  s.opt<double>("y0", 0.0, "initial position");
  s.opt<double>("z0", 0.0, "initial position");
  s.opt<double>("U", 1.0, "velocity scale");
  s.opt<double>("dump_intv", 100.0, "dump interval");
  s.opt<double>("stat_intv", 100.0, "statistics interval");
  s.opt<double>("checkpoint_intv", 1000.0, "checkpoint interval");
  s.opt<int>("seed", 0, "random seed");
  s.opt<bool>("random", true, "draw the seed randomly");
  s.opt<int>("num_threads", 0, "OpenMP threads, 0 = leave alone");
  s.opt<int>("dump_chunk_size", 0, "particles per dump chunk");
  s.opt<bool>("minimal_output", false, "dump less");
  s.opt<bool>("output_all_props", true, "dump all properties");
  s.opt<bool>("verbose", false, "print the parameters");
  s.opt<std::string>("tag", "", "appended to the folder name");
  s.opt<std::string>("restart_folder", "", "folder to restart from");
  s.runtime<std::string>("folder", "", "output folder");
  s.runtime<double>("t", 0.0, "current time, or path length in a march");
  s.runtime<Uint>("it", 0, "current step");
  s.runtime<double>("Lx", 0.0, "domain size, from the interpolator");
  s.runtime<double>("Ly", 0.0, "domain size, from the interpolator");
  s.runtime<double>("Lz", 0.0, "domain size, from the interpolator");
  s.runtime<Uint>("Nrw_init", 0, "particles the initializer placed");
  s.runtime<Uint>("Nrw_current", 0, "particles in the set when this was written");
  // Cloud: no edges, refinement or injection
  s.runtime<double>("ds_init", 0.0, "a cloud: no initial edges");
  s.runtime<double>("ds_max", 0.0, "a cloud: no edges");
  s.runtime<double>("ds_min", 0.0, "a cloud: no edges");
  s.runtime<bool>("inject", false, "a cloud: no injection");
  s.runtime<bool>("inject_edges", true, "a cloud: no injection");
  s.runtime<bool>("clear_initial_edges", false, "a cloud: no initial edges");
  s.runtime<bool>("cut_if_stuck", true, "a cloud: no edges");
  s.runtime<double>("curv_refine_factor", 0.0, "a cloud: no refinement");
  s.runtime<int>("filter_target", 0, "a cloud: no filtering");

  s.choices("mode", {"analytic", "structured", "lbm", "felbm", "fenics",
                     "tet", "triangle", "trianglefreq", "tetfreq", "xdmftriangle", "xdmftet"});
  s.token_choices("init_mode", "_", {"points"});
  s.check([](const partrac::Params& p){ return p.get<int>("int_order") <= 2; },
          "int_order must be 1 or 2");
  // Intervals: 0 is off, negative is an error
  s.check([](const partrac::Params& p){
            for (const auto& key : {"checkpoint_intv", "dump_intv", "stat_intv"})
              if (p.get<double>(key) < 0.) return false;
            return true;
          },
          "intervals cannot be negative");
}

// Floor output intervals at one step
inline void floor_intervals(partrac::Params& p, const double step){
  if (p.get<double>("dump_intv") > 0.)
    p.set<double>("dump_intv", std::max(p.get<double>("dump_intv"), step));
  if (p.get<double>("stat_intv") > 0.)
    p.set<double>("stat_intv", std::max(p.get<double>("stat_intv"), step));
  p.set<Uint>("Nrw_max", std::max(p.get<Uint>("Nrw_max"), p.get<Uint>("Nrw")));
}

}  // namespace tracer_params_detail

// Time stepping tracers
inline void add_tracer_params(partrac::Schema& s, const TracerDefaults& d){
  tracer_params_detail::add_common(s, d);
  s.require<double>("Dm", "molecular diffusivity");
  s.require<double>("dt", "timestep");
  s.require<double>("T", "final time");
  s.opt<std::string>("scheme", "RK4", "ODE integration scheme; explicit is the diffusive step");
  s.opt<std::string>("outside", d.outside, "a particle that cannot take its step: ignore (it stays), reinject at a random offset, or mark (c = 2)");
  s.choices("scheme", {"explicit", "RK4"});
  s.choices("outside", {"ignore", "reinject", "mark"});
  // RK4 ignores Dm
  s.warn([](const partrac::Params& p){
           return p.get<std::string>("scheme") == "RK4" && p.get<double>("Dm") != 0.0;
         },
         "scheme=RK4 ignores Dm");
  s.finalize([](partrac::Params& p){
    tracer_params_detail::floor_intervals(p, p.get<double>("dt"));
  });
}

// Marching tracers
inline void add_march_tracer_params(partrac::Schema& s, const TracerDefaults& d){
  tracer_params_detail::add_common(s, d);
  s.require<double>("dxn", "path length of a step");
  s.opt<double>("Ln", 0.0, "path length the march ends at");
  s.opt<double>("xn0", 0.0, "path length the march starts from");
  s.opt<double>("T", 1e9, "integration time after which a particle stops");
  s.opt<double>("u_eps", 1e-7, "a particle this slow or slower cannot move");
  s.opt<double>("dx_max", 1e9, "longest step accepted");
  s.opt<double>("Dm", 0.0, "diffusivity, enters the folder name only");
  s.opt<double>("dt", 1.0, "timestep, enters the folder name only");
  s.opt<std::string>("outside", d.outside, "a particle that cannot take its step (outside, too slow, done, or a step too long): ignore (it stays), mark (c = 2), or remove");
  s.choices("outside", {"ignore", "mark", "remove"});
  s.finalize([](partrac::Params& p){
    tracer_params_detail::floor_intervals(p, p.get<double>("dxn"));
  });
}

#endif
