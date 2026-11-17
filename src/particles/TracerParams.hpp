#ifndef __TRACERPARAMS_HPP
#define __TRACERPARAMS_HPP

#include <string>
#include "typedefs.hpp"
#include "Params.hpp"
#include "interpol_factory.hpp"
#include "AppParams.hpp"
#include "Initializer.hpp"

// Tracer app parameters
struct TracerDefaults {
  std::string outside;     // default for outside
  int int_order;           // fixed order, or 0 if required
  bool output_phi;         // what it dumped by default
  bool output_cell_type;
  bool output_J;
};

// The time scheme, def by default; reads Dm
inline void add_scheme(partrac::Schema& s, const std::string& def){
  s.opt<std::string>("scheme", def, "ODE integration scheme; explicit is the diffusive step, RK4cells RK4 cut at the cells' facets");
  s.choices("scheme", {"explicit", "RK4", "RK4cells"});
  s.check([](const partrac::Params& p){
            return p.get<std::string>("scheme") != "RK4cells" || p.get<double>("Dm") == 0.0;
          },
          "scheme=RK4cells has no diffusion: Dm must be 0");
  // RK4 ignores Dm
  s.warn([](const partrac::Params& p){
           return p.get<std::string>("scheme") == "RK4" && p.get<double>("Dm") != 0.0;
         },
         "scheme=RK4 ignores Dm");
}

// First crossings of the levels of n.(x - origin) (crossings.hpp); with
// needs_n, cross_intv needs cross_n (partrac's exit_plane=x|y|z names it too)
inline void add_crossing_params(partrac::Schema& s, const bool needs_n = true){
  s.opt<std::string>("cross_n", "", "plane normal: x, y, z or nx,ny,nz");
  s.opt<std::string>("cross_intv", "0", "record first crossings of the levels of n.(x - cross_origin): every "
                     "X, or a schedule of X, X:U, logN:S:U, logN:U segments, as distances; 0: none");
  s.opt<std::string>("cross_origin", "0,0,0", "where the levels are measured from: x,y,z, or start, each "
                     "point's position when first seen");
  s.opt<std::string>("cross_file", "crossings.bin",
                     "crossing records, in the run folder: a text line naming the float64 columns, then "
                     "the records");
  s.check([](const partrac::Params& p){
            Schedule sc;
            return schedule_message("cross_intv", parse_schedule(p.get<std::string>("cross_intv"), sc));
          }, "cross_intv's schedule");
  s.check([](const partrac::Params& p){
            bool start;
            Vector3d o;
            return try_parse_origin(p.get<std::string>("cross_origin"), start, o);
          },
          "cross_origin must be start or x,y,z");
  s.check([](const partrac::Params& p){
            Vector3d n;
            const std::string v = p.get<std::string>("cross_n");
            return v.empty() || (v != "none" && try_parse_normal(v, n));
          },
          "cross_n must be x, y, z or nx,ny,nz, not zero");
  if (needs_n)
    s.check([](const partrac::Params& p){
              return !crossings_on(p) || !p.get<std::string>("cross_n").empty();
            },
            "cross_intv needs cross_n");
}

namespace tracer_params_detail {

inline void add_common(partrac::Schema& s, const TracerDefaults& d, const std::string& step_key){
  s.require<std::string>("mode", "interpolator type");
  s.require<std::string>("init_mode", "initial distribution: points_<directions>, as points_xy");
  s.require<Uint>("Nrw", "number of particles");
  s.require<Uint>("Nrw_max", "max number of particles");
  add_int_order(s, "integration order", d.int_order);
  s.opt<std::string>("init_weight", "uniform", "what the initial points are sampled by: uniform, u, ux, uy or uz");
  s.choices("init_weight", {"uniform", "u", "ux", "uy", "uz"});
  s.opt<int>("sort_every", 0, "reorder particles by cell every this many steps, 0 = never");
  s.opt<bool>("output_J", d.output_J, "dump the velocity gradient at each particle");
  s.opt<bool>("output_S", false, "dump the stretching rates along the frame (tracertensors)");
  s.opt<bool>("output_phi", d.output_phi, "dump the phase field at each particle");
  s.opt<bool>("output_cell_type", d.output_cell_type, "dump the cell marker at each particle");
  s.opt<double>("t0", 0.0, "start time");
  s.opt<double>("x0", 0.0, "initial position");
  s.opt<double>("y0", 0.0, "initial position");
  s.opt<double>("z0", 0.0, "initial position");
  add_run_params(s, "current time, or path length in a march");
  add_intervals(s, step_key);
  s.opt<bool>("output_all_props", true, "dump all properties");
  s.runtime<Uint>("Nrw_init", 0, "particles the initializer placed");
  s.runtime<Uint>("Nrw_current", 0, "particles in the set when this was written");
  // Cloud: no edges, refinement or injection
  s.runtime<double>("ds_init", 0.0, "a cloud: no initial edges");
  s.runtime<double>("ds_max", 0.0, "a cloud: no edges");
  s.runtime<double>("ds_min", 0.0, "a cloud: no edges");
  s.runtime<bool>("inject", false, "a cloud: no injection");
  s.runtime<bool>("clear_initial_edges", false, "a cloud: no initial edges");
  // Written by older runs; a cloud's Topology takes none of them
  s.retired({"inject_edges", "cut_if_stuck", "curv_refine_factor", "filter_target"});

  s.choices("mode", interpol_modes());
  s.token_choices("init_mode", "_", {"points"});
}

}  // namespace tracer_params_detail

// Time stepping tracers
inline void add_tracer_params(partrac::Schema& s, const TracerDefaults& d){
  tracer_params_detail::add_common(s, d, "dt");
  s.require<double>("Dm", "molecular diffusivity");
  s.require<double>("dt", "timestep");
  s.require<double>("T", "final time");
  s.opt<bool>("frozen_fields", false, "freeze the velocity field");
  s.optional<double>("t_frozen", "time to freeze the fields at, clamped to the fields' times; default t0");
  add_scheme(s, "RK4");
  add_outside(s, d.outside, "a particle that cannot take its step", false);
  add_crossing_params(s);
}

// Marching tracers
inline void add_march_tracer_params(partrac::Schema& s, const TracerDefaults& d){
  tracer_params_detail::add_common(s, d, "dxn");
  s.require<double>("dxn", "path length of a step");
  s.opt<double>("Ln", 0.0, "path length the march ends at");
  s.opt<double>("xn0", 0.0, "path length the march starts from");
  s.opt<double>("T", 1e9, "integration time after which a particle stops");
  s.opt<double>("u_eps", 1e-7, "a particle this slow or slower cannot move");
  s.opt<double>("dx_max", 1e9, "longest step accepted");
  s.opt<double>("Dm", 0.0, "diffusivity, enters the folder name only");
  s.opt<double>("dt", 1.0, "timestep, enters the folder name only");
  add_outside(s, d.outside, "a particle that cannot take its step (outside, too slow, done, or a step too long)",
              true);
}

#endif
