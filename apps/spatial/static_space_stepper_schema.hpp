#ifndef __STATIC_SPACE_STEPPER_SCHEMA_HPP
#define __STATIC_SPACE_STEPPER_SCHEMA_HPP

#include "typedefs.hpp"
#include "Params.hpp"
#include "interpol_factory.hpp"
#include "AppParams.hpp"
#include "Initializer.hpp"

// Parameters accepted by this app
inline partrac::Schema spatial_schema(){
  partrac::Schema s("static_space_stepper");
  add_initializer_params(s);
  s.require<std::string>("mode", "interpolator type");
  s.require<double>("T", "integration time after which a node stops");
  s.require<int>("int_order", "interpolation order");
  s.require<double>("dx_max", "max step length");
  s.require<double>("dxn", "normal step length");
  s.opt<double>("Dm", 0.0, "diffusivity, enters the folder name only");
  s.opt<double>("dt", 1.0, "timestep, enters the folder name only");
  add_run_params(s, "path length marched");
  s.opt<double>("Ln", 0.0, "path length the march ends at");
  s.opt<double>("xn0", 0.0, "path length the march starts from");
  s.opt<double>("u_eps", 1e-7, "a node this slow or slower cannot move");
  s.opt<std::string>("outside", "remove", "a node that cannot take its step (outside, too slow, done, or a step longer than dx_max): ignore (it stays), mark (c = 2), or remove");
  s.opt<int>("sort_every", 5, "reorder particles by cell every this many steps, 0 = never");
  add_intervals(s, "dxn", {"coarsen_intv", "refine_intv"});
  s.opt<double>("refine_intv", 100.0, "refinement interval");
  s.opt<double>("coarsen_intv", 1000.0, "coarsening interval");
  s.opt<double>("curv_refine_factor", 0.0, "curvature refinement factor");
  s.opt<int>("filter_target", 0, "filter target");
  s.opt<bool>("refine", false, "refine the mesh");
  s.opt<bool>("coarsen", false, "coarsen the mesh");
  s.opt<bool>("inject_edges", true, "inject edges too");
  s.opt<bool>("local_dt", false, "use a local timestep");
  s.opt<bool>("cut_if_stuck", true, "cut edges that get stuck");
  s.opt<bool>("output_all_props", true, "dump all properties");
  s.choices("mode", interpol_modes());
  s.choices("outside", {"ignore", "mark", "remove"});
  s.check([](const partrac::Params& p){ return p.get<int>("int_order") <= 2; },
          "int_order must be 1 or 2");
  return s;
}

#endif
