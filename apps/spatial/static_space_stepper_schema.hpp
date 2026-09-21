#ifndef __STATIC_SPACE_STEPPER_SCHEMA_HPP
#define __STATIC_SPACE_STEPPER_SCHEMA_HPP

#include "typedefs.hpp"
#include "Params.hpp"
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
  s.opt<double>("U", 1.0, "velocity scale");
  s.opt<double>("Ln", 0.0, "path length the march ends at");
  s.opt<double>("xn0", 0.0, "path length the march starts from");
  s.opt<double>("u_eps", 1e-7, "a node this slow or slower cannot move");
  s.opt<std::string>("outside", "remove", "a node that cannot take its step (outside, too slow, done, or a step longer than dx_max): ignore (it stays), mark (c = 2), or remove");
  s.opt<int>("sort_every", 5, "reorder particles by cell every this many steps, 0 = never");
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
  s.opt<bool>("inject_edges", true, "inject edges too");
  s.opt<bool>("local_dt", false, "use a local timestep");
  s.opt<bool>("cut_if_stuck", true, "cut edges that get stuck");
  s.opt<bool>("output_all_props", true, "dump all properties");
  s.opt<bool>("minimal_output", false, "dump less");
  s.opt<std::string>("tag", "", "appended to the folder name");
  s.opt<std::string>("restart_folder", "", "folder to restart from");
  s.runtime<std::string>("folder", "", "output folder");
  s.runtime<double>("t", 0.0, "path length marched");
  s.runtime<Uint>("it", 0, "current step");
  s.runtime<double>("Lx", 0.0, "domain size, from the interpolator");
  s.runtime<double>("Ly", 0.0, "domain size, from the interpolator");
  s.runtime<double>("Lz", 0.0, "domain size, from the interpolator");
  s.choices("mode", {"analytic", "structured", "lbm", "felbm", "fenics",
                     "tet", "triangle", "trianglefreq", "tetfreq", "xdmftriangle", "xdmftet"});
  s.choices("outside", {"ignore", "mark", "remove"});
  s.check([](const partrac::Params& p){ return p.get<int>("int_order") <= 2; },
          "int_order must be 1 or 2");
  // Floor output intervals at one step; 0 is off, negative an error
  s.check([](const partrac::Params& p){
            for (const auto& key : {"checkpoint_intv", "coarsen_intv", "dump_intv", "refine_intv", "stat_intv"})
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
