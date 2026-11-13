#ifndef __APPPARAMS_HPP
#define __APPPARAMS_HPP

// Parameter blocks the apps' schemas share

#include <algorithm>
#include <string>
#include <vector>
#include "typedefs.hpp"
#include "Params.hpp"
#include "mesh.hpp"
#include "schedule.hpp"

// Every app: field scale, seed, threads and folder naming
inline void add_app_params(partrac::Schema& s){
  s.opt<double>("U", 1.0, "velocity scale");
  s.opt<int>("seed", 0, "random seed");
  s.opt<bool>("random", true, "draw the seed randomly");
  s.opt<int>("num_threads", 0, "OpenMP threads, 0 = leave alone");
  s.opt<std::string>("tag", "", "appended to the folder name");
  s.runtime<std::string>("folder", "", "output folder");
}

// int_order, 1 or 2; required unless def is
inline void add_int_order(partrac::Schema& s, const std::string& doc, const int def = 0){
  if (def > 0)
    s.opt<int>("int_order", def, doc);
  else
    s.require<int>("int_order", doc);
  s.check([](const partrac::Params& p){ const int k = p.get<int>("int_order"); return k == 1 || k == 2; },
          "int_order must be 1 or 2");
}

// Every particle run: app keys, restart, output, checkpoint state
inline void add_run_params(partrac::Schema& s, const std::string& t_doc){
  add_app_params(s);
  s.opt<std::string>("restart_folder", "", "folder to restart from");
  s.opt<int>("dump_chunk_size", 0, "dumps per dump file, 0 = one file");
  s.opt<bool>("minimal_output", false, "dump less");
  s.opt<bool>("verbose", false, "print the parameters");
  s.runtime<double>("t", 0.0, t_doc);
  s.runtime<Uint>("it", 0, "current step");
  s.runtime<double>("Lx", 0.0, "domain size, from the interpolator");
  s.runtime<double>("Ly", 0.0, "domain size, from the interpolator");
  s.runtime<double>("Lz", 0.0, "domain size, from the interpolator");
}

// The curvature limit on a strip's edge lengths (CurvLimit); no limit with factor 0
inline void add_curv_limit(partrac::Schema& s){
  s.opt<double>("curv_refine_factor", CurvLimit{}.factor, "curvature refinement factor, strips only");
  s.opt<double>("curv_refine_cap", CurvLimit{}.cap, "largest factor the curvature shortens an edge by");
  s.check([](const partrac::Params& p){ return p.get<double>("curv_refine_factor") >= 0.; },
          "curv_refine_factor cannot be negative");
  s.check([](const partrac::Params& p){ return p.get<double>("curv_refine_cap") >= 1.; },
          "curv_refine_cap must be at least 1: below, it lengthens edges");
}

// dump_intv, stat_intv and checkpoint_intv: each a number or a schedule
// (schedule.hpp). A number is written back as a number, and for dump_intv and
// stat_intv floored at one step of step_key. The keys in more are numbers. 0 is
// off, a negative interval an error. Nrw_max is raised to Nrw.
inline void add_intervals(partrac::Schema& s, const std::string& step_key,
                          const std::vector<std::string>& more = {}){
  const std::string forms = ", or a schedule of X, X:U, logN:S:U, logN:U segments";
  s.opt<std::string>("dump_intv", "100", "dump interval" + forms);
  s.opt<std::string>("stat_intv", "100", "statistics interval" + forms);
  s.opt<std::string>("checkpoint_intv", "1000", "checkpoint interval" + forms);
  const std::vector<std::string> scheduled = {"checkpoint_intv", "dump_intv", "stat_intv"};
  if (!more.empty())
    s.check([more](const partrac::Params& p){
              for (const auto& key : more)
                if (p.get<double>(key) < 0.) return false;
              return true;
            },
            "intervals cannot be negative");
  // Which segment breaks which rule
  for (const auto& key : scheduled)
    s.check([key](const partrac::Params& p){
      Schedule sc;
      return schedule_message(key, parse_schedule(p.get<std::string>(key), sc));
    });
  s.finalize([step_key, scheduled](partrac::Params& p){
    const double step = p.get<double>(step_key);
    // A number as a number; dump and statistics floored at a step
    for (const auto& key : scheduled){
      Schedule sc;
      if (parse_schedule(p.get<std::string>(key), sc) || !sc.plain()) continue;
      double intv = sc.segs[0].X;
      if (key != "checkpoint_intv" && intv > 0.)
        intv = std::max(intv, step);
      p.set<std::string>(key, partrac::value_to_string(partrac::Value(intv)));
    }
    p.set<Uint>("Nrw_max", std::max(p.get<Uint>("Nrw_max"), p.get<Uint>("Nrw")));
  });
}

#endif
