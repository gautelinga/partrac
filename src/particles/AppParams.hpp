#ifndef __APPPARAMS_HPP
#define __APPPARAMS_HPP

// Parameter blocks the apps' schemas share

#include <algorithm>
#include <string>
#include <vector>
#include "typedefs.hpp"
#include "Params.hpp"

// Every app: field scale, seed, threads, folder naming and restart
inline void add_app_params(partrac::Schema& s){
  s.opt<double>("U", 1.0, "velocity scale");
  s.opt<int>("seed", 0, "random seed");
  s.opt<bool>("random", true, "draw the seed randomly");
  s.opt<int>("num_threads", 0, "OpenMP threads, 0 = leave alone");
  s.opt<std::string>("tag", "", "appended to the folder name");
  s.opt<std::string>("restart_folder", "", "folder to restart from");
  s.runtime<std::string>("folder", "", "output folder");
}

// Every particle run: the app's keys, output, and the state a checkpoint keeps;
// t is what the run counts, a time or a path length
inline void add_run_params(partrac::Schema& s, const std::string& t_doc){
  add_app_params(s);
  s.opt<int>("dump_chunk_size", 0, "particles per dump chunk");
  s.opt<bool>("minimal_output", false, "dump less");
  s.opt<bool>("verbose", false, "print the parameters");
  s.runtime<double>("t", 0.0, t_doc);
  s.runtime<Uint>("it", 0, "current step");
  s.runtime<double>("Lx", 0.0, "domain size, from the interpolator");
  s.runtime<double>("Ly", 0.0, "domain size, from the interpolator");
  s.runtime<double>("Lz", 0.0, "domain size, from the interpolator");
}

// dump_intv, stat_intv and checkpoint_intv. These and the keys in more: 0 is off,
// negative an error; dump and statistics floored at one step of step_key, and
// Nrw_max raised to Nrw
inline void add_intervals(partrac::Schema& s, const std::string& step_key,
                          const std::vector<std::string>& more = {}){
  s.opt<double>("dump_intv", 100.0, "dump interval");
  s.opt<double>("stat_intv", 100.0, "statistics interval");
  s.opt<double>("checkpoint_intv", 1000.0, "checkpoint interval");
  std::vector<std::string> keys = {"checkpoint_intv", "dump_intv", "stat_intv"};
  keys.insert(keys.end(), more.begin(), more.end());
  s.check([keys](const partrac::Params& p){
            for (const auto& key : keys)
              if (p.get<double>(key) < 0.) return false;
            return true;
          },
          "intervals cannot be negative");
  s.finalize([step_key](partrac::Params& p){
    const double step = p.get<double>(step_key);
    if (p.get<double>("dump_intv") > 0.)
      p.set<double>("dump_intv", std::max(p.get<double>("dump_intv"), step));
    if (p.get<double>("stat_intv") > 0.)
      p.set<double>("stat_intv", std::max(p.get<double>("stat_intv"), step));
    p.set<Uint>("Nrw_max", std::max(p.get<Uint>("Nrw_max"), p.get<Uint>("Nrw")));
  });
}

#endif
