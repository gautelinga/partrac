#ifndef __RUNLOOP_HPP
#define __RUNLOOP_HPP

// Shared run loop: setup and per-step cadence

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <vector>
#include <omp.h>
#include "H5Cpp.h"

#include "typedefs.hpp"
#include "Params.hpp"
#include "strings.hpp"
#include "intervals.hpp"
#include "perf_window.hpp"
#include "io.hpp"
#include "rng.hpp"
#include "run_folders.hpp"
#include "interpol_factory.hpp"
#include "Initializer.hpp"
#include "stepping.hpp"
#include "ParticleSet.hpp"
#include "Topology.hpp"
#include "Integrator.hpp"
#include "stats_columns.hpp"

// Run setup
struct Run {
  partrac::Params& prm;
  std::shared_ptr<Interpol> intp;
  RunFolders out;
  std::vector<std::mt19937> gens;
  double t0;
  double T;
  bool frozen_fields;
  // March: t0 and T are path lengths, fields evaluated at t_fields
  bool marches;
  double t_fields;
};

// Threads, interpolator, folders, generators and time window
inline Run start_run(partrac::Params& prm, const std::string& name,
                     const unsigned folder_opts = DefaultLayout, const bool marches = false){
  if (prm.get<int>("num_threads") > 0){
    omp_set_dynamic(0);
    omp_set_num_threads(prm.get<int>("num_threads"));
  }

  std::cout << "Setting interpolator..." << std::endl;
  std::shared_ptr<Interpol> intp;
  set_interpolate_mode(intp, prm.get<std::string>("mode"), prm.input_file());
  intp->set_U0(prm.get<double>("U"));
  intp->set_int_order(prm.get<int>("int_order"));

  // Parallel generators; before the folders, which name the seed
  std::vector<std::mt19937> gens = make_generators(prm);

  std::cout << "Creating folders..." << std::endl;
  RunFolders out = make_run_folders(intp->get_folder(), name, prm,
                                    prm.check_only() ? folder_opts | DryRun : folder_opts);

  if (prm.get<bool>("verbose"))
    prm.print();

  // Domain size for the dumped parameters: runtime entries, refused as input
  prm.set<double>("Lx", intp->get_Lx());
  prm.set<double>("Ly", intp->get_Ly());
  prm.set<double>("Lz", intp->get_Lz());

  const double t0 = std::max(intp->get_t_min(), prm.get<double>("t0"));
  prm.set<double>("t0", t0);
  if (marches){
    // Fields frozen at t0; T is the integration time
    intp->update(t0);
    return Run{prm, intp, out, std::move(gens), prm.get<double>("xn0"), prm.get<double>("Ln"),
               true, true, t0};
  }
  const bool frozen_fields = prm.has("frozen_fields") && prm.get<bool>("frozen_fields");
  double T = std::min(intp->get_t_max(), prm.get<double>("T"));
  if (frozen_fields)
    T = prm.get<double>("T");
  prm.set<double>("T", T);

  if (frozen_fields)
    intp->freeze(prm.get<double>("t_frozen"));
  else
    intp->update(t0);

  // Walls for the explicit diffusive step
  const bool diffuses = prm.has("Dm") && prm.get<double>("Dm") > 0.;
  if (diffuses && (!prm.has("scheme") || prm.get<std::string>("scheme") == "explicit")){
    intp->enable_reflection();
    // Warn when the noise step exceeds a cell
    const double sigma = sqrt(2*prm.get<double>("Dm")*prm.get<double>("dt"));
    const double h = intp->hmin();
    if (intp->can_reflect && h > 0. && sigma > h)
      std::cout << "Note: diffusive step " << sigma << " exceeds the smallest cell " << h << std::endl;
  }

  return Run{prm, intp, out, std::move(gens), t0, T, frozen_fields, false, t0};
}

// Load checkpoint or initial state
inline bool load_or_initialize(Run& run, Topology& mesh){
  const bool restarting = run.prm.get<std::string>("restart_folder") != "";
  if (restarting){
    mesh.load_checkpoint(run.prm.get<std::string>("restart_folder") + "/Checkpoints", run.prm);
  }
  else {
    mesh.load_initial_state(set_initial_state(run.intp, run.prm, run.gens[0]), run.prm);
  }
  mesh.compute_maps();
  return restarting;
}

// --check: loaded and initialized, nothing run or written
inline bool check_only(const Run& run, const ParticleSet& ps, Topology& mesh){
  if (!run.prm.check_only())
    return false;
  std::cout << "Check OK: " << ps.N() << " particles, dim = " << mesh.dim() << std::endl;
  return true;
}

// Particles that could not step: ignore, reinject, mark or remove
inline void handle_outside(Run& run, Topology& mesh, ParticleSet& ps, const std::string& outside,
                           const std::vector<Uint>& nodes, const double t, const bool verbose){
  if (nodes.empty())
    return;
  if (outside == "reinject"){
    // Along init_mode's directions; positions from a file have none: all three
    const std::string init_mode = run.prm.get<std::string>("init_mode");
    const auto key = split_string(init_mode, "_");
    const std::string dirs = init_mode_is_file(init_mode) || key.size() < 2 ? "xyz" : key[1];
    const Vector3d Dx_max = 0.5*(run.intp->get_x_max() - run.intp->get_x_min());
    std::uniform_real_distribution<> ux(-Dx_max[0], Dx_max[0]), uy(-Dx_max[1], Dx_max[1]), uz(-Dx_max[2], Dx_max[2]);
    const bool rx = contains(dirs, "x"), ry = contains(dirs, "y"), rz = contains(dirs, "z");
    const Uint max_draws = 1000000;
    for (const Uint i : nodes){
      Vector3d Dx = {0., 0., 0.};
      Uint draws = 0;
      do {
        if (++draws > max_draws)
          partrac::fail("outside=reinject: no position inside the domain along ", dirs, " in ", max_draws, " draws");
        if (rx) Dx[0] = ux(run.gens[0]);
        if (ry) Dx[1] = uy(run.gens[0]);
        if (rz) Dx[2] = uz(run.gens[0]);
      } while (!run.intp->locate(ps.x(i) + Dx));
      ps.move(i, ps.x(i) + Dx);
    }
  }
  if (outside == "mark")
    for (const Uint i : nodes) ps.set_c(i, 2.0);
  if (verbose){
    Vector3d x_stuck = {0., 0., 0.};
    for (const Uint i : nodes)
      x_stuck += ps.x(i);
    x_stuck /= nodes.size();
    std::cout << nodes.size() << " nodes could not move at t = " << t
              << ", centred on (" << x_stuck[0] << ", " << x_stuck[1] << ", "
              << x_stuck[2] << ")" << std::endl;
  }
  // Remove last: slots shift
  if (outside == "remove"){
    std::vector<bool> node_isactive(ps.N(), true);
    for (const Uint i : nodes)
      node_isactive[i] = false;
    mesh.remove_nodes_safe(node_isactive);
  }
}

// A step over several cells, once, from the fields just updated: it cuts across
// streamlines and the walls' layers
inline void note_step_size(Run& run, const ParticleSet& ps, const double t_fields, const double dt){
  double worst = 0.;
  Uint over = 0, counted = 0;
  #pragma omp parallel for reduction(max:worst) reduction(+:over,counted)
  for (Uint i = 0; i < ps.N(); ++i){
    CellPos pos;
    pos.id = ps.get_cell_id(i);
    if (!run.intp->locate(ps.x(i), t_fields, pos))
      continue;
    const double h = run.intp->cell_size(pos.id);
    if (!(h > 0.))
      continue;
    // A march steps in path length
    const double cells = (run.marches ? 1. : ps.u(i).norm())*dt/h;
    worst = std::max(worst, cells);
    ++counted;
    if (cells > 1.) ++over;
  }
  if (worst > 1.)
    std::cout << "Note: a step crosses up to " << worst << " cells, more than one for "
              << over << " of " << counted << " particles (dt = " << dt << ")" << std::endl;
}

// App hooks
struct RunHooks {
  // After field update: injection, remeshing, resizing, removal
  std::function<void(int it, double t)> reshape = [](int, double){};
  // After each piece of a step: the particles that could not step
  std::function<void(const std::vector<Uint>& nodes, double t)> outside =
      [](const std::vector<Uint>&, double){};
  // After step, once
  std::function<void(int it, double t)> after_step = [](int, double){};
  // Statistics columns; unset: mesh statistics
  std::function<std::vector<StatsColumn>(double t, Integrator& counters)> statistics;
  // After statistics row
  std::function<void(int it, double t)> after_statistics = [](int, double){};
  // Continue?
  std::function<bool()> keep_going = []{ return true; };
};

// Simulation loop
template<typename Stepper>
void run_loop(Run& run, ParticleSet& ps, Topology& mesh, Stepper& stepper,
              std::map<std::string, bool>& output_fields, const double dt,
              const RunHooks& hooks = RunHooks()){
  partrac::Params& prm = run.prm;
  const bool restarting = prm.get<std::string>("restart_folder") != "";
  // Resume step count
  int it = restarting ? static_cast<int>(prm.get<Uint>("it")) : 0;
  // Time from the step count since the run's start, else since the restart
  double t_start = run.t0;
  int it_start = 0;
  if (restarting && t_start + it*dt != prm.get<double>("t")){
    t_start = prm.get<double>("t");
    it_start = it;
  }
  const auto time_at = [&](const int k){ return t_start + (k - it_start)*dt; };
  double t = time_at(it);

  const std::string& newfolder = run.out.run;
  prm.dump(newfolder, t);

  // Loop constants
  const double dump_intv = prm.get<double>("dump_intv");
  const double stat_intv = prm.get<double>("stat_intv");
  const double checkpoint_intv = prm.get<double>("checkpoint_intv");
  const double chunk_intv = dump_intv*prm.get<int>("dump_chunk_size");
  const double ds_max = prm.get<double>("ds_max");
  const int sort_every = prm.has("sort_every") ? prm.get<int>("sort_every") : 0;
  // S is computed in the refresh
  const bool refresh_S = ps.carries() == TransportElement::Vector;

  std::string h5fname = newfolder + "/data_from_t" + std::to_string(t) + ".h5";
  // No dump file when dumping is off
  H5::H5File h5f;
  if (dump_intv > 0.){
    try {
      { H5::H5File create(h5fname.c_str(), H5F_ACC_TRUNC); }
      h5f.openFile(h5fname.c_str(), H5F_ACC_RDWR);
    } catch (const H5::Exception&){
      partrac::fail("cannot create the dump file ", h5fname);
    }
  }

  std::ofstream statfile;
  if (stat_intv > 0.){
    statfile.open(newfolder + "/tdata_from_t" + std::to_string(t) + ".dat");
    write_stats_header(statfile, hooks.statistics ? hooks.statistics(0., stepper.counters())
                                                  : mesh.stats_header_columns(ds_max));
  }

  bool step_noted = false;

  // Steps cut at stamps; a stamp within snap of a piece's end makes no piece
  const bool cuts = !run.frozen_fields;
  const double snap = 1e-9*dt;
  run.intp->set_stamp_snap(snap);
  const double t_last = run.intp->get_t_max();
  const auto next_stamp = [&](const double ta){
    const double s = run.intp->next_stamp_after(ta);
    return s < t_last ? s : std::numeric_limits<double>::infinity();
  };

  // Counted from the second step: the first refreshes and checkpoints
  PerfWindow counters;
  const int it_counted = it + 1;
  const double wall_0 = omp_get_wtime();   // wall clock, not the threads' summed CPU time
  while (t < run.T + dt/2 && hooks.keep_going()){
    if (it == it_counted)
      counters.enable();
    // Sort by cell
    if (sort_every > 0 && it % sort_every == 0 && it > 0)
      mesh.sort_by_cell();

    if (!run.frozen_fields)
      run.intp->update(t);

    hooks.reshape(it, t);

    // Update fields
    if (at_interval(it, dump_intv, dt) || at_interval(it, stat_intv, dt)
        || (refresh_S && at_interval(it, checkpoint_intv, dt))){
      const double t_fields = run.marches ? run.t_fields : t;
      update_fields(ps, *run.intp, t_fields, output_fields);
      if (!step_noted){
        note_step_size(run, ps, t_fields, dt);
        step_noted = true;
      }
    }

    // Statistics
    if (at_interval(it, stat_intv, dt)){
      std::cout << "Time = " << t << std::endl;
      if (hooks.statistics)
        write_stats_row(statfile, hooks.statistics(t, stepper.counters()));
      else
        mesh.write_statistics(statfile, t, ds_max, stepper.counters());
      hooks.after_statistics(it, t);
    }

    // Checkpoint
    if (at_interval(it, checkpoint_intv, dt)){
      prm.set<Uint>("it", it);
      mesh.write_checkpoint(run.out.checkpoints, t, prm);
    }

    // Dump detailed data
    if (at_interval(it, dump_intv, dt)){
      std::string groupname = std::to_string(t);
      try {
        // Clear file if it exists, otherwise create
        if (at_interval(it, chunk_intv, dt) && it > 0){
          h5fname = newfolder + "/data_from_t" + std::to_string(t) + ".h5";
          h5f.openFile(h5fname.c_str(), H5F_ACC_TRUNC);
        }
        else {
          h5f.openFile(h5fname.c_str(), H5F_ACC_RDWR);
        }
        h5f.createGroup(groupname + "/");
        mesh.dump_hdf5(h5f, groupname, output_fields);
        h5f.close();
      } catch (const H5::Exception&){
        partrac::fail("cannot write the dump at t = ", t, " to ", h5fname);
      }
    }

    // Pieces between stamps
    const double t_end = time_at(it + 1);
    double ta = t;
    while (true){
      double tb = t_end;
      if (cuts){
        // Stamps just after the piece's start: taken as its start
        double s = next_stamp(ta), t_upd = ta;
        while (s <= ta + snap){
          t_upd = s;
          s = next_stamp(s);
        }
        if (t_upd != t)
          run.intp->update(t_upd);
        if (s < t_end - snap)
          tb = s;
      }
      // A whole step is dt
      const double h = (ta == t && tb == t_end) ? dt : tb - ta;
      const auto outside_nodes = stepper.step(*run.intp, ps, ta, h);
      hooks.outside(outside_nodes, ta);
      if (tb == t_end)
        break;
      ta = tb;
    }

    hooks.after_step(it, t);

    it += 1;
    t = time_at(it);
  }
  const double duration = omp_get_wtime() - wall_0;
  counters.disable();
  std::cout << "Total simulation time: " << duration << " seconds" << std::endl;

  if (refresh_S){
    const double t_fields = run.marches ? run.t_fields : t;
    update_fields(ps, *run.intp, t_fields, output_fields);
  }
  prm.set<Uint>("it", it);
  mesh.write_checkpoint(run.out.checkpoints, t, prm);

  statfile.close();
}

#endif
