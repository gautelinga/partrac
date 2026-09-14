#ifndef __RUNLOOP_HPP
#define __RUNLOOP_HPP

// Shared run loop: setup and per-step cadence

#include <algorithm>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
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
#include "utils.hpp"
#include "io.hpp"
#include "rng.hpp"
#include "run_folders.hpp"
#include "helpers.hpp"
#include "interpol_dispatch.hpp"
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

  std::cout << "Creating folders..." << std::endl;
  RunFolders out = make_run_folders(intp->get_folder(), name, prm, folder_opts);

  if (prm.get<bool>("verbose"))
    prm.print();

  // Parallel generators
  std::vector<std::mt19937> gens = make_generators(prm);

  // TODO: These should not be stored in particle tracker parameters.
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

  intp->update(frozen_fields ? prm.get<double>("t_frozen") : t0);

  return Run{prm, intp, out, std::move(gens), t0, T, frozen_fields, false, t0};
}

// Load checkpoint or initial state
inline bool load_or_initialize(Run& run, Topology& mesh){
  const bool restarting = run.prm.get<std::string>("restart_folder") != "";
  if (restarting){
    mesh.load_checkpoint(run.prm.get<std::string>("restart_folder") + "/Checkpoints", run.prm);
  }
  else {
    std::shared_ptr<Initializer> init_state;
    set_initial_state(init_state, run.intp, run.prm, run.gens[0]);
    mesh.load_initial_state(init_state, run.prm);
  }
  mesh.compute_maps();
  return restarting;
}

// Particles that could not step: ignore, reinject, mark or remove
inline void handle_outside(Run& run, Topology& mesh, ParticleSet& ps, const std::string& outside,
                           const std::vector<Uint>& nodes, const double t, const bool verbose){
  if (nodes.empty())
    return;
  if (outside == "reinject"){
    const auto key = split_string(run.prm.get<std::string>("init_mode"), "_");
    const std::string dirs = key.size() > 1 ? key[1] : "xyz";
    const Vector3d Dx_max = 0.5*(run.intp->get_x_max() - run.intp->get_x_min());
    std::uniform_real_distribution<> ux(-Dx_max[0], Dx_max[0]), uy(-Dx_max[1], Dx_max[1]), uz(-Dx_max[2], Dx_max[2]);
    const bool rx = contains(dirs, "x"), ry = contains(dirs, "y"), rz = contains(dirs, "z");
    for (const Uint i : nodes){
      Vector3d Dx = {0., 0., 0.};
      do {
        if (rx) Dx[0] = ux(run.gens[0]);
        if (ry) Dx[1] = uy(run.gens[0]);
        if (rz) Dx[2] = uz(run.gens[0]);
      } while (!run.intp->locate(ps.x(i) + Dx));
      ps.set_x(i, ps.x(i) + Dx);
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

// App hooks
struct RunHooks {
  // After field update: injection, remeshing, resizing, removal
  std::function<void(int it, double t)> reshape = [](int, double){};
  // After step
  std::function<void(int it, double t, const std::vector<Uint>& outside)> after_step =
      [](int, double, const std::vector<Uint>&){};
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
  double t = restarting ? prm.get<double>("t") : run.t0;
  // Resume step count
  int it = restarting ? static_cast<int>(prm.get<Uint>("it")) : 0;

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
    { H5::H5File create(h5fname.c_str(), H5F_ACC_TRUNC); }
    h5f.openFile(h5fname.c_str(), H5F_ACC_RDWR);
  }

  std::ofstream statfile;
  if (stat_intv > 0.){
    statfile.open(newfolder + "/tdata_from_t" + std::to_string(t) + ".dat");
    write_stats_header(statfile, hooks.statistics ? hooks.statistics(0., stepper.counters())
                                                  : mesh.stats_header_columns(ds_max));
  }

  std::clock_t clock_0 = std::clock();
  while (t < run.T + dt/2 && hooks.keep_going()){
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
      with_concrete(*run.intp, [&](auto& ip){ ps.update_fields(ip, t_fields, output_fields); });
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
    }

    auto outside_nodes = stepper.step(*run.intp, ps, t, dt);

    hooks.after_step(it, t, outside_nodes);

    t += dt;
    it += 1;
  }
  std::clock_t clock_1 = std::clock();
  double duration = (clock_1-clock_0) / (double) CLOCKS_PER_SEC;
  std::cout << "Total simulation time: " << duration << " seconds" << std::endl;

  if (refresh_S){
    const double t_fields = run.marches ? run.t_fields : t;
    with_concrete(*run.intp, [&](auto& ip){ ps.update_fields(ip, t_fields, output_fields); });
  }
  prm.set<Uint>("it", it);
  mesh.write_checkpoint(run.out.checkpoints, t, prm);

  statfile.close();
}

#endif
