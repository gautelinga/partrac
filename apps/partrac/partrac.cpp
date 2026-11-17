#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>

#include "Error.hpp"
#include "wait_policy.hpp"
#include "param_print.hpp"
#include "RunLoop.hpp"
#include "TimeScheme.hpp"

#include "partrac_schema.hpp"
#include "crossings.hpp"

// A generator's state as a text attribute of a checkpoint
static void write_generator(H5::H5File& h5f, const std::string& name, const std::mt19937& gen){
  std::ostringstream os;
  os << gen;
  const H5::StrType text(H5::PredType::C_S1, os.str().size());
  h5f.createAttribute(name, text, H5::DataSpace(H5S_SCALAR)).write(text, os.str());
}

// The state written so, if the checkpoint at path has it
static void read_generator(const std::string& path, const std::string& name, std::mt19937& gen){
  if (!std::ifstream(path))
    return;
  try {
    H5::H5File h5f(path, H5F_ACC_RDONLY);
    if (H5Aexists(h5f.getId(), name.c_str()) <= 0)
      return;
    H5::Attribute a = h5f.openAttribute(name);
    std::string state;
    a.read(a.getStrType(), state);
    std::istringstream is(state);
    is >> gen;
    if (!is)
      partrac::fail(path, ": '", name, "' is no generator state");
  } catch (const H5::Exception&){
    partrac::fail("cannot read '", name, "' in ", path);
  }
}

static int run(int argc, char* argv[])
{

    std::cout << "Initialized Partrac." << std::endl;

  partrac::Params prm = partrac::parse_or_exit(partrac_schema(), argc, argv);

  Run run = start_run(prm, "RandomWalkers");
  TimeScheme scheme(prm, run.gens);

  ParticleSet ps(run.intp, prm.get<Uint>("Nrw_max"));
  scheme.prepare(*run.intp);
  if (prm.get<bool>("output_phi"))
    ps.record_phi();
  Topology mesh(ps, mesh_options(prm));

  if (prm.get<bool>("inject"))
    std::cout << "Injection activated!" << std::endl;

  const bool restarting = load_or_initialize(run, mesh);

  const bool filter = prm.get<bool>("filter");
  if (filter && !mesh.faces.empty())
    partrac::fail("filter=true filters a strip; a sheet is not supported");
  // The filter's draws: a sequence of their own, so the thread count and the
  // diffusive noise leave them alone; checkpointed, so a restart goes on with
  // them (an older checkpoint has none: seeded afresh)
  std::mt19937 filter_gen = stream_generator(prm, 1);
  if (filter && restarting)
    read_generator(run.out.checkpoints + "checkpoint.h5", "filter_rng", filter_gen);
  const bool crossing = crossings_on(prm);
  if (crossing && !mesh.faces.empty())
    partrac::fail("cross_intv records strips and particles; a sheet is not supported");
  if (crossing && prm.get<std::string>("cross_origin") == "start" && !mesh.edges.empty())
    partrac::fail("cross_origin=start follows points and line elements: a strip's new nodes have no start of "
                  "their own");

  if (check_only(run, ps, mesh))
    return 0;

  const bool inject = prm.get<bool>("inject");
  // A restart from a checkpoint written after its step's reshaping repeats
  // that step: refinement, coarsening and filtering change nothing then, an
  // injection of points would inject again
  const int injected_at = restarting && prm.get<bool>("reshaped") ? static_cast<int>(prm.get<Uint>("it")) : -1;
  const bool verbose = prm.get<bool>("verbose");

  const double dt = prm.get<double>("dt");
  Remesher remesher(mesh, run, dt);
  remesher.initial(restarting, inject);

  if (verbose){
    print_param("sqrt(2*Dm*dt)", sqrt(2*prm.get<double>("Dm")*dt));
    print_param("U*dt         ", prm.get<double>("U")*dt);
  }

  OutputFields output_fields = mesh_output_fields(prm, mesh);
  // rho whatever output_all_props says
  output_fields.rho = !prm.get<bool>("minimal_output");
  output_fields.tau = prm.get<bool>("integrate_tau");
  output_fields.phi = prm.get<bool>("output_phi");

  // Hook constants
  const bool integrate_tau = prm.get<bool>("integrate_tau");
  const double inject_intv = prm.get<double>("inject_intv");
  const double T_inject = prm.get<double>("T_inject");
  const double filter_intv = prm.get<double>("filter_intv");
  const double tau_intv = prm.get<double>("tau_intv");
  const double tau_max = prm.get<double>("tau_max");
  const double Ln = prm.get<double>("Ln");
  const std::string outside = prm.get<std::string>("outside");

  const std::string exit_plane = prm.get<std::string>("exit_plane");
  const bool any_exit_plane = exit_plane != "none";
  // One normal for the exit plane and the crossing planes: an axis, or cross_n
  const std::string cross_n = prm.get<std::string>("cross_n");
  const Vector3d plane_n = parse_normal(cross_n.empty() ? exit_plane : cross_n);
  std::unique_ptr<Crossings> crossings;
  if (crossing)
    crossings = std::make_unique<Crossings>(ps, mesh.edges.empty() ? CrossingRecords::Points : CrossingRecords::Edges,
                                            plane_n, prm.get<std::string>("cross_intv"),
                                            prm.get<std::string>("cross_origin"),
                                            run.out.run + "/" + prm.get<std::string>("cross_file"),
                                            restarting ? run.out.checkpoints + "checkpoint.h5" : "");

  RunHooks hooks;

  // Reshaping
  hooks.reshape = [&](const int it, const double t, const bool dumping){
    // Injection
    if (inject && it > 0 && it != injected_at && at_interval(it, inject_intv, dt) && t <= T_inject){
      mesh.inject();
    }
    remesher.reshape(it, dumping);
    // Filtering
    if (filter && at_interval(it, filter_intv, dt)){
      const Uint before = mesh.edges.size();
      bool filtered = mesh.filter(filter_gen);
      if (crossings)
        crossings->filtered(t, before, mesh.edges.size());
      if (verbose && filtered)
        std::cout << "Filtered edges." << std::endl;
    }
    // Removal
    if (any_exit_plane && at_interval(it, filter_intv, dt)){
      Uint n_rem = mesh.remove_beyond(plane_n, Ln);
      if (verbose && n_rem > 0)
        std::cout << "Removed " << n_rem << " nodes that were beyond." << std::endl;
    }
    // The starts of the nodes injected
    if (crossings)
      crossings->place(ps);
  };

  // Outside nodes
  hooks.outside = [&](const std::vector<Uint>& outside_nodes, const double t){
    handle_outside(mesh, ps, outside, outside_nodes, t, verbose);
  };
  // Tau integration
  hooks.after_step = [&](const int it, const double t){
    if (integrate_tau && at_interval(it + 1, tau_intv, dt)){
      mesh.integrate_tau(dt * steps_per(tau_intv, dt), tau_max);
    }
    if (crossings)
      crossings->record(ps, mesh.edges, mesh.node2edges, t, dt);
  };
  if (crossings || filter)
    hooks.checkpoint = [&](H5::H5File& h5f){
      if (crossings)
        crossings->checkpoint(h5f, ps);
      if (filter)
        write_generator(h5f, "filter_rng", filter_gen);
    };

  OutsideTotals totals{outside, *run.intp, [&]{ return scheme.counters().get_declined(); }};
  totals.attach(hooks, verbose);

  run_loop(run, ps, mesh, scheme, output_fields, dt, hooks);
  scheme.report(std::cout);
  totals.finish(std::cout);

  return 0;
}

int main(int argc, char* argv[])
{
  partrac::wait_passively(argv);
  return partrac::report_errors([&]{ return run(argc, argv); });
}
