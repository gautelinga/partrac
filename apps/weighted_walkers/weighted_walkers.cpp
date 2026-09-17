#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <vector>
#include <omp.h>
#include "files.hpp"
#include "H5Cpp.h"

#include "RunLoop.hpp"
#include "ExplicitIntegrator.hpp"
#include "Initializer.hpp"
#include "stats.hpp"

#include "weighted_walkers_schema.hpp"

// Exit plane axes: normal, tangent, tangent; -1 for none
inline std::array<int, 3> exit_axes(const std::string& exit_plane){
    if (exit_plane == "x") return {0, 1, 2};
    if (exit_plane == "y") return {1, 2, 0};
    if (exit_plane == "z") return {2, 0, 1};
    return {-1, -1, -1};
}

// Walkers beyond the exit plane, ascending
inline void get_exited_nodes(std::vector<Uint>& exited_nodes, std::vector<std::vector<Uint>>& buffers,
                             const ParticleSet& ps, const std::array<int, 3>& axes, const double Ln, const double Lt){
    exited_nodes.clear();
    if (axes[0] < 0)
        return;
    const int dn = axes[0], dt1 = axes[1], dt2 = axes[2];

    #pragma omp parallel
    {
        auto nthreads = omp_get_num_threads();
        auto id = omp_get_thread_num();

        #pragma omp single
        {
            buffers.resize( nthreads );
            for ( auto & buffer : buffers )
                buffer.clear();
        }

        #pragma omp for
        for ( Uint i = 0; i < ps.N(); ++i ){
            const Vector3d x = ps.x(i);
            if ( (x[dn] > Ln) || (Lt > 0. and (x[dt1]*x[dt1] + x[dt2]*x[dt2] > Lt*Lt)) ){
                buffers[id].push_back(i);
            }
        }

        #pragma omp single
        {
            for ( auto & buffer : buffers )
                exited_nodes.insert(exited_nodes.end(), buffer.begin(), buffer.end());
        }
    }
    std::sort(exited_nodes.begin(), exited_nodes.end());
}

// Replace exited walkers by splitting survivors (serial: a parent may be drawn twice)
inline bool split_random_nodes(const std::vector<Uint>& nodes_to_replace, std::vector<double>& weights,
                               ParticleSet& ps, std::vector<std::mt19937>& gens){
    weights.resize(ps.N());

    #pragma omp parallel for
    for ( Uint i = 0; i < ps.N(); ++i){
        weights[i] = std::ldexp(1.0, -static_cast<int>(ps.generation(i)));
    }
    for ( auto & i : nodes_to_replace ){
        weights[i] = 0.;
    }

    double total = 0.;
    #pragma omp parallel for reduction(+:total)
    for ( Uint i = 0; i < weights.size(); ++i)
        total += weights[i];
    if (total <= 0.)
        return false;

    std::discrete_distribution<std::mt19937::result_type> discrete_dist(weights.begin(), weights.end());
    auto & gen = gens[0];
    for ( auto & i : nodes_to_replace ){
        Uint j = discrete_dist(gen);
        ps.set_x(i, ps.x(j));
        ps.set_cell_id(i, ps.get_cell_id(j));
        const double generation = ps.generation(j) + 1;
        ps.set_generation(i, generation);
        ps.set_generation(j, generation);
    }
    return true;
}

// Separation data selection
struct SepdataSelection {
    Uint dim;
    int a = -1, b = -1;       // tangent axes, ascending; -1 for none
    bool use_a = false, use_b = false;
    double ds_max;
    Vector3d x0;
};

inline SepdataSelection sepdata_selection(const partrac::Params& prm, const std::vector<std::string>& key, const Uint dim){
    SepdataSelection sel;
    sel.dim = dim;
    sel.ds_max = prm.get<double>("ds_max");
    sel.x0 = {prm.get<double>("x0"), prm.get<double>("y0"), prm.get<double>("z0")};
    const std::string exit_plane = prm.get<std::string>("exit_plane");
    if (exit_plane == "x"){ sel.a = 1; sel.b = 2; }
    if (exit_plane == "y"){ sel.a = 0; sel.b = 2; }
    if (exit_plane == "z"){ sel.a = 0; sel.b = 1; }
    const char* axis = "xyz";
    if (sel.a >= 0){
        sel.use_a = contains(key[1], std::string(1, axis[sel.a]));
        sel.use_b = contains(key[1], std::string(1, axis[sel.b]));
    }
    return sel;
}

// Separation data
inline void write_separation_data(const std::string& folder, const double t, const ParticleSet& ps,
                                  const SepdataSelection& sel){
    const std::string sepdatafname = folder + "/sepdata_from_t" + std::to_string(t) + ".h5";
    H5::H5File sepdata_h5f(sepdatafname.c_str(), H5F_ACC_TRUNC);

    const std::string groupname = std::to_string(t);
    sepdata_h5f.createGroup(groupname + "/");

    std::vector<Vector3d> xyz_;
    std::vector<double> w_;
    if (sel.a >= 0){
        const double ds_max = sel.ds_max;
        for ( Uint i = 0; i < ps.N(); ++i ){
            const Vector3d x = ps.x(i);
            const double da = std::abs(x[sel.a]-sel.x0[sel.a]), db = std::abs(x[sel.b]-sel.x0[sel.b]);
            const bool pick = sel.dim == 2
                ? ((sel.use_a && da < ds_max) || (sel.use_b && db < ds_max))
                : (da*da + db*db < ds_max*ds_max);
            if (pick){
                xyz_.push_back(x);
                w_.push_back(ps.generation(i));
            }
        }
    }
    vector2hdf5(sepdata_h5f, groupname + "/x", xyz_, xyz_.size());
    scalar2hdf5(sepdata_h5f, groupname + "/w", w_, w_.size());

    sepdata_h5f.close();
}

int main(int argc, char* argv[])
{

    {
        std::cout << "======================================================================\n"
                  << "||  Initialized weighted walkers.                                   ||\n"
                  << "======================================================================" << std::endl;
    }

    // Input parameters
    if (argc < 2) {
        std::cout << "Please specify an input file." << std::endl;
        return 1;
    }
    partrac::Params prm = partrac::parse_or_exit(weighted_walkers_schema(), argc, argv);

    Run run = start_run(prm, "WeightedWalkers");
    const std::string sepdatafolder = run.out.run + "Sepdata/";
    create_folder(sepdatafolder);

    ExplicitIntegrator integrator(prm.get<double>("Dm"), prm.get<int>("int_order"), run.gens);

    // Generation per walker
    ParticleSet ps(run.intp, prm.get<Uint>("Nrw_max"));
    ps.record_generation();
    Topology mesh(ps, prm);

    // Gaussian strip or circle
    const std::vector<std::string> key = split_string(prm.get<std::string>("init_mode"), "_");
    const Uint dim = contains(key[0], "strip") ? 2 : 3;
    if (prm.get<std::string>("restart_folder") != ""){
        mesh.load_checkpoint(prm.get<std::string>("restart_folder") + "/Checkpoints", prm);
    }
    else {
        std::shared_ptr<Initializer> init_state = make_gaussian_initializer(dim, key, run.intp, prm, run.gens[0]);
        mesh.load_initial_state(init_state, prm);
    }
    mesh.compute_maps();

    std::map<std::string, bool> output_fields;
    output_fields["u"] = false;
    output_fields["c"] = !prm.get<bool>("minimal_output");
    output_fields["p"] = false;
    output_fields["rho"] = false;
    output_fields["H"] = false;
    output_fields["n"] = false;

    const double dt = prm.get<double>("dt");
    const double refine_intv = prm.get<double>("refine_intv");
    const std::array<int, 3> axes = exit_axes(prm.get<std::string>("exit_plane"));
    const SepdataSelection selection = sepdata_selection(prm, key, dim);
    const double Ln = prm.get<double>("Ln");
    const double Lt = prm.get<double>("Lt");

    // Stepper
    struct Stepper {
        ExplicitIntegrator& integrator;
        Integrator& counters(){ return integrator; }
        std::vector<Uint> step(Interpol& intp, ParticleSet& ps, const double t, const double dt){
            return explicit_step<TransportElement::Point>(integrator, intp, ps, t, dt);
        }
    } stepper{integrator};

    bool nothing_left = false;
    std::vector<Uint> exited_nodes;
    std::vector<std::vector<Uint>> exit_buffers;
    std::vector<double> weights;
    RunHooks hooks;
    hooks.statistics = [&](const double t, Integrator& counters){
        return cloud_stats_columns(t, ps, counters.get_declined());
    };
    hooks.after_statistics = [&](const int, const double t){
        write_separation_data(sepdatafolder, t, ps, selection);
    };
    // Resampling
    hooks.after_step = [&](const int it, const double t, const std::vector<Uint>&){
        if (!at_interval(it, refine_intv, dt))
            return;
        get_exited_nodes(exited_nodes, exit_buffers, ps, axes, Ln, Lt);
        if (exited_nodes.size() > 0 && !split_random_nodes(exited_nodes, weights, ps, run.gens)){
            std::cout << "Every walker has crossed the exit plane at t = " << t
                      << "; nothing is left to copy from. Stopping." << std::endl;
            nothing_left = true;
        }
    };
    hooks.keep_going = [&]{ return !nothing_left; };

    run_loop(run, ps, mesh, stepper, output_fields, dt, hooks);

    return 0;
}
