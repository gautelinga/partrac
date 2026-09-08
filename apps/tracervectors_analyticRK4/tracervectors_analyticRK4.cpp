#include <iostream>
#include <vector>
#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <sstream>
#include <random>
#include <cmath>
#include <set>
#include <iterator>
#include "H5Cpp.h"
//#include "hdf5.h"
#include <ctime>

#include "experimental/particles.hpp"
#include "utils.hpp"
#include "Params.hpp"
#include "rng.hpp"
#include "AnalyticInterpol.hpp"
#include "experimental/integrator_explicit.hpp"
#include "experimental/initializer.hpp"
#include "experimental/statistics.hpp"

std::string get_newfoldername(const std::string& rwfolder, const partrac::Params& prm){
  std::ostringstream ss_Dm, ss_dt, ss_Nrw, ss_seed;
  ss_Dm << std::scientific << std::setprecision(7) << prm.get<double>("Dm");
  ss_dt << std::scientific << std::setprecision(7) << prm.get<double>("dt");
  ss_Nrw << prm.get<Uint>("Nrw");
  ss_seed << prm.get<int>("seed");
  std::string newfoldername = rwfolder +
                            "/Dm" + ss_Dm.str() + // "_U" + std::to_string(prm.U0) +
                            "_dt" + ss_dt.str() +
                            "_Nrw" + ss_Nrw.str() +
                            "_seed" + ss_seed.str() +
                            prm.get<std::string>("tag") +
                            "/";
  return newfoldername;
}

#include "tracervectors_analyticRK4_schema.hpp"

int main(int argc, char* argv[])
{

    {
        std::cout << "======================================================================\n"
                  << "||  Initialized experimental tracer vectors.                        ||\n"
                  << "======================================================================" << std::endl;
    }

    // Input parameters
    if (argc < 2) {
        std::cout << "Please specify an input file." << std::endl;
        return 1;
    }
    partrac::Params prm = partrac::parse_or_exit(tracervectors_analyticRK4_schema(), argc, argv);
    if (prm.get<int>("num_threads") > 0){
        omp_set_dynamic(0);
        omp_set_num_threads(prm.get<int>("num_threads"));
    }

    std::string infilename = prm.input_file();

    std::cout << "Initializing AnalyticInterpol." << std::endl;
    AnalyticInterpol intp(infilename);

    intp.set_U0(prm.get<double>("U"));
    intp.set_int_order(2);  // To evaluate gradients

    std::string folder = intp.get_folder();
    std::string rwfolder = folder + "/TracerVectors/";
    
        create_folder(rwfolder);
    
    std::string newfolder;
    if (prm.get<std::string>("restart_folder") != ""){
        newfolder = prm.get<std::string>("folder");
    }
    else {
        newfolder = get_newfoldername(rwfolder, prm);
            create_folder(newfolder);
    }
    newfolder = newfolder + "" + "0" + "/";
    std::string posfolder = newfolder + "Positions/";
    std::string checkpointsfolder = newfolder + "Checkpoints/";
    {
        create_folder(newfolder);
        create_folder(posfolder);
        create_folder(checkpointsfolder);
    }
    prm.set<std::string>("folder", newfolder);

        if (prm.get<bool>("verbose")) prm.print();

    // Parallel generators
    std::vector<std::mt19937> gens = make_generators(prm);

    Real dt = prm.get<double>("dt");
    Real t0 = std::max(intp.get_t_min(), prm.get<double>("t0"));
    Real T = std::min(intp.get_t_max(), prm.get<double>("T"));
    prm.set<double>("t0", t0);
    prm.set<double>("T", T);

    // This part is unique
    std::cout << "Initializing Integrator..." << std::endl;
    // Integrator_RK4 integrator;
    Integrator_Explicit integrator(prm.get<double>("Dm"), 2, gens);

    std::cout << "Initializing ParticleSet..." << std::endl;
    Particles<Particle> ps(prm.get<Uint>("Nrw_max"));

    auto key = split_string(prm.get<std::string>("init_mode"), "_");
    if (key.size() == 0){
        std::cout << "init_mode not specified." << std::endl;
        exit(1);
    }

    experimental::RandomPointsInitializer init_state(key, prm, gens[0]);
    init_state.probe(intp);
    init_state.initialize(ps);
    spin_all(key, ps, gens[0]);

    // Check mesh connectivity: should be uneccessary
    ps.edges().clear();
    ps.faces().clear();

    int it = 0;
    Real t = t0;
    if (prm.get<std::string>("restart_folder") != ""){
        t = prm.get<double>("t");
    }
    prm.dump(newfolder, t);

    std::ofstream statfile;
    if (prm.get<double>("stat_intv") > 0.){
      statfile.open(newfolder + "/tdata_from_t" + std::to_string(t) + ".dat");
      write_stats_header(statfile, ps.dim());
    }

    std::string h5fname = newfolder + "/data_from_t" + std::to_string(t) + ".h5";
    // no dump file at all when dumping is off
    H5::H5File h5f;
    if (prm.get<double>("dump_intv") > 0.){
      { H5::H5File create(h5fname.c_str(), H5F_ACC_TRUNC); }
      h5f.openFile(h5fname.c_str(), H5F_ACC_RDWR);
    }

  const double chunk_intv = prm.get<double>("dump_intv")*prm.get<int>("dump_chunk_size");

    std::map<std::string, bool> output_fields;
    output_fields["u"] = !prm.get<bool>("minimal_output");
    output_fields["c"] = !prm.get<bool>("minimal_output");
    output_fields["p"] = !prm.get<bool>("minimal_output") && prm.get<bool>("output_all_props");
    output_fields["rho"] = !prm.get<bool>("minimal_output") && prm.get<bool>("output_all_props");        
    output_fields["H"] = !prm.get<bool>("minimal_output") && ps.dim() > 0;
    output_fields["n"] = true;
    output_fields["w"] = true;
    output_fields["S"] = true;
    output_fields["tau"] = true;
    output_fields["J"] = true;

    intp.update(t);
    intp.assign_fields(ps, output_fields);

    // Simulation start
    std::clock_t clock_0 = std::clock();

    while (t <= T){
        intp.update(t);

        // Update fields for output
        if (at_interval(it, prm.get<double>("dump_intv"), dt) || at_interval(it, prm.get<double>("stat_intv"), dt)){
            intp.assign_fields(ps, output_fields);
        }

        // Statistics
        if (at_interval(it, prm.get<double>("stat_intv"), dt)){
            std::cout << "Time = " << t << std::endl;
            write_stats(statfile, t, ps, integrator.get_declined());

            // intp.print_found();
        }
        // Checkpoint
        if (at_interval(it, prm.get<double>("checkpoint_intv"), dt)){
            //mesh.write_checkpoint(checkpointsfolder, t, prm);
        }

        // Dump detailed data
        if (at_interval(it, prm.get<double>("dump_intv"), dt)){
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
            ps.dump_hdf5(h5f, groupname, output_fields);
            h5f.close();
        }

        auto outside_nodes = integrator.step_vec(intp, ps, t, dt);

        if (outside_nodes.size() > 0){
            std::cout << outside_nodes.size() << " nodes are outside." << std::endl;
            //reinject_nodes(outside_nodes, key, ps, intp, gens[0]);
        }

        t += dt;
        ++it;
    }

    std::clock_t clock_1 = std::clock();
    Real duration = (clock_1-clock_0) / (Real) CLOCKS_PER_SEC;
    std::cout << "Total simulation time: " << duration << " seconds" << std::endl;

    return 0;
}
