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
#include "run_folders.hpp"
#include "StructuredInterpol.hpp"
#include "TriangleInterpol.hpp"
#include "experimental/integrator_RK.hpp"
#include "experimental/initializer.hpp"
#include "experimental/statistics.hpp"


template<typename ParticleType, typename InterpolatorType>
void reinject_nodes( const std::set<Uint>& outside_node_ids
                   , const std::vector<std::string>& key
                   , Particles<ParticleType>& ps
                   , InterpolatorType& intp
                   , std::mt19937 &gen){
    assert(outside_node_ids.size() > 0);

    Vector Dx_max = 0.5*(intp.get_x_max()-intp.get_x_min());
    std::uniform_real_distribution<> uni_dist_x(-Dx_max[0], Dx_max[0]);
    std::uniform_real_distribution<> uni_dist_y(-Dx_max[1], Dx_max[1]);
    std::uniform_real_distribution<> uni_dist_z(-Dx_max[2], Dx_max[2]);
    for ( auto node_id : outside_node_ids ){
        auto & node = ps.particles()[node_id];

        bool outside = true;
        Vector Dx = {0., 0., 0.};
        while (outside)
        { 
            if (contains(key[2], "x")){
                Dx[0] = uni_dist_x(gen);
            }
            if (contains(key[2], "y")){
                Dx[1] = uni_dist_y(gen);
            }
            if (contains(key[2], "z")){
                Dx[2] = uni_dist_z(gen);
            }
            Vector x0 = node.x();
            outside = !intp.locate(x0 + Dx);
        }
        node.x() += Dx;
    }
}

#include "tracers_triangleRK4_schema.hpp"

int main(int argc, char* argv[])
{

    {
        std::cout << "======================================================================\n"
                  << "||  Initialized experimental tracers.\t\t\t\t\t ||\n"
                  << "======================================================================" << std::endl;
    }
    

    // Input parameters
    if (argc < 2) {
        std::cout << "Please specify an input file." << std::endl;
        return 1;
    }
    partrac::Params prm = partrac::parse_or_exit(tracers_triangleRK4_schema(), argc, argv);

    std::string infilename = prm.input_file();

    TriangleInterpol intp(infilename);

    intp.set_U0(prm.get<double>("U"));
    intp.set_int_order(prm.get<int>("int_order"));

    std::string folder = intp.get_folder();
    RunFolders out = make_run_folders(folder, "Tracers", prm);
    const std::string& newfolder = out.run;

        if (prm.get<bool>("verbose")) prm.print();

    std::mt19937 gen;
    if (prm.get<bool>("random")) {
        std::random_device rd;
        gen.seed(rd());
    }
    else {
        std::seed_seq rd{prm.get<int>("seed") + 0};
        gen.seed(rd);
    }

    Real dt = prm.get<double>("dt");
    Real t0 = std::max(intp.get_t_min(), prm.get<double>("t0"));
    Real T = std::min(intp.get_t_max(), prm.get<double>("T"));
    prm.set<double>("t0", t0);
    prm.set<double>("T", T);

    // This part is unique
    std::cout << "initializing Integrator..." << std::endl;
    Integrator_RK4 integrator;
    //Integrator_Explicit integrator(prm.Dm, prm.int_order, gen);

    std::cout << "Initializing ParticleSet..." << std::endl;
    Particles<Particle> ps(prm.get<Uint>("Nrw_max"));

    auto key = split_string(prm.get<std::string>("init_mode"), "_");
    if (key.size() == 0){
        std::cout << "init_mode not specified." << std::endl;
        exit(1);
    }

    experimental::RandomPointsInitializer init_state(key, prm, gen);
    init_state.probe(intp);
    init_state.initialize(ps);

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
    output_fields["n"] = !prm.get<bool>("minimal_output") && ps.dim() > 1;

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
        }
        // Checkpoint
        if (at_interval(it, prm.get<double>("checkpoint_intv"), dt)){
            //mesh.write_checkpoint(out.checkpoints, t, prm);
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

        auto outside_nodes = integrator.step(intp, ps, t, dt);

        if (outside_nodes.size() > 0){
            std::cout << outside_nodes.size() << " nodes are outside." << std::endl;
            reinject_nodes(outside_nodes, key, ps, intp, gen);
        }

        t += dt;
        ++it;
    }

    std::clock_t clock_1 = std::clock();
    Real duration = (clock_1-clock_0) / (Real) CLOCKS_PER_SEC;
    std::cout << "Total simulation time: " << duration << " seconds" << std::endl;

    return 0;
}
