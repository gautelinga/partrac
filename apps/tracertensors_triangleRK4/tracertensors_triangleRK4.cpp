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
#include "MPIwrap.hpp"
#include "Parameters.hpp"
#include "StructuredInterpol.hpp"
#include "TetInterpol.hpp"
#include "experimental/integrator_RK.hpp"
#include "experimental/initializer.hpp"
#include "experimental/statistics.hpp"

std::string get_newfoldername(const std::string& rwfolder, const Parameters& prm){
  std::ostringstream ss_Dm, ss_dt, ss_Nrw, ss_seed;
  ss_Dm << std::scientific << std::setprecision(7) << prm.Dm;
  ss_dt << std::scientific << std::setprecision(7) << prm.dt;
  ss_Nrw << prm.Nrw;
  ss_seed << prm.seed;
  std::string newfoldername = rwfolder +
                            "/Dm" + ss_Dm.str() + // "_U" + std::to_string(prm.U0) +
                            "_dt" + ss_dt.str() +
                            "_Nrw" + ss_Nrw.str() +
                            "_seed" + ss_seed.str() +
                            prm.tag +
                            "/";
  return newfoldername;
}

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
            if (contains(key[1], "x")){
                Dx[0] = uni_dist_x(gen);
            }
            if (contains(key[1], "y")){
                Dx[1] = uni_dist_y(gen);
            }
            if (contains(key[1], "z")){
                Dx[2] = uni_dist_z(gen);
            }
            Vector x0 = node.x();
            intp.probe(x0 + Dx);
            outside = !intp.inside_domain();
        }
        node.x() += Dx;
    }
}

int main(int argc, char* argv[])
{
    MPIwrap mpi(argc, argv);

    if (mpi.rank() == 0)
    {
        std::cout << "======================================================================\n"
                  << "||  Initialized experimental tracer tensors                         ||\n"
                  << "======================================================================" << std::endl;
    }
    mpi.barrier();

    // Input parameters
    if (argc < 2 && mpi.rank() == 0) {
        std::cout << "Please specify an input file." << std::endl;
        return 0;
    }
    Parameters prm(argc, argv);
    if (prm.restart_folder != ""){
        prm.parse_file(prm.restart_folder + "/Checkpoints/params.dat");
        prm.parse_cmd(argc, argv);
    }

    std::string infilename = std::string(argv[1]);

    std::cout << "Initializing TetInterpol." << std::endl;
    TetInterpol intp(infilename);
    // std::cout << "Initialized TriangleInterpol." << std::endl;

    intp.set_U0(prm.U0);
    intp.set_int_order(2);  // To evaluate gradients

    std::string folder = intp.get_folder();
    std::string rwfolder = folder + "/TracerTensors/";
    
    if (mpi.rank() == 0)
        create_folder(rwfolder);
    
    std::string newfolder;
    if (prm.restart_folder != ""){
        newfolder = prm.folder;
    }
    else {
        newfolder = get_newfoldername(rwfolder, prm);
        mpi.barrier();
        if (mpi.rank() == 0)
            create_folder(newfolder);
        mpi.barrier();
    }
    newfolder = newfolder + "" + std::to_string(mpi.rank()) + "/";
    std::string posfolder = newfolder + "Positions/";
    std::string checkpointsfolder = newfolder + "Checkpoints/";
    {
        create_folder(newfolder);
        create_folder(posfolder);
        create_folder(checkpointsfolder);
    }
    prm.folder = newfolder;

    if (mpi.rank() == 0)
        prm.print();

    std::mt19937 gen;
    if (prm.random) {
        std::random_device rd;
        gen.seed(rd());
    }
    else {
        std::seed_seq rd{prm.seed + mpi.rank()};
        gen.seed(rd);
    }

    Real dt = prm.dt;
    Real t0 = std::max(intp.get_t_min(), prm.t0);
    Real T = std::min(intp.get_t_max(), prm.T);
    prm.t0 = t0;
    prm.T = T;

    // This part is unique
    std::cout << "Initializing Integrator..." << std::endl;
    Integrator_RK4 integrator;

    std::cout << "Initializing ParticleSet..." << std::endl;
    Particles<Particle> ps(prm.Nrw_max);

    auto key = split_string(prm.init_mode, "_");
    if (key.size() == 0){
        std::cout << "init_mode not specified." << std::endl;
        exit(0);
    }

    RandomPointsInitializer init_state(key, prm, gen);
    init_state.probe(intp);
    init_state.initialize(ps);

    // Check mesh connectivity: should be uneccessary
    ps.edges().clear();
    ps.faces().clear();

    int it = 0;
    Real t = t0;
    if (prm.restart_folder != ""){
        t = prm.t;
    }
    prm.dump(newfolder, t);

    std::ofstream statfile(newfolder + "/tdata_from_t" + std::to_string(t) + ".dat");
    write_stats_header(mpi, statfile, ps.dim());

    std::string h5fname = newfolder + "/data_from_t" + std::to_string(t) + ".h5";
    H5::H5File h5f(h5fname.c_str(), H5F_ACC_TRUNC);

    Uint int_stat_intv = int(prm.stat_intv/dt);
    Uint int_dump_intv = int(prm.dump_intv/dt);
    Uint int_checkpoint_intv = int(prm.checkpoint_intv/dt);
    Uint int_chunk_intv = int_dump_intv*prm.dump_chunk_size;

    std::map<std::string, bool> output_fields;
    output_fields["u"] = !prm.minimal_output;
    output_fields["c"] = !prm.minimal_output;
    output_fields["p"] = !prm.minimal_output && prm.output_all_props;
    output_fields["rho"] = !prm.minimal_output && prm.output_all_props;        
    output_fields["H"] = !prm.minimal_output && ps.dim() > 0;
    output_fields["n"] = false;
    output_fields["w"] = false;
    output_fields["S"] = false;
    output_fields["F"] = true;

    intp.update(t);
    intp.assign_fields(ps, output_fields);

    // Simulation start
    std::clock_t clock_0 = std::clock();

    while (t <= T){
        intp.update(t);

        // Update fields for output
        if (it % int_dump_intv == 0 || it % int_stat_intv == 0){
            intp.assign_fields(ps, output_fields);
        }

        // Statistics
        if (it % int_stat_intv == 0){
            std::cout << "Time = " << t << std::endl;
            write_stats(mpi, statfile, t, ps, integrator.get_declined());

            intp.print_found();
        }
        // Checkpoint
        if (it % int_checkpoint_intv == 0){
            //mesh.write_checkpoint(checkpointsfolder, t, prm);
        }

        // Dump detailed data
        if (it % int_dump_intv == 0){
            std::string groupname = std::to_string(t);

            // Clear file if it exists, otherwise create
            if (int_chunk_intv > 0 && it % int_chunk_intv == 0 && it > 0){
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

        // auto outside_nodes = integrator.step_vec(intp, ps, t, dt);
        
        auto outside_nodes = integrator.step_tensor(intp, ps, t, dt);

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
