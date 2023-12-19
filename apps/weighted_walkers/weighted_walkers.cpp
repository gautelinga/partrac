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
#include <omp.h>
#include <chrono>

#include "experimental/particles.hpp"
#include "utils.hpp"
#include "MPIwrap.hpp"
#include "Parameters.hpp"
#include "AnalyticInterpol.hpp"
#include "experimental/integrator_explicit.hpp"
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

template<typename T>
std::vector<Uint> get_exited_nodes(T& ps, const std::string& exit_plane, const double Ln, const double Lt){
    std::vector<Uint> exited_nodes;
    Uint dn, dt1, dt2;
    if (exit_plane == "x"){
        dn = 0;
        dt1 = 1;
        dt2 = 2;
    }
    else if (exit_plane == "y"){
        dn = 1;
        dt1 = 2;
        dt2 = 0;
    }
    else if (exit_plane == "z"){
        dn = 2;
        dt1 = 0;
        dt2 = 1;
    }
    else {
        return exited_nodes;
    }

    std::vector<std::vector<Uint>> buffers;
    
    #pragma omp parallel
    {
        auto nthreads = omp_get_num_threads();
        auto id = omp_get_thread_num();

        #pragma omp single
        {
            buffers.resize( nthreads );
        }

        #pragma omp for
        for ( Uint i = 0; i < ps.particles().size(); ++i ){
            Vector x = ps.particles()[i].x();
            if ( (x[dn] > Ln) || (Lt > 0. and (pow(x[dt1], 2) + pow(x[dt2], 2) > Lt*Lt)) ){
                //exited_nodes.insert(i);
                buffers[id].push_back(i);
            }
        }

        #pragma omp single
        {
            for ( auto & buffer : buffers ) {
                move(buffer.begin(), buffer.end(), std::back_inserter(exited_nodes));
                //exited_nodes.insert(buffer.begin(), buffer.end());
            }
            //exited_nodes.insert(vec.begin(), vec.end());
        }
    }
    std::sort(exited_nodes.begin(), exited_nodes.end());
    return exited_nodes;
}

template<typename T>
void split_random_nodes(std::vector<Uint>& nodes_to_replace, T& ps, std::vector<std::mt19937>& gens, Parameters& prm){

    //auto t0 = std::chrono::high_resolution_clock::now();

    //std::vector<double> nodes_to_replace_vec(nodes_to_replace.begin(), nodes_to_replace.end());

    /*
    std::set<Uint> all_nodes;
    for (Uint i = 0; i < ps.particles().size(); ++i)
    {
        all_nodes.emplace_hint(all_nodes.end(), i);
    }
    */

    //auto t1 = std::chrono::high_resolution_clock::now();

    /*
    std::set<Uint> good_nodes(all_nodes.begin(), all_nodes.end());
    for ( auto & i : nodes_to_replace ){
        good_nodes.erase(i);
    }
    */

    /*
    std::set_difference(all_nodes.begin(), all_nodes.end(), 
                        nodes_to_replace.begin(), nodes_to_replace.end(),
                        std::inserter(good_nodes, good_nodes.end()));
    */


    /*
    std::vector<Uint> good_nodes_vec(good_nodes.begin(), good_nodes.end());
    */

    std::vector<double> weights(ps.particles().size());

    //auto t2 = std::chrono::high_resolution_clock::now();

    //for ( auto &i : good_nodes_vec ){
    #pragma omp parallel for
    for ( Uint i = 0; i < ps.particles().size(); ++i){
        //double weight = pow(2, -ps.particles()[i].w()); // 1
        double weight = pow(2, -ps.particles()[i].w()); // 1
        // double weight = ps.particles()[i].w();
        weights[i] = weight;
    }

    //auto t3 = std::chrono::high_resolution_clock::now();
    
    #pragma omp parallel for
    for ( auto & i : nodes_to_replace ){
        weights[i] = 0.;
    }
    
    //std::cout << "Weights: " << weights.size() << std::endl;

    //auto t4 = std::chrono::high_resolution_clock::now();

    #pragma omp parallel 
    {
        std::discrete_distribution<std::mt19937::result_type> discrete_dist(weights.begin(), weights.end());
        auto & gen = gens[omp_get_thread_num()];

        #pragma omp for
        for ( auto & i : nodes_to_replace ){
            //Uint j = good_nodes_vec[discrete_dist(gen)];
            Uint j = discrete_dist(gen);
            ps.particles()[i].x() = ps.particles()[j].x();
            double w_new = ps.particles()[j].w()+1; // ps.particles()[j].w()/2;
            ps.particles()[i].w() = w_new;
            ps.particles()[j].w() = w_new;
        }
    }

    //auto t5 = std::chrono::high_resolution_clock::now();

    //auto dt1 = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0);
    //auto dt2 = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1);
    //auto dt3 = std::chrono::duration_cast<std::chrono::microseconds>(t3-t2);
    //auto dt4 = std::chrono::duration_cast<std::chrono::microseconds>(t4-t3);
    //auto dt5 = std::chrono::duration_cast<std::chrono::microseconds>(t5-t4);

    //std::cout << dt1.count() << " " << dt2.count() << " " << dt3.count() << " " << dt4.count() << " " << dt5.count() << std::endl;
}

int main(int argc, char* argv[])
{
    MPIwrap mpi(argc, argv);

    if (mpi.rank() == 0)
    {
        std::cout << "======================================================================\n"
                  << "||  Initialized weighted walkers.                                   ||\n"
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

    if (prm.num_threads > 0){
        omp_set_dynamic(0);
        omp_set_num_threads(prm.num_threads);
    }

    std::string infilename = std::string(argv[1]);

    AnalyticInterpol intp(infilename);

    intp.set_U0(prm.U0);
    intp.set_int_order(prm.int_order);

    std::string folder = intp.get_folder();
    std::string rwfolder = folder + "/WeightedWalkers/";
    
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
    std::string sepdatafolder = newfolder + "Sepdata/";
    {
        create_folder(newfolder);
        create_folder(posfolder);
        create_folder(checkpointsfolder);
        create_folder(sepdatafolder);
    }
    prm.folder = newfolder;

    if (mpi.rank() == 0)
        prm.print();

    // Parallel generators
    std::vector<std::mt19937> gens;
    for (int i=0, N=omp_get_max_threads(); i<N; ++i) {
        std::mt19937 gen;
        if (prm.random) {
            std::random_device rd;
            gen.seed(rd());
        }
        else {
            std::seed_seq rd{prm.seed + omp_get_thread_num() };
            gen.seed(rd);
        }
        gens.emplace_back(gen);
    }

    std::uniform_int_distribution<std::mt19937::result_type> uniform_dist(0, prm.Nrw);

    Real dt = prm.dt;
    Real t0 = std::max(intp.get_t_min(), prm.t0);
    Real T = std::min(intp.get_t_max(), prm.T);
    prm.t0 = t0;
    prm.T = T;

    // This part is unique
    std::cout << "initializing Integrator..." << std::endl;
    Integrator_Explicit integrator(prm.Dm, prm.int_order, gens);

    std::cout << "Initializing ParticleSet..." << std::endl;
    Particles<Particle> ps(prm.Nrw_max);

    auto key = split_string(prm.init_mode, "_");
    if (key.size() == 0){
        std::cout << "init_mode not specified." << std::endl;
        exit(0);
    }

    //std::shared_ptr<Interpol> intp_ptr (&intp);
    //std::shared_ptr<Initializer> init_state;
    //init_state = std::make_shared<RandomPointsInitializer>(key, intp_ptr, prm, mpi, gen);    
    //init_state->initialize(ps);
    //RandomPointsInitializer init_state(key, prm, mpi, gen);
    //std::shared_ptr<Initializer> init_state;
    Uint dim;
    if (contains(key[0], "strip")){
        dim = 2;
        RandomGaussianStripInitializer init_state(key, prm, gens[0]);
        init_state.probe(intp);
        init_state.initialize(ps);
    }
    else if (contains(key[0], "circle")){
        dim = 3;
        RandomGaussianCircleInitializer init_state(key, prm, gens[0]);
        init_state.probe(intp);
        init_state.initialize(ps);
    }
    else {
        std::cout << "Unknown initial state: " << key[0] << "." << std::endl;
        exit(0);
    }

    //RandomGaussianCircleInitializer init_state(key, prm, mpi, gen);
    
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
    Uint int_reinject_intv = int(prm.refine_intv/dt);
    Uint int_checkpoint_intv = int(prm.checkpoint_intv/dt);
    Uint int_chunk_intv = int_dump_intv*prm.dump_chunk_size;

    std::map<std::string, bool> output_fields;
    output_fields["u"] = false; // !prm.minimal_output;
    output_fields["c"] = !prm.minimal_output;
    output_fields["p"] = false; // !prm.minimal_output && prm.output_all_props;
    output_fields["rho"] = false;  // !prm.minimal_output && prm.output_all_props;        
    output_fields["H"] = false;  //& !prm.minimal_output && ps.dim() > 0;
    output_fields["n"] = false; // !prm.minimal_output && ps.dim() > 1;
    output_fields["w"] = true;

    intp.update(t);
    intp.assign_fields(ps, output_fields);
    for ( auto & particle : ps.particles() ){
        //particle.w() = 1.0;
        particle.w() = 0.0;
    }

    // Simulation start
    auto clock_0 = std::clock();

    double duration_par = 0;
    double duration_split = 0;

    while (t <= T){
        intp.update(t);

        // Update fields for output
        if (it % int_dump_intv == 0 || it % int_stat_intv == 0){
            intp.assign_fields(ps, output_fields);
        }

        // Statistics
        if (it % int_stat_intv == 0){
            std::cout << "Time = " << t << " [" << duration_par << " + " << duration_split << "]" << std::endl;
            duration_par = 0;
            duration_split = 0;
            write_stats(mpi, statfile, t, ps, integrator.get_declined());

            std::string sepdatafname = sepdatafolder + "/sepdata_from_t" + std::to_string(t) + ".h5";
            H5::H5File sepdata_h5f(sepdatafname.c_str(), H5F_ACC_TRUNC);

            std::string groupname = std::to_string(t);
            sepdata_h5f.createGroup(groupname + "/");

            std::vector<double> xyz_;
            std::vector<double> w_;
            xyz_.reserve(ps.particles().size() * 3);
            w_.reserve(ps.particles().size());
            for ( auto & particle : ps.particles() ){
                Vector3d x = particle.x();
                Vector3d ds = {abs(x[0]-prm.x0), abs(x[1]-prm.y0), abs(x[2]-prm.z0)};
                if (dim == 2 && (
                    (prm.exit_plane == "x" && (
                     (contains(key[1], "y") && ds[1] < prm.ds_max) || 
                     (contains(key[1], "z") && ds[2] < prm.ds_max)
                    )) ||
                    (prm.exit_plane == "y" && (
                     (contains(key[1], "x") && ds[0] < prm.ds_max) ||
                     (contains(key[1], "z") && ds[0] < prm.ds_max)
                    ))) || 
                   dim == 3 && (
                    (prm.exit_plane == "x" && (std::pow(ds[1], 2) + std::pow(ds[2], 2) < std::pow(prm.ds_max, 2))) ||
                    (prm.exit_plane == "y" && (std::pow(ds[0], 2) + std::pow(ds[2], 2) < std::pow(prm.ds_max, 2))) ||
                    (prm.exit_plane == "z" && (std::pow(ds[0], 2) + std::pow(ds[1], 2) < std::pow(prm.ds_max, 2)))
                   )){
                    xyz_.push_back(x[0]);
                    xyz_.push_back(x[1]);
                    xyz_.push_back(x[2]);
                    w_.push_back(particle.w());
                }
            }
            vector_to_h5(sepdata_h5f, groupname + "/x", xyz_, 3);
            scalar_to_h5(sepdata_h5f, groupname + "/w", w_);

            sepdata_h5f.close();
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

        auto ct0 = std::chrono::high_resolution_clock::now();
        //auto outside_nodes = integrator.step(intp, ps, t, dt);
        integrator.step_parallel(intp, ps, t, dt);
        auto ct1 = std::chrono::high_resolution_clock::now();
        auto dct10 = std::chrono::duration_cast<std::chrono::microseconds>(ct1-ct0);
        duration_par += dct10.count();

        //if (outside_nodes.size() > 0){
        //    //std::cout << outside_nodes.size() << " nodes are outside." << std::endl;
        //}
        if (it % int_reinject_intv == 0){
            auto exited_nodes = get_exited_nodes(ps, prm.exit_plane, prm.Ln, prm.Lt);
            //std::cout << exited_nodes.size() << " nodes have crossed the " << prm.exit_plane << " plane." << std::endl;
            if (exited_nodes.size() > 0){
                //
                split_random_nodes(exited_nodes, ps, gens, prm);
            }
            auto ct2 = std::chrono::high_resolution_clock::now();
            auto dct21 = std::chrono::duration_cast<std::chrono::microseconds>(ct2-ct1);
            duration_split += dct21.count();
        }

        t += dt;
        ++it;
    }

    auto clock_1 = std::clock();
    Real duration = (clock_1-clock_0) / (Real) CLOCKS_PER_SEC;
    std::cout << "Total simulation time: " << duration << " seconds" << std::endl;

    return 0;
}
