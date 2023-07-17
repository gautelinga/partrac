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
std::set<Uint> get_exited_nodes(T& ps, const std::string& exit_plane, const double Ln, const double Lt){
    std::set<Uint> exited_nodes;
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
    Uint i = 0;
    for ( auto & particle : ps.particles() ){
        Vector x = particle.x();
        if ( (x[dn] > Ln) || (Lt > 0. and (pow(x[dt1], 2) + pow(x[dt2], 2) > Lt*Lt)) ){
            exited_nodes.insert(i);
        }
        ++i;
    }
    return exited_nodes;
}

template<typename T>
void split_random_nodes(std::set<Uint>& nodes_to_replace, T& ps, std::mt19937& gen, Parameters& prm){
    std::set<Uint> good_nodes;
    for (Uint i = 0; i < ps.particles().size(); ++i)
    {
        good_nodes.emplace_hint(good_nodes.end(), i);
    }
    for ( auto & i : nodes_to_replace ){
        good_nodes.erase(i);
    }
    std::vector<Uint> good_nodes_vec(good_nodes.begin(), good_nodes.end());
    std::vector<double> weights;
    for ( auto &i : good_nodes_vec ){
        double weight = pow(2, -ps.particles()[i].w()); // 1
        // double weight = ps.particles()[i].w();
        weights.push_back(weight);
    }
    std::discrete_distribution<std::mt19937::result_type> discrete_dist(weights.begin(), weights.end());
    for ( auto & i : nodes_to_replace ){
        Uint j = good_nodes_vec[discrete_dist(gen)];
        ps.particles()[i].x() = ps.particles()[j].x();
        double w_new = ps.particles()[j].w()+1; // ps.particles()[j].w()/2;
        ps.particles()[i].w() = w_new;
        ps.particles()[j].w() = w_new;
    }
}

int main(int argc, char* argv[])
{
    MPIwrap mpi(argc, argv);

    if (mpi.rank() == 0)
    {
        std::cout << "======================================================================\n"
                  << "||  Initialized experimental tracers with " << mpi.size() << " processes. \t\t\t ||\n"
                  << "======================================================================" << std::endl;
    }
    mpi.barrier();
    
    std::cout << "This is process " << mpi.rank() << " out of " << mpi.size() << "." << std::endl;
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

    std::mt19937 gen;
    if (prm.random) {
        std::random_device rd;
        gen.seed(rd());
    }
    else {
        std::seed_seq rd{prm.seed + mpi.rank()};
        gen.seed(rd);
    }
    std::uniform_int_distribution<std::mt19937::result_type> uniform_dist(0, prm.Nrw);

    Real dt = prm.dt;
    Real t0 = std::max(intp.get_t_min(), prm.t0);
    Real T = std::min(intp.get_t_max(), prm.T);
    prm.t0 = t0;
    prm.T = T;

    // This part is unique
    std::cout << "initializing Integrator..." << std::endl;
    Integrator_Explicit integrator(prm.Dm, prm.int_order, gen);

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
        RandomGaussianStripInitializer init_state(key, prm, mpi, gen);
        init_state.probe(intp);
        init_state.initialize(ps);
    }
    else if (contains(key[0], "circle")){
        dim = 3;
        RandomGaussianCircleInitializer init_state(key, prm, mpi, gen);
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

        auto outside_nodes = integrator.step(intp, ps, t, dt);

        if (outside_nodes.size() > 0){
            //std::cout << outside_nodes.size() << " nodes are outside." << std::endl;
        }

        auto exited_nodes = get_exited_nodes(ps, prm.exit_plane, prm.Ln, prm.Lt);
        if (exited_nodes.size() > 0){
            //std::cout << exited_nodes.size() << " nodes have crossed the " << prm.exit_plane << " plane." << std::endl;
            split_random_nodes(exited_nodes, ps, gen, prm);
        }

        t += dt;
        ++it;
    }

    std::clock_t clock_1 = std::clock();
    Real duration = (clock_1-clock_0) / (Real) CLOCKS_PER_SEC;
    std::cout << "Total simulation time: " << duration << " seconds" << std::endl;

    return 0;
}
