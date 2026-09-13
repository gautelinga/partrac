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

#include "experimental/particles.hpp"
#include "utils.hpp"
#include "Params.hpp"
#include "run_folders.hpp"
#include "StructuredInterpol.hpp"
#include "TriangleInterpol.hpp"
#include "XDMFTriangleInterpol.hpp"

#include "experimental/integrator_RK.hpp"
#include "experimental/initializer.hpp"
#include "stats_columns.hpp"
//#include "experimental/statistics.hpp"

template<typename T>
std::vector<StatsColumn> phase_stats_columns( const Real t
                                             , T& ps
                                             , const unsigned long int n_declined
                                             )
{
  std::vector<StatsColumn> cols;
  Vector x_mean = {0., 0., 0.};
  Vector x_var = {0., 0., 0.};

  Vector u_mean = {0., 0., 0.};
  Vector u_var = {0., 0., 0.};

  Uint Nrw = ps.particles().size();

  Real w_mean = 0.;
  Real w_var = 0.;

  Real S_mean = 0.;
  Real rho_mean = 0.;
  
  // Phase quantities 1,2
  Uint Nrw1 = 0;

  Vector u1_mean = {0., 0., 0.};
  Vector u2_mean = {0., 0., 0.};
  
  Vector u1_var = {0., 0., 0.};
  Vector u2_var = {0., 0., 0.};

  Real w1_mean = 0.;
  Real w2_mean = 0.;
  
  Real w1_var = 0.;
  Real w2_var = 0.;

  Real S1_mean = 0.;
  Real S2_mean = 0.;

  #pragma omp parallel
  {
    Vector x_mean_loc = {0., 0., 0.};
    Vector u_mean_loc = {0., 0., 0.};

    Real w_mean_loc = 0.;
    Real S_mean_loc = 0.;
    Real rho_mean_loc = 0.;

    // Phase quantities 1,2
    Real w1_mean_loc = 0.;
    Real w2_mean_loc = 0.;

    Uint Nrw1_loc = 0;
    Vector u1_mean_loc = {0., 0., 0.};
    Vector u2_mean_loc = {0., 0., 0.};

    Real S1_mean_loc = 0.;
    Real S2_mean_loc = 0.;

    #pragma omp for
    for ( auto & particle : ps.particles() )
    {
        //auto & particle = ps.particles()[i];
        // Sample mean
        x_mean_loc += particle.get_x(); // /Nrw;
        u_mean_loc += particle.get_u(); // /Nrw;
        rho_mean_loc += particle.get_rho();
        w_mean_loc += particle.get_w();
        S_mean_loc += particle.get_S();

        double phi = particle.get_rho(); // consider renaming
        if (phi > 0){
            ++Nrw1_loc;
            u1_mean_loc += particle.get_u();
            w1_mean_loc += particle.get_w();
            S1_mean_loc += particle.get_S();
        }
        else { // if (phi < 0){
            u2_mean_loc += particle.get_u();
            w2_mean_loc += particle.get_w();
            S2_mean_loc += particle.get_S();
        }
    }
    #pragma omp critical
    {
        Nrw1 += Nrw1_loc;

        x_mean += x_mean_loc;
        u_mean += u_mean_loc;
        rho_mean += rho_mean_loc;
        w_mean += w_mean_loc;
        S_mean += S_mean_loc;

        u1_mean += u1_mean_loc;
        w1_mean += w1_mean_loc;
        S1_mean += S1_mean_loc;
        
        u2_mean += u2_mean_loc;
        w2_mean += w2_mean_loc;
        S2_mean += S2_mean_loc;
    } 
  }
  x_mean /= Nrw;
  u_mean /= Nrw;
  rho_mean /= Nrw;
  w_mean /= Nrw;
  S_mean /= Nrw;

  u1_mean /= Nrw1;
  w1_mean /= Nrw1;
  S1_mean /= Nrw1;
  
  Uint Nrw2 = Nrw - Nrw1;
  u2_mean /= Nrw2;
  w2_mean /= Nrw2;
  S2_mean /= Nrw2;

  #pragma omp parallel
  {
    Vector x_var_loc = {0., 0., 0.};
    Vector u_var_loc = {0., 0., 0.};

    Real w_var_loc = 0.;
    
    Vector u1_var_loc = {0., 0., 0.};
    Real w1_var_loc = 0.;

    Vector u2_var_loc = {0., 0., 0.};
    Real w2_var_loc = 0.;

    #pragma omp for
    for ( auto & particle : ps.particles() )
    {
        // Sample variance
        Vector dx = particle.get_x()-x_mean;
        x_var_loc += dx.cwiseProduct(dx);

        Vector du = particle.get_u()-u_mean;
        u_var_loc += du.cwiseProduct(du);

        w_var_loc += pow(particle.get_w()-w_mean, 2);

        double phi = particle.get_rho(); // consider renaming
        if (phi > 0){
            Vector du1 = particle.get_u()-u1_mean;
            u1_var_loc += du1.cwiseProduct(du1);
            w1_var_loc += pow(particle.get_w()-w1_mean, 2);
        }
        else {
            Vector du2 = particle.get_u()-u2_mean;
            u2_var_loc += du2.cwiseProduct(du2);
            w2_var_loc += pow(particle.get_w()-w2_mean, 2);
        }
    }
    #pragma omp critical
    {
        x_var += x_var_loc;
        u_var += u_var_loc;

        w_var += w_var_loc;

        u1_var += u1_var_loc;
        w1_var += w1_var_loc;

        u2_var += u2_var_loc;
        w2_var += w2_var_loc;
    }
  }
  // Unbiased sample variance
  x_var /= (Nrw-1);
  u_var /= (Nrw-1);
  w_var /= (Nrw-1);

  u1_var /= (Nrw1-1);
  w1_var /= (Nrw1-1);

  u2_var /= (Nrw2-1);
  w2_var /= (Nrw2-1);
  
  cols.push_back({"t", t});
  cols.push_back({"x_mean", x_mean[0]});
  cols.push_back({"y_mean", x_mean[1]});
  cols.push_back({"z_mean", x_mean[2]});
  cols.push_back({"x_var", x_var[0]});
  cols.push_back({"y_var", x_var[1]});
  cols.push_back({"z_var", x_var[2]});
  cols.push_back({"ux_mean", u_mean[0]});
  cols.push_back({"uy_mean", u_mean[1]});
  cols.push_back({"uz_mean", u_mean[2]});
  cols.push_back({"ux_var", u_var[0]});
  cols.push_back({"uy_var", u_var[1]});
  cols.push_back({"uz_var", u_var[2]});
  cols.push_back({"w_mean", w_mean});
  cols.push_back({"w_var", w_var});
  cols.push_back({"S_mean", S_mean});
  cols.push_back({"rho_mean", rho_mean});
  cols.push_back({"Nrw", double(Nrw), true});
  cols.push_back({"n_declined", double(n_declined), true});
  cols.push_back({"Nrw1", double(Nrw1), true});
  cols.push_back({"u1x_mean", u1_mean[0]});
  cols.push_back({"u1y_mean", u1_mean[1]});
  cols.push_back({"u1z_mean", u1_mean[2]});
  cols.push_back({"u1x_var", u1_var[0]});
  cols.push_back({"u1y_var", u1_var[1]});
  cols.push_back({"u1z_var", u1_var[2]});
  cols.push_back({"w1_mean", w1_mean});
  cols.push_back({"w1_var", w1_var});
  cols.push_back({"S1_mean", S1_mean});
  cols.push_back({"Nrw2", double(Nrw2), true});
  cols.push_back({"u2x_mean", u2_mean[0]});
  cols.push_back({"u2y_mean", u2_mean[1]});
  cols.push_back({"u2z_mean", u2_mean[2]});
  cols.push_back({"u2x_var", u2_var[0]});
  cols.push_back({"u2y_var", u2_var[1]});
  cols.push_back({"u2z_var", u2_var[2]});
  cols.push_back({"w2_mean", w2_mean});
  cols.push_back({"w2_var", w2_var});
  cols.push_back({"S2_mean", S2_mean});
  return cols;
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
            outside = !intp.locate(x0 + Dx);
        }
        node.x() += Dx;
    }
}

template<typename ParticleType>
void align_all(const std::vector<std::string>& key, Particles<ParticleType>& ps, std::mt19937 &gen){
    std::normal_distribution<Real> rnd_normal(0.0, 1.0);
    double delta = 1e-8;
    for ( auto & particle : ps.particles() ){
        Vector Dn = particle.get_u();
        if (contains(key[1], "x")){
            Dn[0] += delta * rnd_normal(gen);
        }
        if (contains(key[1], "y")){
            Dn[1] += delta * rnd_normal(gen);
        }
        if (contains(key[1], "z")){
            Dn[2] += delta * rnd_normal(gen);
        }
        Dn /= Dn.norm();
        particle.n() = Dn;
    }
}

#include "tracervectors_triangleRK4_schema.hpp"

int main(int argc, char* argv[])
{

    {
        std::cout << "======================================================================\n"
                  << "||  Initialized experimental tracer vectors.                        ||\n"
                  << "======================================================================" << std::endl;
    }
   // mpi.barrier();
    
    // Input parameters
    if (argc < 2){
        std::cout << "Please specify an input file." << std::endl;
        return 0;
    }

    partrac::Params prm = partrac::parse_or_exit(tracervectors_triangleRK4_schema(), argc, argv);

    if (prm.get<int>("num_threads") > 0){
        omp_set_dynamic(0);
        omp_set_num_threads(prm.get<int>("num_threads"));
    }

    std::string infilename = prm.input_file();

    std::cout << "Initializing TriangleInterpol." << std::endl;
    //TriangleInterpol intp(infilename);
    XDMFTriangleInterpol intp(infilename);
    // std::cout << "Initialized TriangleInterpol." << std::endl;

    intp.set_U0(prm.get<double>("U"));
    intp.set_int_order(2);  // To evaluate gradients

    std::string folder = intp.get_folder();
    RunFolders out = make_run_folders(folder, "TracerVectors", prm, NoRunIndex);
    const std::string& newfolder = out.run;

        if (prm.get<bool>("verbose")) prm.print();

    std::mt19937 gen;
    if (prm.get<bool>("random")) {
        std::random_device rd;
        gen.seed(rd());
    }
    else {
        std::seed_seq rd{prm.get<int>("seed")}; // + 0};
        gen.seed(rd);
    }

    Real dt = prm.get<double>("dt");
    Real t0 = std::max(intp.get_t_min(), prm.get<double>("t0"));
    Real T = std::min(intp.get_t_max(), prm.get<double>("T"));
    prm.set<double>("t0", t0);
    prm.set<double>("T", T);

    // This part is unique
    std::cout << "Initializing Integrator..." << std::endl;
    Integrator_RK4 integrator;

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
      write_stats_header(statfile, phase_stats_columns(0., ps, 0));
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
    output_fields["rho"] = true; // !prm.minimal_output && prm.output_all_props;        
    output_fields["H"] = !prm.get<bool>("minimal_output") && ps.dim() > 0;
    output_fields["n"] = true;
    output_fields["w"] = true;
    output_fields["S"] = true;

    // cell_type (for debugging)
    output_fields["cell_type"] = true;

    intp.update(t);
    intp.assign_fields(ps, output_fields);

    spin_all(key, ps, gen);
    //align_all(key, ps, gen);

    // Simulation start
    std::clock_t clock_0 = std::clock();

    double duration_step = 0.;
    double duration_other = 0.;

    const int sort_every = prm.get<int>("sort_every");

    while (t < T + dt/2){

        if (sort_every > 0 && it % sort_every == 0 && it > 0)

            ps.sort_by_cell();
        auto ct0 = std::chrono::high_resolution_clock::now();

        intp.update(t);

        // Update fields for output
        if (at_interval(it, prm.get<double>("dump_intv"), dt) || at_interval(it, prm.get<double>("stat_intv"), dt)){
            intp.assign_fields(ps, output_fields);
        }

        // Statistics
        if (at_interval(it, prm.get<double>("stat_intv"), dt)){
            std::cout << "Time = " << t << std::endl;
            std::cout << "(step: " << duration_step << ", other stuff: " << duration_other << ")" << std::endl;

            duration_step = 0;
            duration_other = 0;

            write_stats_row(statfile, phase_stats_columns(t, ps, integrator.get_declined()));

            intp.print_found();
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
        auto ct1 = std::chrono::high_resolution_clock::now();

        auto dct_other = std::chrono::duration_cast<std::chrono::microseconds>(ct1-ct0);
        duration_other += dct_other.count();

        auto outside_nodes = integrator.step_vec(intp, ps, t, dt);
        auto ct2 = std::chrono::high_resolution_clock::now();
        auto dct_step = std::chrono::duration_cast<std::chrono::microseconds>(ct2-ct1);
        duration_step += dct_step.count();
        // std::cout << dct10.count() << std::endl;

        if (outside_nodes.size() > 0){
            std::cout << outside_nodes.size() << " nodes are outside." << std::endl;
            // reinject_nodes(outside_nodes, key, ps, intp, gen);
            for ( auto & node_id : outside_nodes )
            {
                ps.particles()[node_id].c() = 2.0;
            }
        }

        t += dt;
        ++it;
    }

    std::clock_t clock_1 = std::clock();
    Real duration = (clock_1-clock_0) / (Real) CLOCKS_PER_SEC;
    std::cout << "Total simulation time: " << duration << " seconds" << std::endl;

    return 0;
}
