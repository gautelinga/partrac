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

#include "StructuredInterpol.hpp"
#include "TriangleInterpol.hpp"
#include "experimental/integrator_RK.hpp"
#include "experimental/particles.hpp"
#include "experimental/initializer.hpp"
#include "Parameters.hpp"
// #include "Integrator.hpp"


/*
class Integrator_Omp : public Integrator {
public:
  Integrator_Omp(const int int_order);
  ~Integrator_Omp() {};
  template<typename InterpolType, typename T>
  std::set<Uint> step(InterpolType&, T&, double t, double s);
protected:
  int      m_int_order;
};

template<typename InterpolType, typename T>
std::set<Uint> Integrator_Omp::step(InterpolType& intp, T& ps, const double t, const double dt) {
    std::set<Uint> outside_nodes;

    for (Uint i=0; i < ps.N(); ++i){
        Vector3d x = ps.x(i);

        int cell_id = ps.get_cell_id(i);

        intp.probe(x, t, cell_id);

        Vector3d u_1 = intp.get_u();

        bool is_inside = false;
        if (uabs_est > m_u_min && ps.t_loc(i) < m_T){

            Vector3d dx = u_1 * dt;

            // Second-order terms
            if (m_int_order >= 2){
                dx += 0.5 * (intp.get_a() + intp.get_Ju()) * dt * dt;
            }
        }
        // count things
        if (is_inside){
            ++n_accepted;
            ps.set_x(i, x + dx);
            ps.set_t_loc(i, ps.t_loc(i) + dt);
            ps.set_cell_id(i, cell_id);
        }
        else {
            outside_nodes.insert(i);
            ++n_declined;
        }
    }
    return outside_nodes;
}
*/

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

int main(int argc, char* argv[])
{

    {
        std::cout << "======================================================================\n"
                  << "   Initialized omp_test ...                                           \n"
                  << "======================================================================" << std::endl;
    }    

    // Input parameters
    if (argc < 2) {
        std::cout << "Please specify an input file." << std::endl;
        return 0;
    }

    Parameters prm(argc, argv);

    double dt = prm.dt;
    Uint it = 0;
    double t = prm.t0;
    double T = prm.T;

    if (prm.num_threads > 0){
        omp_set_dynamic(0);
        omp_set_num_threads(prm.num_threads);
    }

    std::string infilename = std::string(argv[1]);
    TriangleInterpol intp(infilename);
    intp.set_U0(prm.U0);
    intp.set_int_order(prm.int_order);

    std::string folder = intp.get_folder();
    std::string rwfolder = folder + "/OMPTest/";
    create_folder(rwfolder);
    std::string newfolder = get_newfoldername(rwfolder, prm);
    create_folder(newfolder);
    prm.dump(newfolder, t);

    std::cout << "Initializing ParticleSet..." << std::endl;
    Particles<Particle> ps(prm.Nrw_max);

    std::random_device rd;
    std::vector<std::mt19937> gens;
    for (int i=0, N=omp_get_max_threads(); i<N; ++i) {
        gens.emplace_back(std::mt19937(rd()));
    }

    auto key = split_string(prm.init_mode, "_");

    RandomPointsInitializer init_state(key, prm, gens[0]);
    init_state.probe(intp);
    init_state.initialize(ps);

    prm.print();

    #pragma omp parallel
    {
        printf("Hello from process: %d\n", omp_get_thread_num());
    }

    // This part is unique
    //std::cout << "initializing Integrator..." << std::endl;
    //Integrator_RK4 integrator;
    //Integrator_Explicit integrator(prm.Dm, prm.int_order, gen);

    std::string h5fname = newfolder + "/data_from_t" + std::to_string(t) + ".h5";
    H5::H5File h5f(h5fname.c_str(), H5F_ACC_TRUNC);

    Uint int_stat_intv = int(prm.stat_intv/dt);
    Uint int_dump_intv = int(prm.dump_intv/dt);
    Uint int_checkpoint_intv = int(prm.checkpoint_intv/dt);
    Uint int_chunk_intv = int_dump_intv*prm.dump_chunk_size;

    std::map<std::string, bool> output_fields;
    output_fields["u"] = true; // !prm.minimal_output;
    output_fields["c"] = !prm.minimal_output;
    output_fields["p"] = true; // !prm.minimal_output && prm.output_all_props;
    output_fields["rho"] = false;  // !prm.minimal_output && prm.output_all_props;        
    output_fields["H"] = false;  //& !prm.minimal_output && ps.dim() > 0;
    output_fields["n"] = false; // !prm.minimal_output && ps.dim() > 1;
    output_fields["w"] = true;

    intp.update(t);
    intp.assign_fields(ps, output_fields);

    double sqrt2Dmdt = sqrt(2 * prm.Dm * dt);

    std::normal_distribution<double> rnd_normal(0.0, 1.0);

    while (t <= T){
        intp.update(t);
       
        // Statistics
        if (it % int_stat_intv == 0){
            std::cout << "Time = " << t << std::endl;
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

        auto clock_0 = std::chrono::high_resolution_clock::now();
        #pragma omp parallel for
        for ( auto & particle : ps.particles() ){
            Vector3d x = particle.x();
            int cell_id = particle.cell_id();
            PointValues ptvals(prm.U0);
            bool is_inside = intp.probe_light(x, t, cell_id);
            intp.probe_heavy(x, t, cell_id, ptvals);
            
            Vector3d dx = ptvals.get_u() * dt;
            if (prm.int_order > 1) dx += 0.5*(ptvals.get_Ju() + ptvals.get_a()) * dt * dt;
            if (prm.Dm > 0) {
                std::mt19937& gen = gens[omp_get_thread_num()];
                Vector eta = {rnd_normal(gen), rnd_normal(gen), rnd_normal(gen)};
                dx += sqrt2Dmdt * eta;
            }

            is_inside = intp.probe_light(x+dx, t+dt, cell_id);
            if (is_inside) {
                particle.x() = x + dx;
                particle.cell_id() = cell_id;
            }
        }
        auto clock_1 = std::chrono::high_resolution_clock::now();
        auto duration_1_0 = std::chrono::duration_cast<std::chrono::microseconds>(clock_1-clock_0);
        std::cout << "time elapsed = " << duration_1_0.count() << std::endl;

        // auto outside_nodes = integrator.step(intp, ps, t, dt);

        t += dt;
        ++it;
    }

    return EXIT_SUCCESS;
}
