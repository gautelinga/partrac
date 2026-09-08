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
#include "Params.hpp"
#include "rng.hpp"
// #include "Integrator.hpp"



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

#include "omp_test_schema.hpp"

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
        return 1;
    }

    partrac::Params prm = partrac::parse_or_exit(omp_test_schema(), argc, argv);

    double dt = prm.get<double>("dt");
    Uint it = 0;
    double t = prm.get<double>("t0");
    double T = prm.get<double>("T");

    if (prm.get<int>("num_threads") > 0){
        omp_set_dynamic(0);
        omp_set_num_threads(prm.get<int>("num_threads"));
    }

    std::string infilename = prm.input_file();
    TriangleInterpol intp(infilename);
    intp.set_U0(prm.get<double>("U"));
    intp.set_int_order(prm.get<int>("int_order"));

    std::string folder = intp.get_folder();
    std::string rwfolder = folder + "/OMPTest/";
    create_folder(rwfolder);
    std::string newfolder = get_newfoldername(rwfolder, prm);
    create_folder(newfolder);
    prm.dump(newfolder, t);

    std::cout << "Initializing ParticleSet..." << std::endl;
    Particles<Particle> ps(prm.get<Uint>("Nrw_max"));

    std::vector<std::mt19937> gens = make_generators(prm);

    auto key = split_string(prm.get<std::string>("init_mode"), "_");

    experimental::RandomPointsInitializer init_state(key, prm, gens[0]);
    init_state.probe(intp);
    init_state.initialize(ps);

    if (prm.get<bool>("verbose")) prm.print();

    #pragma omp parallel
    {
        printf("Hello from process: %d\n", omp_get_thread_num());
    }

    // This part is unique
    //std::cout << "initializing Integrator..." << std::endl;
    //Integrator_RK4 integrator;
    //Integrator_Explicit integrator(prm.Dm, prm.int_order, gen);

    std::string h5fname = newfolder + "/data_from_t" + std::to_string(t) + ".h5";
    // no dump file at all when dumping is off
    H5::H5File h5f;
    if (prm.get<double>("dump_intv") > 0.){
      { H5::H5File create(h5fname.c_str(), H5F_ACC_TRUNC); }
      h5f.openFile(h5fname.c_str(), H5F_ACC_RDWR);
    }

  const double chunk_intv = prm.get<double>("dump_intv")*prm.get<int>("dump_chunk_size");

    std::map<std::string, bool> output_fields;
    output_fields["u"] = true; // !prm.minimal_output;
    output_fields["c"] = !prm.get<bool>("minimal_output");
    output_fields["p"] = true; // !prm.minimal_output && prm.output_all_props;
    output_fields["rho"] = false;  // !prm.minimal_output && prm.output_all_props;        
    output_fields["H"] = false;  //& !prm.minimal_output && ps.dim() > 0;
    output_fields["n"] = false; // !prm.minimal_output && ps.dim() > 1;
    output_fields["w"] = true;

    intp.update(t);
    intp.assign_fields(ps, output_fields);

    double sqrt2Dmdt = sqrt(2 * prm.get<double>("Dm") * dt);

    std::normal_distribution<double> rnd_normal(0.0, 1.0);

    while (t <= T){
        intp.update(t);
       
        // Statistics
        if (at_interval(it, prm.get<double>("stat_intv"), dt)){
            std::cout << "Time = " << t << std::endl;
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

        auto clock_0 = std::chrono::high_resolution_clock::now();
        // the distribution caches a value between calls, so give each thread its own
        #pragma omp parallel for firstprivate(rnd_normal)
        for ( auto & particle : ps.particles() ){
            Vector3d x = particle.x();
            int cell_id = particle.cell_id();
            PointValues ptvals(prm.get<double>("U"));
            bool is_inside = intp.locate(x, t, cell_id);
            intp.evaluate(x, t, cell_id, ptvals);
            
            Vector3d dx = ptvals.get_u() * dt;
            if (prm.get<int>("int_order") > 1) dx += 0.5*(ptvals.get_Ju() + ptvals.get_a()) * dt * dt;
            if (prm.get<double>("Dm") > 0) {
                std::mt19937& gen = gens[omp_get_thread_num()];
                Vector eta = {rnd_normal(gen), rnd_normal(gen), rnd_normal(gen)};
                dx += sqrt2Dmdt * eta;
            }

            is_inside = intp.locate(x+dx, t+dt, cell_id);
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
