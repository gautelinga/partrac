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
//#include "H5Cpp.h"
//#include "hdf5.h"
#include <ctime>

#include "Error.hpp"
#include "io.hpp"
#include "rng.hpp"
#include "PointValues.hpp"
#include "Params.hpp"

#include "ParticleSet.hpp"
#include "Topology.hpp"
#include "Integrator.hpp"
#include "ExplicitIntegrator.hpp"
#include "RKIntegrator.hpp"
#include "Initializer.hpp"
#include "interpol_factory.hpp"
#include "morton.hpp"
#include "h5part.hpp"
#include "run_folders.hpp"

#include "interpol_schema.hpp"

inline void test_interpolation(Uint num_points, std::shared_ptr<Interpol> intp,
                        const std::string &newfolder, const double t0,
                        std::vector<std::mt19937>& gens){

  Vector3d x_min = intp->get_x_min();
  Vector3d x_max = intp->get_x_max();

  double Lx = intp->get_Lx();
  double Ly = intp->get_Ly();
  double Lz = intp->get_Lz();

  std::cout << "Lx=" << Lx << ", Ly=" << Ly << ", Lz=" << Lz << std::endl;

  intp->update(t0);

  // std::ofstream ofile(newfolder + "/interpolation.txt");

  std::vector<std::vector<double>> ptdata_threads_;


  std::vector<std::string> ptheader = {
    "x", "y", "z",
    "ux", "uy", "uz",
    "rho", "p", "divu", "vortz",
    "uxx", "uxy", "uxz",
    "uyx", "uyy", "uyz",
    "uzx", "uzy", "uzz"
  };

  // The points first, one generator a thread as before, so the set is the same
  std::vector<double> xs(3 * std::size_t(num_points));
  #pragma omp parallel
  {
    std::mt19937 &gen = gens[omp_get_thread_num()];

    std::uniform_real_distribution<> uni_dist_x(x_min[0], x_max[0]);
    std::uniform_real_distribution<> uni_dist_y(x_min[1], x_max[1]);
    std::uniform_real_distribution<> uni_dist_z(x_min[2], x_max[2]);

    #pragma omp for
    for (Uint i = 0; i < num_points; ++i){
      // one expression, so the three draws keep the order they had
      Vector3d x(uni_dist_x(gen), uni_dist_y(gen), uni_dist_z(gen));
      xs[3*std::size_t(i)]     = x[0];
      xs[3*std::size_t(i) + 1] = x[1];
      xs[3*std::size_t(i) + 2] = x[2];
    }
  }

  // Probed in Morton order, a contiguous range a thread: consecutive queries
  // descend the same subtree and gather the same dofs. The rows come out in
  // that order, which no consumer of this file depends on.
  const int qdim = (x_max[2] > x_min[2]) ? 3 : 2;
  const partrac::MortonBox box(x_min, x_max, qdim);
  const std::vector<std::uint32_t> order =
    partrac::morton_order(xs.data(), std::size_t(num_points), 3, box);

  #pragma omp parallel
  {
    #pragma omp single
    ptdata_threads_.resize(omp_get_num_threads());

    auto& ptdata_loc_ = ptdata_threads_[omp_get_thread_num()];
    ptdata_loc_.reserve(num_points * ptheader.size() / omp_get_num_threads());

    #pragma omp for schedule(static)
    for (Uint q = 0; q < num_points; ++q){
      CellPos pos;
      const double* p = xs.data() + 3*std::size_t(order[q]);
      Vector3d x(p[0], p[1], p[2]);

      bool inside = intp->locate(x, t0, pos);
      if (inside){
        PointValues ptvals(intp->get_U0());
        intp->evaluate(x, t0, pos, ptvals);

        Vector3d u = ptvals.get_u();
        Matrix3d gradu = ptvals.get_J();

        ptdata_loc_.insert(ptdata_loc_.end(), {
          x[0], x[1], x[2],
          u[0], u[1], u[2],
          ptvals.get_rho(), ptvals.get_p(),
          gradu(0,0) + gradu(1,1) + gradu(2,2), 
          gradu(1,0) - gradu(0,1),
          gradu(0,0), gradu(0,1), gradu(0,2),
          gradu(1,0), gradu(1,1), gradu(1,2),
          gradu(2,0), gradu(2,1), gradu(2,2)
        });
      }
    }
  }
  std::cout << "Done probing." << std::endl;

  Uint n_inside = 0;
  for (auto& v : ptdata_threads_)
    n_inside += v.size() / ptheader.size();

  std::cout << "Inside:             " << n_inside << "/" << num_points << std::endl;
  std::cout << "Approximate volume: " << (n_inside*Lx*Ly*Lz)/num_points << std::endl;
  std::cout << "Approximate area:   " << (n_inside*Lx*Ly)/num_points << std::endl;

  std::vector<double> ptdata_;
  ptdata_.reserve(n_inside * ptheader.size());

  for (auto& v : ptdata_threads_)
    ptdata_.insert(ptdata_.end(), v.begin(), v.end());

  write_h5part(newfolder + "/interpolation.h5part", ptheader, ptdata_);

  std::cout << "Done writing." << std::endl;

  // ofile.close();
}

static int run(int argc, char* argv[])
{

    std::cout << "Initialized Interpolator." << std::endl;

  // Input parameters
  if (argc < 2) {
    std::cout << "Specify an input file." << std::endl;
    return 1;
  }
  partrac::Params prm = partrac::parse_or_exit(interpol_schema(), argc, argv);

  if (prm.get<int>("num_threads") > 0){
      omp_set_dynamic(0);
      omp_set_num_threads(prm.get<int>("num_threads"));
  }

  std::string infilename = prm.input_file();

  std::shared_ptr<Interpol> intp;
  set_interpolate_mode(intp, prm.get<std::string>("mode"), infilename);
  intp->set_U0(prm.get<double>("U"));
  intp->set_int_order(prm.get<int>("int_order"));

  std::string folder = intp->get_folder();
  const std::string newfolder = make_run_folders(folder, "Interpolation", prm,
                                                NoSubfolders | NoRunIndex).run;

  // Parallel generators
  std::vector<std::mt19937> gens = make_generators(prm);

  double t0 = std::max(intp->get_t_min(), prm.get<double>("t0"));
  std::cout << "Testing interpolation..." << std::endl;
  test_interpolation(prm.get<Uint>("Nrw"), intp, newfolder, t0, gens);

  return 0;
}

int main(int argc, char* argv[])
{
  return partrac::report_errors([&]{ return run(argc, argv); });
}
