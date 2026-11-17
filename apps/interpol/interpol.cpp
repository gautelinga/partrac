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
#include <ctime>
#include <limits>

#include "Error.hpp"
#include "wait_policy.hpp"
#include "io.hpp"
#include "rng.hpp"
#include "PointValues.hpp"
#include "Params.hpp"

#include "ParticleSet.hpp"
#include "Topology.hpp"
#include "Integrator.hpp"
#include "ExplicitIntegrator.hpp"
#include "RKIntegrator.hpp"
#include "interpol_factory.hpp"
#include "morton.hpp"
#include "h5part.hpp"
#include "h5direct.hpp"
#include "geometry.hpp"
#include "run_folders.hpp"

#include "interpol_schema.hpp"

// The columns of a probe at x, in the header's order up to inside, into row;
// false and row untouched outside the fluid
inline bool probe(Interpol& intp, const Vector3d& x, const double t0, const bool phase,
                  const bool phase_gradient, double* row){
  CellPos pos;
  if (!intp.locate(x, t0, pos))
    return false;
  PointValues ptvals(intp.get_U0());
  intp.evaluate(x, t0, pos, ptvals);

  Vector3d u = ptvals.get_u();
  Matrix3d gradu = ptvals.get_J();

  const double cols[] = {
    x[0], x[1], x[2],
    u[0], u[1], u[2],
    ptvals.get_rho(), ptvals.get_p(),
    gradu(0,0) + gradu(1,1) + gradu(2,2),
    gradu(1,0) - gradu(0,1),
    gradu(0,0), gradu(0,1), gradu(0,2),
    gradu(1,0), gradu(1,1), gradu(1,2),
    gradu(2,0), gradu(2,1), gradu(2,2)
  };
  row = std::copy(std::begin(cols), std::end(cols), row);
  if (phase) *row++ = ptvals.get_phi();
  if (phase_gradient){
    Vector3d g;
    intp.evaluate_phase_gradient(x, t0, pos, g);
    row = std::copy(g.data(), g.data() + 3, row);
  }
  return true;
}

// Random points in the box, one generator a thread
inline std::vector<double> random_points(Uint num_points, const Interpol& intp,
                                         std::vector<std::mt19937>& gens){
  const Vector3d x_min = intp.get_x_min();
  const Vector3d x_max = intp.get_x_max();
  std::vector<double> xs(3 * std::size_t(num_points));
  #pragma omp parallel
  {
    std::mt19937 &gen = gens[omp_get_thread_num()];

    std::uniform_real_distribution<> uni_dist_x(x_min[0], x_max[0]);
    std::uniform_real_distribution<> uni_dist_y(x_min[1], x_max[1]);
    std::uniform_real_distribution<> uni_dist_z(x_min[2], x_max[2]);

    #pragma omp for
    for (Uint i = 0; i < num_points; ++i){
      // The three draws in one expression: the point set depends on their order
      Vector3d x(uni_dist_x(gen), uni_dist_y(gen), uni_dist_z(gen));
      xs[3*std::size_t(i)]     = x[0];
      xs[3*std::size_t(i) + 1] = x[1];
      xs[3*std::size_t(i) + 2] = x[2];
    }
  }
  return xs;
}

// Points to probe from a file: positions, and the ids beside them if any
struct GivenPoints {
  std::vector<double> xs;    // x, y, z a point
  std::vector<double> ids;   // one a point, or empty
};

// The N x 3 float positions of a dataset, in its order, and the integer ids of
// the dataset id beside it (a dump's <t>/id), one a row, when it is there
inline GivenPoints read_points(const std::string& path, const std::string& dataset){
  const partrac::H5Id h = partrac::h5_open_read(path);
  const partrac::H5DatasetInfo info = partrac::h5_dataset_info(h, dataset);
  if (info.rank() != 2 || info.cols() != 3){
    std::ostringstream shape;
    for (std::size_t k = 0; k < info.rank(); ++k) shape << (k ? " x " : "") << info.shape[k];
    partrac::fail(path, ":", dataset, " is ", info.rank() ? shape.str() : std::string("a scalar"),
                  ", not N x 3 positions");
  }
  if (info.type_class != H5T_FLOAT)
    partrac::fail(path, ":", dataset, " holds no floating-point positions");
  GivenPoints g;
  partrac::h5_read(h, dataset, g.xs);

  const std::size_t slash = dataset.rfind('/');
  const std::string idpath = (slash == std::string::npos ? std::string() : dataset.substr(0, slash + 1)) + "id";
  if (idpath != dataset && H5Lexists(h, idpath.c_str(), H5P_DEFAULT) > 0){
    const partrac::H5DatasetInfo ii = partrac::h5_dataset_info(h, idpath);
    if (ii.type_class != H5T_INTEGER || ii.rows() != info.rows() || ii.cols() != 1)
      partrac::fail(path, ":", idpath, " is not one integer id a point of ", dataset);
    partrac::h5_read(h, idpath, g.ids);   // exact as doubles up to 2^53
  }
  return g;
}

// The time a dataset path names in its first group, as a dump's <t>/points, or NaN
inline double time_in_path(const std::string& dataset){
  const std::size_t slash = dataset.find('/');
  if (slash == std::string::npos || slash == 0) return std::numeric_limits<double>::quiet_NaN();
  const std::string head = dataset.substr(0, slash);
  char* end = nullptr;
  const double t = std::strtod(head.c_str(), &end);
  return (end == head.c_str() + head.size()) ? t : std::numeric_limits<double>::quiet_NaN();
}

// Probes the fields at the points xs (x, y, z a point). Random points: the
// rows of those inside the fluid. Given points: a row a point in their order,
// the id beside it when given, inside 1 or 0; an outside row (a point outside
// the fluid, or not finite) is NaN but for x, y, z and the id.
inline void probe_points(std::shared_ptr<Interpol> intp, const std::string &newfolder,
                         const double t0, const std::vector<double>& xs, const bool given,
                         const std::vector<double>& ids = {}){

  const Uint num_points = xs.size() / 3;
  Vector3d x_min = intp->get_x_min();
  Vector3d x_max = intp->get_x_max();

  double Lx = intp->get_Lx();
  double Ly = intp->get_Ly();
  double Lz = intp->get_Lz();

  std::cout << "Lx=" << Lx << ", Ly=" << Ly << ", Lz=" << Lz << std::endl;

  intp->update(t0);

  std::vector<std::string> ptheader = {
    "x", "y", "z",
    "ux", "uy", "uz",
    "rho", "p", "divu", "vortz",
    "uxx", "uxy", "uxz",
    "uyx", "uyy", "uyz",
    "uzx", "uzy", "uzz"
  };
  const bool phase = intp->has_phase_field();
  if (phase) ptheader.push_back("phi");
  const bool phase_gradient = intp->has_phase_gradient();
  if (phase_gradient) ptheader.insert(ptheader.end(), {"dphi_dx", "dphi_dy", "dphi_dz"});
  const bool with_ids = given && !ids.empty();
  if (with_ids) ptheader.push_back("id");
  if (given) ptheader.push_back("inside");
  const std::size_t ncols = ptheader.size();

  // Probed in Morton order, a contiguous range a thread: consecutive queries
  // descend the same subtree and gather the same dofs. Random rows come out
  // in that order, which no consumer of this file depends on. Given points
  // are keyed as folded into the box (locate wraps them itself, as it does
  // the tracers' unwrapped positions); a point not finite keys as x_min.
  const int qdim = (x_max[2] > x_min[2]) ? 3 : 2;
  const partrac::MortonBox box(x_min, x_max, qdim);
  const std::vector<std::uint32_t> order =
      given ? partrac::morton_order(xs.data(), std::size_t(num_points), 3, box, x_min, x_max)
            : partrac::morton_order(xs.data(), std::size_t(num_points), 3, box);

  std::vector<double> ptdata_;
  Uint n_inside = 0;
  if (given){
    // Written in place, a row a point; filled where probed, first touch
    // by the thread that probes it
    ptdata_.resize(std::size_t(num_points) * ncols);
    #pragma omp parallel for schedule(static) reduction(+:n_inside)
    for (std::ptrdiff_t q = 0; q < std::ptrdiff_t(num_points); ++q){
      const std::size_t i = order[q];
      const double* p = xs.data() + 3*i;
      double* r = ptdata_.data() + i*ncols;
      std::fill(r, r + ncols, std::numeric_limits<double>::quiet_NaN());
      const Vector3d x(p[0], p[1], p[2]);
      const bool inside = x.allFinite() && probe(*intp, x, t0, phase, phase_gradient, r);
      n_inside += inside;
      r[0] = x[0]; r[1] = x[1]; r[2] = x[2];
      if (with_ids) r[ncols - 2] = ids[i];
      r[ncols - 1] = inside;
    }
    std::cout << "Done probing." << std::endl;
    std::cout << "Inside:             " << n_inside << "/" << num_points << std::endl;
  }
  else {
    std::vector<std::vector<double>> ptdata_threads_;
    #pragma omp parallel reduction(+:n_inside)
    {
      #pragma omp single
      ptdata_threads_.resize(omp_get_num_threads());

      auto& ptdata_loc_ = ptdata_threads_[omp_get_thread_num()];
      ptdata_loc_.reserve(num_points * ncols / omp_get_num_threads());
      std::vector<double> row(ncols);

      #pragma omp for schedule(static)
      for (Uint q = 0; q < num_points; ++q){
        const double* p = xs.data() + 3*std::size_t(order[q]);
        const Vector3d x(p[0], p[1], p[2]);
        if (probe(*intp, x, t0, phase, phase_gradient, row.data())){
          ++n_inside;
          ptdata_loc_.insert(ptdata_loc_.end(), row.begin(), row.end());
        }
      }
    }
    std::cout << "Done probing." << std::endl;
    std::cout << "Inside:             " << n_inside << "/" << num_points << std::endl;
    std::cout << "Approximate volume: " << (n_inside*Lx*Ly*Lz)/num_points << std::endl;
    std::cout << "Approximate area:   " << (n_inside*Lx*Ly)/num_points << std::endl;

    ptdata_.reserve(n_inside * ncols);
    for (auto& v : ptdata_threads_)
      ptdata_.insert(ptdata_.end(), v.begin(), v.end());
  }

  write_h5part(newfolder + "/interpolation.h5part", ptheader, ptdata_, t0);

  std::cout << "Done writing." << std::endl;
}

static int run(int argc, char* argv[])
{
  std::cout << "Initialized interpol." << std::endl;

  partrac::Params prm = partrac::parse_or_exit(interpol_schema(), argc, argv);
  if (prm.check_only())
    partrac::fail("interpol has no --check: it only probes the field");

  if (prm.get<int>("num_threads") > 0){
      omp_set_dynamic(0);
      omp_set_num_threads(prm.get<int>("num_threads"));
  }

  // Given points, before the field is loaded: they set Nrw, and name the
  // output folder, so that probes of other files or datasets do not land on it
  const std::string points = prm.get<std::string>("points");
  const std::string points_dataset = prm.get<std::string>("points_dataset");
  GivenPoints given;
  if (!points.empty()){
    given = read_points(points, points_dataset);
    const Uint n = given.xs.size() / 3;
    if (prm.was_set("Nrw") && prm.get<Uint>("Nrw") != n)
      partrac::fail("Nrw=", prm.get<Uint>("Nrw"), " but ", points, " has ", n, " points");
    prm.set<Uint>("Nrw", n);
    std::string name = std::filesystem::path(points).stem().string() + "-" + points_dataset;
    for (char& c : name)
      if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '.' || c == '-' || c == '_')) c = '-';
    prm.set<std::string>("tag", prm.get<std::string>("tag") + "_pts-" + name);
    std::cout << "Probing " << n << " points from " << points << ":" << points_dataset
              << (given.ids.empty() ? "" : ", with their ids") << std::endl;
  }

  std::string infilename = prm.input_file();

  std::shared_ptr<Interpol> intp;
  set_interpolate_mode(intp, prm.get<std::string>("mode"), infilename);
  intp->set_U0(prm.get<double>("U"));
  intp->set_int_order(prm.get<int>("int_order"));

  // Parallel generators; before the folders, which name the seed
  std::vector<std::mt19937> gens = make_generators(prm);

  std::vector<double> xs = points.empty() ? random_points(prm.get<Uint>("Nrw"), *intp, gens)
                                          : std::move(given.xs);

  std::string folder = intp->get_folder();
  if (!points.empty()){
    // A probe of given points is not redone over an earlier one
    const std::string out = get_newfoldername(folder + "/Interpolation/", prm) + "interpolation.h5part";
    if (std::filesystem::exists(out))
      partrac::fail(out, " exists: remove it, or give tag= to probe these points again");
  }
  const std::string newfolder = make_run_folders(folder, "Interpolation", prm,
                                                NoSubfolders | NoRunIndex).run;

  double t0 = std::max(intp->get_t_min(), prm.get<double>("t0"));
  const double t_points = time_in_path(points_dataset);
  if (!points.empty() && std::isfinite(t_points) && std::abs(t_points - t0) > 1e-9*std::max(1., std::abs(t0)))
    std::cout << "Note: " << points_dataset << " names t = " << t_points << ", the field is probed at t0 = "
              << t0 << " (set t0= for a field that changes in time)" << std::endl;
  std::cout << "Probing the fields..." << std::endl;
  probe_points(intp, newfolder, t0, xs, !points.empty(), given.ids);

  return 0;
}

int main(int argc, char* argv[])
{
  partrac::wait_passively(argv);
  return partrac::report_errors([&]{ return run(argc, argv); });
}
