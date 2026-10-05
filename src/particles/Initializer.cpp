#include <array>
#include <cmath>
#include <iostream>
#include <functional>
#include <numeric>
#include <random>
#include "Error.hpp"
#include "Initializer.hpp"
#include "mesh.hpp"
#include "files.hpp"
#include "geometry.hpp"

namespace {

// Lexicographic order of positions
struct less_than_op {
  bool operator() (const Vector3d &a, const Vector3d &b) const {
    return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]) || (a[0] == b[0] && a[1] == b[1] && a[2] < b[2]);
  }
};

// (x0, y0, z0)
Vector3d initial_position(const partrac::Params& prm){
  return {prm.get<double>("x0"), prm.get<double>("y0"), prm.get<double>("z0")};
}

// 0, 1 or 2 for x, y or z; -1 for anything else
int axis_of(const std::string& dir){
  if (dir == "x") return 0;
  if (dir == "y") return 1;
  if (dir == "z") return 2;
  return -1;
}

// Unit vector along x, y or z; zero for anything else
Vector3d unit_vector(const std::string& dir){
  Vector3d n(0., 0., 0.);
  const int axis = axis_of(dir);
  if (axis >= 0)
    n[axis] = 1.;
  return n;
}

// Which of x, y and z the directions name
std::array<bool, 3> axis_mask(const std::string& dirs){
  return {contains(dirs, "x"), contains(dirs, "y"), contains(dirs, "z")};
}

// The points inside, in order; an edge joins two that follow each other
InitialState polyline_inside(const std::vector<Vector3d>& points, const std::vector<bool>& inside){
  InitialState state;
  for (Uint i=0; i < points.size(); ++i){
    if (!inside[i])
      continue;
    state.nodes.push_back(points[i]);
    const Uint inode = state.nodes.size()-1;
    if (i > 0 && inside[i-1])
      state.edges.push_back({{inode-1, inode}, dist(state.nodes[inode-1], points[i])});
  }
  return state;
}

// Draws in a row that miss the domain before a mode gives up
const Uint max_misses = 1000000;

// Nrw draws inside the domain; fails after max_misses misses in a row
std::vector<Vector3d> sample_inside(const std::function<Vector3d()>& draw, std::shared_ptr<Interpol> intp, const partrac::Params& prm){
  const Uint Nrw = prm.get<Uint>("Nrw");
  std::vector<Vector3d> points;
  Uint misses = 0;
  while (points.size() < Nrw){
    const Vector3d x = draw();
    if (intp->locate(x)){
      points.push_back(x);
      misses = 0;
    }
    else if (++misses == max_misses){
      partrac::fail("init_mode ", prm.get<std::string>("init_mode"), ": no points inside the domain in ",
                    max_misses, " draws in a row, with ", points.size(), " of ", Nrw, " placed");
    }
  }
  return points;
}

// Triangulated surface being remeshed; no inlet
struct Surface : Connectivity {
  Surface(std::shared_ptr<Interpol> intp, const std::vector<Vector3d>& nodes, const Uint Nrw_max) : ps(intp, Nrw_max) {
    ps.add(nodes, 0);
  }
  ParticleSet ps;
};

void compute_maps(Surface& s){
  compute_edge2faces(s.edge2faces, s.faces, s.edges);
  compute_node2edges(s.node2edges, s.edges, s.ps.N());
}

// Split edges longer than ds_max
Uint refine_surface(Surface& s, const double ds_max, const bool check_if_inside){
  return sheet_refinement(s, s.ps, ds_max, 0.0, StuckEdge::Keep, check_if_inside);
}

// Collapse edges shorter than ds_min
Uint coarsen_surface(Surface& s, const double ds_min){
  return sheet_coarsening(s, s.ps, ds_min, 0.0);
}

// Edge lengths and face areas from the positions
void settle_lengths(Surface& s){
  for (auto & edge : s.edges)
    edge.second = s.ps.dist(edge.first[0], edge.first[1]);
  for (auto & face : s.faces)
    face.second = s.ps.triangle_area(face.first[0], face.first[1], s.edges);
}

InitialState surface_state(const Surface& s){
  InitialState state;
  for (Uint irw=0; irw < s.ps.N(); ++irw)
    state.nodes.push_back(s.ps.x(irw));
  state.edges = s.edges;
  state.faces = s.faces;
  return state;
}

// Cell centres and init_weight weights on a grid over the sampled axes
struct WeightedGrid {
  std::vector<Vector3d> centres;
  std::vector<double> weights;
  Vector3d h;                          // cell size
};

WeightedGrid weighted_grid(const std::array<bool, 3>& sampled, std::shared_ptr<Interpol> intp, const partrac::Params& prm){
  const Vector3d x0 = initial_position(prm);
  const Vector3d x_min = intp->get_x_min();
  const Vector3d L = intp->get_x_max() - x_min;
  // About N_est cells
  const Uint N_est = 1000000;
  double extent = 1.;
  Uint n_sampled = 0;
  for (Uint d=0; d < 3; ++d){
    if (sampled[d]){
      extent *= L[d];
      ++n_sampled;
    }
  }
  const double dx_est = n_sampled == 1 ? extent/N_est : pow(extent/N_est, 1./n_sampled);
  std::array<Uint, 3> N = {1, 1, 1};
  WeightedGrid grid;
  for (Uint d=0; d < 3; ++d){
    if (sampled[d])
      N[d] = L[d]/dx_est+1;
    grid.h[d] = L[d]/N[d];
  }

  const std::string init_weight = prm.get<std::string>("init_weight");
  const Uint n_cells = N[0]*N[1]*N[2];
  grid.weights.resize(n_cells);
  grid.centres.resize(n_cells);
  #pragma omp parallel for
  for (Uint n=0; n < n_cells; ++n){
    const Uint idx[3] = {n/(N[1]*N[2]), (n/N[2]) % N[1], n % N[2]};
    Vector3d x = x0;
    for (Uint d=0; d < 3; ++d)
      if (sampled[d])
        x[d] = x_min[d]+(idx[d]+0.5)*grid.h[d];
    PointValues ptvals(intp->get_U0());
    // Uniform weights need no field
    if (init_weight != "uniform")
      intp->evaluate(x, ptvals);
    double w = 1.;
    if (init_weight == "ux") w = std::abs(ptvals.U[0]);
    else if (init_weight == "uy") w = std::abs(ptvals.U[1]);
    else if (init_weight == "uz") w = std::abs(ptvals.U[2]);
    else if (init_weight == "u") w = ptvals.U.norm();
    grid.weights[n] = w;
    grid.centres[n] = x;
  }
  return grid;
}

}  // namespace

InitialState init_point(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen){
  const Vector3d x0 = initial_position(prm);
  if (!intp->locate(x0))
    partrac::fail("init_mode ", prm.get<std::string>("init_mode"), ": (x0, y0, z0) is not inside the domain");
  InitialState state;
  state.nodes.assign(prm.get<Uint>("Nrw"), x0);
  return state;
}

InitialState init_uniform(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen){
  const Uint Nrw = prm.get<Uint>("Nrw");
  const int axis = axis_of(key[1]);
  if (axis < 0)
    partrac::fail("init_mode ", prm.get<std::string>("init_mode"), ": unknown direction ", key[1]);
  if (intp->get_x_max()[axis] - intp->get_x_min()[axis] <= 1e-12)
    partrac::fail("init_mode ", prm.get<std::string>("init_mode"), ": the domain has no extent along ", key[1]);
  // Ends on the domain's walls
  Vector3d x_a = initial_position(prm);
  Vector3d x_b = x_a;
  x_a[axis] = intp->get_x_min()[axis];
  x_b[axis] = intp->get_x_max()[axis];
  const Vector3d Dx = (x_b - x_a) / (Nrw-1);
  std::vector<Vector3d> points;
  std::vector<bool> inside;
  for (Uint irw=0; irw < Nrw; ++irw){
    points.push_back(x_a + Dx * irw);
    inside.push_back(intp->locate(points.back()));
  }
  return polyline_inside(points, inside);
}

InitialState init_strip(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen){
  const Uint Nrw = prm.get<Uint>("Nrw");
  const double La = prm.get<double>("La");
  const Vector3d x0 = initial_position(prm);
  const Vector3d n = unit_vector(key[1]);
  // Ends
  const Vector3d x00 = x0 - La/2*n;
  const Vector3d x01 = x0 + La/2*n;
  std::vector<Vector3d> points;
  std::vector<bool> inside;
  for (Uint i=0; i < Nrw; ++i){
    const double alpha = double(i)/(Nrw-1);
    points.push_back(alpha * x00 + (1.-alpha) * x01);
    inside.push_back(intp->locate(points.back()));
  }
  InitialState state = polyline_inside(points, inside);
  if (state.nodes.empty())
    partrac::fail("init_mode ", prm.get<std::string>("init_mode"), ": no point of the strip inside the domain");
  return state;
}

InitialState init_sheet(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen){
  const double La = prm.get<double>("La");
  const double Lb = prm.get<double>("Lb");
  const Vector3d x0 = initial_position(prm);
  // In-plane axes
  Vector3d ta(0., 0., 0.);
  Vector3d tb(0., 0., 0.);
  if (key[1] == "xy"){ ta[0] = 1.; tb[1] = 1.; }
  if (key[1] == "xz"){ ta[0] = 1.; tb[2] = 1.; }
  if (key[1] == "yz"){ ta[1] = 1.; tb[2] = 1.; }

  // Corners, two triangles
  const std::vector<Vector3d> corners = {x0 - La/2*ta - Lb/2*tb,
                                         x0 + La/2*ta - Lb/2*tb,
                                         x0 + La/2*ta + Lb/2*tb,
                                         x0 - La/2*ta + Lb/2*tb};
  Surface s(intp, corners, prm.get<Uint>("Nrw_max"));
  for (const auto & e : std::vector<std::array<Uint, 2>>{{0, 1}, {0, 2}, {1, 2}, {2, 3}, {3, 0}})
    s.edges.push_back({e, dist(corners[e[0]], corners[e[1]])});
  s.faces.push_back({{0, 2, 1}, La*Lb/2});
  s.faces.push_back({{1, 3, 4}, La*Lb/2});
  compute_maps(s);

  // Refine to ds_init
  while (refine_surface(s, prm.get<double>("ds_init"), false) > 0);
  settle_lengths(s);
  compute_maps(s);

  // Drop edges and faces at nodes in the solid
  std::vector<bool> edge_isactive(s.edges.size(), true);
  std::vector<bool> face_isactive(s.faces.size(), true);
  std::vector<bool> node_isactive(s.ps.N(), true);
  for (Uint irw=0; irw < s.ps.N(); ++irw){
    int cell_id = s.ps.get_cell_id(irw);
    if (intp->locate(s.ps.x(irw), 0., cell_id))
      continue;
    for (const auto & iedge : s.node2edges[irw]){
      edge_isactive[iedge] = false;
      for (const auto & iface : s.edge2faces[iedge])
        face_isactive[iface] = false;
    }
  }
  remove_inactive(s, face_isactive, edge_isactive, node_isactive, s.ps);

  // Coarsen to ds_min, refine to ds_max
  Uint n_rem, n_add;
  int attempt = 0;
  do {
    n_rem = coarsen_surface(s, prm.get<double>("ds_min"));
    n_add = refine_surface(s, prm.get<double>("ds_max"), true);
    ++attempt;
  } while ((n_add > 0 || n_rem > 0) && attempt < 100);
  return surface_state(s);
}

InitialState init_ellipsoid(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen){
  const double La = prm.get<double>("La");
  const double Lb = prm.get<double>("Lb");
  const Vector3d x_c = initial_position(prm);
  // Squared semi-axes
  double lx2 = Lb*Lb;
  double ly2 = Lb*Lb;
  double lz2 = Lb*Lb;
  if (key[1] == "xy") lz2 = La*La;
  if (key[1] == "xz") ly2 = La*La;
  if (key[1] == "yz") lx2 = La*La;

  // Tetrahedron
  const double R = sqrt(La*Lb);
  const std::vector<Vector3d> corners = {
    {x_c[0] - R/sqrt(2.), x_c[1] - R/sqrt(6.0),   x_c[2] - R/sqrt(3.0)/2},
    {x_c[0] + R/sqrt(2.), x_c[1] - R/sqrt(6.0),   x_c[2] - R/sqrt(3.0)/2},
    {x_c[0],              x_c[1] + R*sqrt(2./3.), x_c[2] - R/sqrt(3.0)/2},
    {x_c[0],              x_c[1],                 x_c[2] + R*sqrt(3.0)/2}};
  int cell_id = -1;
  for (const auto & x : corners)
    if (!intp->locate(x, 0., cell_id))
      partrac::fail("init_mode ", prm.get<std::string>("init_mode"), ": the ellipsoid is not inside the domain");
  Surface s(intp, corners, prm.get<Uint>("Nrw_max"));
  for (const auto & e : std::vector<std::array<Uint, 2>>{{0, 1}, {1, 2}, {2, 0}, {1, 3}, {2, 3}, {0, 3}})
    s.edges.push_back({e, dist(corners[e[0]], corners[e[1]])});
  s.faces.push_back({{0, 1, 2}, 1.});
  s.faces.push_back({{0, 3, 5}, 1.});
  s.faces.push_back({{1, 4, 3}, 1.});
  s.faces.push_back({{2, 5, 4}, 1.});
  compute_maps(s);

  // Refine to ds_max, project onto the ellipsoid, coarsen to ds_min
  Uint n_add, n_rem;
  int attempt = 0;
  do {
    n_add = refine_surface(s, prm.get<double>("ds_max"), true);
    for (Uint irw=0; irw < s.ps.N(); ++irw){
      const Vector3d x = s.ps.x(irw);
      const Vector3d nn = (x - x_c)/ (x - x_c).norm();
      const double rad = 1./sqrt(nn[0]*nn[0]/lx2 + nn[1]*nn[1]/ly2 + nn[2]*nn[2]/lz2);
      s.ps.move(irw, x_c + rad * nn);
    }
    n_rem = coarsen_surface(s, prm.get<double>("ds_min"));
    ++attempt;
  } while ((n_add > 0 || n_rem > 0) && attempt < 100);
  if (n_add > 0 || n_rem > 0)
    std::cout << "Note: the ellipsoid's remeshing stopped after " << attempt << " passes without settling" << std::endl;
  // Every node inside
  for (Uint irw=0; irw < s.ps.N(); ++irw)
    if (!intp->locate(s.ps.x(irw)))
      partrac::fail("init_mode ", prm.get<std::string>("init_mode"), ": the ellipsoid is not inside the domain");
  settle_lengths(s);
  return surface_state(s);
}

InitialState init_pairs(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen){
  const Vector3d x_min = intp->get_x_min();
  const Vector3d x_max = intp->get_x_max();
  std::vector<std::uniform_real_distribution<>> uni_dist;
  for (Uint d=0; d < 3; ++d)
    uni_dist.emplace_back(x_min[d], x_max[d]);
  std::normal_distribution<double> rnd_normal(0.0, 1.0);
  const bool single = key[0] == "pair";
  const Uint Npairs = single ? 1 : prm.get<Uint>("Nrw")/2;
  const std::array<bool, 3> spread = axis_mask(key[1]);
  // Only pairs_<dirs>_<dirs> redraws the centre
  const bool centre_moves = key[0] == "pairs" && key.size() == 3;
  const std::array<bool, 3> centre_axes = centre_moves ? axis_mask(key[2]) : std::array<bool, 3>{false, false, false};

  Vector3d centre = initial_position(prm);
  if (!centre_moves && !intp->locate(centre))
    partrac::fail("init_mode ", prm.get<std::string>("init_mode"), ": the pair centre is not inside the domain");
  InitialState state;
  Uint ipair = 0;
  Uint misses = 0;
  while (ipair < Npairs && misses < max_misses){
    for (Uint d=0; d < 3; ++d)
      if (centre_axes[d])
        centre[d] = uni_dist[d](gen);
    // Random direction in the spread's axes
    Vector3d dx(0., 0., 0.);
    for (Uint d=0; d < 3; ++d)
      if (spread[d])
        dx[d] = rnd_normal(gen);
    dx *= 0.5*prm.get<double>("ds_init")/dx.norm();
    const Vector3d x_a = centre + dx;
    const bool inside_a = intp->locate(x_a);
    const Vector3d x_b = centre - dx;
    const bool inside_b = intp->locate(x_b);
    if (inside_a && inside_b){
      state.nodes.push_back(x_a);
      state.nodes.push_back(x_b);
      state.edges.push_back({{2*ipair, 2*ipair+1}, (x_a-x_b).norm()});
      ++ipair;
      misses = 0;
    }
    else if (single){
      partrac::fail("init_mode ", prm.get<std::string>("init_mode"), ": the pair is not inside the domain");
    }
    else {
      ++misses;
    }
  }
  if (ipair < Npairs)
    partrac::fail("init_mode ", prm.get<std::string>("init_mode"), ": could not place all pairs inside the domain");
  return state;
}

InitialState init_points(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen){
  const Vector3d L = intp->get_x_max() - intp->get_x_min();
  // Sampled axes: named, with extent
  std::array<bool, 3> sampled = axis_mask(key[1]);
  for (Uint d=0; d < 3; ++d)
    sampled[d] = L[d] > 1e-12 && sampled[d];
  if (!sampled[0] && !sampled[1] && !sampled[2])
    partrac::fail("init_mode ", prm.get<std::string>("init_mode"), ": the domain has no extent along its directions");

  const WeightedGrid grid = weighted_grid(sampled, intp, prm);
  // Some weight to sample by
  const double total = std::accumulate(grid.weights.begin(), grid.weights.end(), 0.);
  if (!(std::isfinite(total) && total > 0.))
    partrac::fail("init_mode ", prm.get<std::string>("init_mode"), ": init_weight ", prm.get<std::string>("init_weight"),
                  " has no finite positive total");

  // A weighted cell, a uniform offset in it
  std::discrete_distribution<Uint> pick_cell(grid.weights.begin(), grid.weights.end());
  std::vector<std::uniform_real_distribution<>> offset;
  for (Uint d=0; d < 3; ++d)
    offset.emplace_back(-0.5*grid.h[d], 0.5*grid.h[d]);
  InitialState state;
  state.nodes = sample_inside([&]{
      Vector3d x = grid.centres[pick_cell(gen)];
      for (Uint d=0; d < 3; ++d)
        if (sampled[d])
          x[d] += offset[d](gen);
      return x;
    }, intp, prm);

  // Sorted; neighbours joined along one axis
  std::sort(state.nodes.begin(), state.nodes.end(), less_than_op());
  if (points_along_one_axis(key)){
    for (Uint irw=1; irw < state.nodes.size(); ++irw){
      const double ds0 = dist(state.nodes[irw-1], state.nodes[irw]);
      if (ds0 < 10*prm.get<double>("ds_init"))
        state.edges.push_back({{irw-1, irw}, ds0});
    }
  }
  return state;
}

InitialState init_gaussian_strip(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen){
  const double La = prm.get<double>("La");
  const double sigma0 = prm.get<double>("Lb");
  const Vector3d x0 = initial_position(prm);
  const Vector3d n = unit_vector(key[1]);
  const std::array<bool, 3> spread = axis_mask(key[2]);
  std::uniform_real_distribution<double> rnd_unit(0.0, 1.0);
  std::normal_distribution<double> rnd_normal(0.0, 1.0);
  // Ends
  const Vector3d x00 = x0 - La/2*n;
  const Vector3d x01 = x0 + La/2*n;

  InitialState state;
  state.nodes = sample_inside([&]{
      const double alpha = rnd_unit(gen);
      Vector3d x = alpha * x00 + (1.-alpha) * x01;
      for (Uint d=0; d < 3; ++d)
        if (spread[d])
          x[d] += sigma0 * rnd_normal(gen);
      return x;
    }, intp, prm);
  return state;
}

InitialState init_gaussian_circle(const std::vector<std::string>& key, std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen){
  const double R = prm.get<double>("La")/2;
  const double sigma0 = prm.get<double>("Lb");
  const Vector3d x0 = initial_position(prm);
  // In-plane axes
  Vector3d t1(0., 0., 0.);
  Vector3d t2(0., 0., 0.);
  if (key[1] == "x"){ t1[1] = 1.; t2[2] = 1.; }
  if (key[1] == "y"){ t1[0] = 1.; t2[2] = 1.; }
  if (key[1] == "z"){ t1[0] = 1.; t2[1] = 1.; }
  // Spread directions (all by default)
  const std::array<bool, 3> spread = key.size() > 2 ? axis_mask(key[2]) : std::array<bool, 3>{true, true, true};
  std::uniform_real_distribution<double> rnd_unit(0.0, 1.0);
  std::normal_distribution<double> rnd_normal(0.0, 1.0);

  InitialState state;
  state.nodes = sample_inside([&]{
      // Uniform on the disc
      double alpha1 = 1.;
      double alpha2 = 1.;
      while (alpha1*alpha1 + alpha2*alpha2 > 1){
        alpha1 = 2*rnd_unit(gen)-1;
        alpha2 = 2*rnd_unit(gen)-1;
      }
      Vector3d x = x0 + R * (alpha1 * t1 + alpha2 * t2);
      for (Uint d=0; d < 3; ++d)
        if (spread[d])
          x[d] += sigma0 * rnd_normal(gen);
      return x;
    }, intp, prm);
  return state;
}

InitialState init_file(const std::string& path, std::shared_ptr<Interpol> intp, const partrac::Params& prm){
  verify_file_exists(path);
  hsize_t dims_nodes[2];
  std::vector<double> nodes_buf;
  try {
    H5::H5File h5file(path, H5F_ACC_RDONLY);
    H5::DataSet dset_nodes = h5file.openDataSet("nodes");
    H5::DataSpace dspace_nodes = dset_nodes.getSpace();
    // One row a point, one to three coordinates
    if (dspace_nodes.getSimpleExtentNdims() != 2)
      partrac::fail(path, ": nodes is not a two-dimensional array, one row a point");
    dspace_nodes.getSimpleExtentDims(dims_nodes, NULL);
    if (dims_nodes[1] < 1 || dims_nodes[1] > 3)
      partrac::fail(path, ": nodes has ", dims_nodes[1], " columns, not one to three coordinates");
    nodes_buf.resize(dims_nodes[0]*dims_nodes[1]);
    dset_nodes.read(nodes_buf.data(), H5::PredType::NATIVE_DOUBLE, dspace_nodes, dspace_nodes);
    h5file.close();
  } catch (const H5::Exception&){
    partrac::fail(path, ": cannot read the dataset 'nodes'");
  }

  // Coordinates the file leaves out from (x0, y0, z0)
  const Vector3d x0 = initial_position(prm);
  const double t0 = prm.get<double>("t0");
  std::vector<Vector3d> points;
  std::vector<bool> inside;
  int cell_id = -1;
  for (Uint i=0; i < dims_nodes[0]; ++i){
    Vector3d x = x0;
    for (Uint j=0; j < dims_nodes[1]; ++j)
      x[j] = nodes_buf[i * dims_nodes[1] + j];
    points.push_back(x);
    inside.push_back(intp->locate(x, t0, cell_id));
  }
  InitialState state = polyline_inside(points, inside);
  if (state.nodes.empty())
    partrac::fail("init_mode ", prm.get<std::string>("init_mode"), ": no points inside the domain");
  return state;
}

InitialState set_initial_state(std::shared_ptr<Interpol> intp, const partrac::Params& prm, std::mt19937& gen){
  const std::string init_mode = prm.get<std::string>("init_mode");
  // file:<path>, whose path may hold anything
  if (init_mode_is_file(init_mode))
    return init_file(init_mode.substr(5), intp, prm);
  const std::vector<std::string> key = split_string(init_mode, "_");
  const InitMode* mode = find_init_mode(key[0]);
  if (!mode)
    partrac::fail("init_mode ", init_mode, ": no such mode");
  return mode->build(key, intp, prm, gen);
}
