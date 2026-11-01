// wall_facets and wall_markers: the no-slip walls of wall_upgrade as a
// meshio facet file names them, matched to the loaded mesh by coordinates, in
// the stamped loaders (a dolfin checkpoint, P1 and P2, and XDMF) and the
// frequency loader, on dolfin's unit square and cube. The facet file is
// written from dolfin's own numbering, which the loaders renumber. Without the
// file, a time series' walls are those at rest in its first stamp with flow.
// The warnings: of walls from the velocity, and of named walls that move. The
// point matching on a strongly graded mesh, and mode openfoam's warnings.
#include <iostream>
#include <sstream>
#include <string>

namespace {

// What was written to cerr while f ran
template<typename F>
std::string captured_cerr(F&& f){
  std::stringstream buf;
  std::streambuf* old = std::cerr.rdbuf(buf.rdbuf());
  try { f(); }
  catch (...){ std::cerr.rdbuf(old); throw; }
  std::cerr.rdbuf(old);
  return buf.str();
}

}  // namespace

#ifdef USE_DOLFIN
#include <catch2/catch.hpp>

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <dolfin.h>
#include "H5Cpp.h"

#include "Error.hpp"
#include "Params.hpp"
#include "SimplexInterpol.hpp"
#include "XDMFInterpol.hpp"
#include "case_dir.hpp"
#include "held_eval.hpp"
#include "interpolators.hpp"
#include "loader_params.hpp"
#include "region_cases.hpp"
#include "wall_facets.hpp"

namespace {

const std::vector<bool> no_axis = {false, false, false};
const std::vector<bool> x_axis = {true, false, false};

// A facet file's contents: points of cols coordinates, facets of D points, a marker each
struct Facets {
  std::vector<double> points;
  std::vector<std::int32_t> facets;
  std::vector<std::int32_t> markers;
  int cols = 3;
};

// The unit box's exterior facets in dolfin's numbering, marked by side: x =
// 0 and 1 are 1 and 2, y 3 and 4, z 5 and 6
template<typename Cell>
Facets box_facets(const std::size_t n, const int cols){
  constexpr int D = dim_of<Cell>;
  auto mesh = unit_mesh<Cell>(n);
  const std::vector<double>& x = mesh->coordinates();
  Facets f;
  f.cols = cols;
  for (std::size_t v = 0; v < x.size()/D; ++v)
    for (int d = 0; d < cols; ++d) f.points.push_back(d < D ? x[v*D + std::size_t(d)] : 0.);
  for (dolfin::FacetIterator e(*mesh); !e.end(); ++e){
    if (!e->exterior()) continue;
    const unsigned* v = e->entities(0);
    int marker = 0;
    for (int d = 0; d < D; ++d){
      bool lo = true, hi = true;
      for (int k = 0; k < D; ++k){
        lo = lo && x[v[k]*D + unsigned(d)] == 0.;
        hi = hi && x[v[k]*D + unsigned(d)] == 1.;
      }
      if (lo) marker = 1 + 2*d;
      if (hi) marker = 2 + 2*d;
    }
    REQUIRE(marker > 0);
    for (int k = 0; k < D; ++k) f.facets.push_back(std::int32_t(v[k]));
    f.markers.push_back(marker);
  }
  return f;
}

// The point nearest p
std::int32_t point_at(const Facets& f, const Vector3d& p){
  std::int32_t best = 0;
  double d_best = 1e300;
  for (std::size_t v = 0; v < f.points.size()/std::size_t(f.cols); ++v){
    double d = 0.;
    for (int k = 0; k < std::min(f.cols, 3); ++k) d += std::pow(f.points[v*std::size_t(f.cols) + std::size_t(k)] - p[k], 2);
    if (d < d_best){ d_best = d; best = std::int32_t(v); }
  }
  return best;
}

// f as meshio writes it: name.xdmf, its data in name.h5 as /data0 to /data2
void write_facets(const CaseDir& c, const std::string& name, const Facets& f, const int D){
  const hsize_t np = f.points.size()/std::size_t(f.cols), nf = f.markers.size();
  {
    H5::H5File h(c.file(name + ".h5"), H5F_ACC_TRUNC);
    const hsize_t dp[2] = {np, hsize_t(f.cols)}, df[2] = {nf, hsize_t(D)}, dm[1] = {nf};
    h.createDataSet("data0", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dp))
      .write(f.points.data(), H5::PredType::NATIVE_DOUBLE);
    h.createDataSet("data1", H5::PredType::NATIVE_INT32, H5::DataSpace(2, df))
      .write(f.facets.data(), H5::PredType::NATIVE_INT32);
    h.createDataSet("data2", H5::PredType::NATIVE_INT32, H5::DataSpace(1, dm))
      .write(f.markers.data(), H5::PredType::NATIVE_INT32);
  }
  std::ofstream(c.file(name + ".xdmf"))
    << "<Xdmf Version=\"3.0\"><Domain><Grid Name=\"Grid\"><Geometry GeometryType=\""
    << (f.cols == 3 ? "XYZ" : "XY") << "\"><DataItem DataType=\"Float\" Dimensions=\"" << np << " " << f.cols
    << "\" Format=\"HDF\" Precision=\"8\">" << name << ".h5:/data0</DataItem></Geometry><Topology TopologyType=\""
    << (D == 3 ? "Triangle" : "Polyline") << "\" NumberOfElements=\"" << nf << "\" NodesPerElement=\"" << D
    << "\"><DataItem DataType=\"Int\" Dimensions=\"" << nf << " " << D << "\" Format=\"HDF\" Precision=\"4\">"
    << name << ".h5:/data1</DataItem></Topology><Attribute Name=\"markers\" AttributeType=\"Scalar\" "
    << "Center=\"Cell\"><DataItem DataType=\"Int\" Dimensions=\"" << nf << "\" Format=\"HDF\" Precision=\"4\">"
    << name << ".h5:/data2</DataItem></Attribute></Grid></Domain></Xdmf>\n";
}

// A copy of parameter file from, as to, with extra appended
std::string with(const CaseDir& c, const std::string& from, const std::string& to, const std::string& extra){
  std::ifstream in(c.file(from));
  std::stringstream s;
  s << in.rdbuf();
  std::ofstream(c.file(to)) << s.str() << extra;
  return c.file(to);
}

// Random points of the box, then points on its sides
template<int D>
std::vector<Vector3d> box_and_sides(const std::size_t n){
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> uni(0.001, 0.999);
  std::vector<Vector3d> pts;
  for (std::size_t i = 0; i < 2*n; ++i){
    Vector3d x = Vector3d::Zero();
    for (int d = 0; d < D; ++d) x[d] = uni(rng);
    if (i >= n) x[int(i % std::size_t(D))] = i % 2 ? 0. : 1.;
    pts.push_back(x);
  }
  return pts;
}

// The fields at the points, at t: located and, where held, held in the cell
template<typename I>
std::vector<PointValues> sample(I& intp, const std::vector<Vector3d>& pts, const double t, const bool hold = false){
  intp.set_int_order(2);
  intp.set_needs_gradient(true);
  intp.update(t);
  std::vector<PointValues> out;
  for (const Vector3d& x : pts){
    CellPos pos;
    REQUIRE(intp.locate(x, t, pos));
    PointValues v(1.);
    if (hold){
      HeldEval<I> ev(intp, v);
      ev.bind(intp.region_in(pos.id, x));
      ev.at(x, t);
    }
    else intp.evaluate_motion(x, t, pos, v);
    out.push_back(v);
  }
  return out;
}

bool same(const PointValues& a, const PointValues& b){
  return a.U == b.U && a.A == b.A && a.gradU == b.gradU && a.gradA == b.gradA;
}

// The distance to the nearest side of those a mask names, by axis
template<int D>
double to_sides(const Vector3d& x, const unsigned axes){
  double d = 1.;
  for (int k = 0; k < D; ++k)
    if (axes >> k & 1) d = std::min({d, x[k], 1. - x[k]});
  return d;
}

// Every interpolator of one case and loader, by its parameter file
template<typename Cell>
using Loaders = std::vector<std::pair<std::string, std::function<std::vector<PointValues>(
                  const std::string&, const std::vector<Vector3d>&, double)>>>;

template<typename Cell>
Loaders<Cell> loaders(const std::string& u_el){
  Loaders<Cell> out;
  const auto stamped = [](const std::string& p, const std::vector<Vector3d>& pts, const double t){
    SimplexInterpol<Cell> i(p); return sample(i, pts, t); };
  const auto held = [](const std::string& p, const std::vector<Vector3d>& pts, const double t){
    SimplexInterpol<Cell> i(p); return sample(i, pts, t, true); };
  const auto freq = [](const std::string& p, const std::vector<Vector3d>& pts, const double t){
    SimplexFreqInterpol<Cell> i(p); return sample(i, pts, t); };
  const auto xdmf = [](const std::string& p, const std::vector<Vector3d>& pts, const double t){
    XDMFInterpol<Cell> i(p); return sample(i, pts, t); };
  out.push_back({"h5_params.dat", stamped});
  out.push_back({"h5_params.dat", held});
  out.push_back({"freq_params.dat", freq});
  if (u_el == "P1") out.push_back({"xdmf_params.dat", xdmf});
  return out;
}

// Both stamps (or components) of a case with field, and its parameter files
template<typename Cell>
void write_case(const CaseDir& c, const std::size_t n, const std::string& u_el, const std::vector<bool>& per,
                const Field& field){
  write_stamped<Cell>(c, n, u_el, per, true, {}, field);
  write_freq_of_stamps(c, u_el, per, "1 0.2 0.7 up_0.h5\n2 0.3 -1.3 up_1.h5\n", "");
}

// At rest on every side of the box
template<int D>
Vector3d at_rest(const Vector3d& x, const int k){
  double b = 1.;
  for (int d = 0; d < D; ++d) b *= x[d]*(1. - x[d]);
  return b*(1. + double(k))*Vector3d(4.*(1. + x[1]), 1. + x[0], D == 3 ? 0.5 + x[0]*x[1] : 0.);
}

// Through the sides x = 0 and 1, at rest on the others; stamp 0 at rest
template<int D>
Vector3d through_x(const Vector3d& x, const int k){
  double b = x[1]*(1. - x[1]);
  if (D == 3) b *= 4.*x[2]*(1. - x[2]);
  return k == 0 ? Vector3d::Zero() : Vector3d(b*(1. + 0.5*x[0]), 0.3*b*x[0], D == 3 ? 0.2*b*x[0] : 0.);
}

// wall_facets naming the sides off the periodic axes reads as the walls found
// at rest, in every loader, bit for bit
template<typename Cell>
void same_as_rest(const std::size_t n, const std::string& u_el, const std::vector<bool>& per, const int cols){
  constexpr int D = dim_of<Cell>;
  CaseDir c("wf_same");
  write_case<Cell>(c, n, u_el, per, at_rest<D>);
  write_facets(c, "facets", box_facets<Cell>(n, cols), D);
  const std::string up = "wall_upgrade=true\n";
  const std::string extra = up + (per[0] ? "wall_facets=facets.xdmf\nwall_markers=3:6\n"
                                         : "wall_facets=facets.xdmf\nwall_markers=1:6\n");
  const std::vector<Vector3d> pts = box_and_sides<D>(300);
  const double t = 0.4;
  for (const auto& [params, run] : loaders<Cell>(u_el)){
    INFO(params);
    const auto a = run(with(c, params, "up_" + params, up), pts, t);
    const auto b = run(with(c, params, "walls_" + params, extra), pts, t);
    const auto p = run(with(c, params, "plain_" + params, "wall_upgrade=false\n"), pts, t);
    std::size_t differ = 0, raised = 0;
    for (std::size_t i = 0; i < pts.size(); ++i){
      differ += !same(a[i], b[i]);
      raised += !same(a[i], p[i]);
    }
    CHECK(differ == 0);
    CHECK(raised > 0);
  }
}

// Walls on the sides off x alone, an inlet and an outlet on x = 0 and 1, the
// first stamp at rest: the cells without a vertex on a wall read the plain
// field bit for bit, those at the inlet and outlet included; on the walls the
// field is zero, and near them it is raised
template<typename Cell>
void inlet_outlet(const std::size_t n, const std::string& u_el){
  constexpr int D = dim_of<Cell>;
  CaseDir c("wf_inlet");
  write_case<Cell>(c, n, u_el, no_axis, through_x<D>);
  write_facets(c, "facets", box_facets<Cell>(n, D), D);
  std::vector<Vector3d> pts = box_and_sides<D>(300);
  // the inlet and outlet away from the walls
  std::mt19937 rng(5);
  std::uniform_real_distribution<double> mid(0.3, 0.7);
  for (int i = 0; i < 100; ++i){
    Vector3d x(i % 2 ? 0. : 1., mid(rng), D == 3 ? mid(rng) : 0.);
    pts.push_back(x);
  }
  const unsigned walls = D == 3 ? 6u : 2u;
  const double reach = 1./double(n) + 1e-9, t = 0.4;
  for (const auto& [params, run] : loaders<Cell>(u_el)){
    INFO(params);
    const auto a = run(with(c, params, "walls_" + params,
                            "wall_upgrade=true\nwall_facets=facets.xdmf\nwall_markers=3:6\n"), pts, t);
    const auto p = run(with(c, params, "plain_" + params, "wall_upgrade=false\n"), pts, t);
    std::size_t off = 0, inlet = 0, on = 0, raised = 0;
    double scale = 0.;
    for (const PointValues& v : p) scale = std::max(scale, v.U.norm());
    REQUIRE(scale > 0.);
    for (std::size_t i = 0; i < pts.size(); ++i){
      const Vector3d& x = pts[i];
      const double d = to_sides<D>(x, walls);
      if (d > reach){
        CHECK(same(a[i], p[i]));
        ++off;
        if (x[0] == 0. || x[0] == 1.) ++inlet;
      }
      else if (d == 0.){
        CHECK(a[i].U.norm() < 1e-14*scale);
        ++on;
      }
      else raised += !same(a[i], p[i]);
    }
    CHECK(inlet >= 100);
    CHECK(off > inlet);
    CHECK(on > 50);
    CHECK(raised > 0);
  }
}

// cell_type from the named walls: with wall_facets naming the sides off x, a
// cell at the inlet's middle is 0; without, its exterior facet makes it 1
template<typename Cell, typename I>
void inlet_cell_type(const std::size_t n, const std::string& u_el, const std::string& params){
  constexpr int D = dim_of<Cell>;
  CaseDir c("wf_type");
  write_case<Cell>(c, n, u_el, no_axis, through_x<D>);
  write_facets(c, "facets", box_facets<Cell>(n, D), D);
  const Vector3d x(0.01, 0.5, D == 3 ? 0.5 : 0.);
  const auto type_at = [&](const std::string& p){
    I intp(p);
    intp.update(0.4);
    CellPos pos;
    REQUIRE(intp.locate(x, 0.4, pos));
    PointValues v(1.);
    intp.evaluate(x, 0.4, pos, v);
    return v.get_cell_type();
  };
  CHECK(type_at(with(c, params, "walls_" + params,
                     "wall_upgrade=true\nwall_facets=facets.xdmf\nwall_markers=3:6\n")) == 0);
  CHECK(type_at(with(c, params, "up_" + params, "wall_upgrade=true\n")) == 1);
}

template<typename M>
double largest(const Eigen::MatrixBase<M>& m){ return m.cwiseAbs().maxCoeff(); }

// Without wall_facets, stamps whose first is at rest (I reads params): their
// walls are fixed by the second, and the field between them is the second's
// raised field times the weight, its rate that field; reversed, the first
// fixes them and the second at rest is kept. Two stamps at rest read zero.
template<typename Cell, typename I>
void first_with_flow(const std::size_t n, const std::string& u_el, const std::string& params){
  constexpr int D = dim_of<Cell>;
  const Field flow = [](const Vector3d& x, int){ return through_x<D>(x, 1); };
  const Field rest_first = [](const Vector3d& x, const int k){ return through_x<D>(x, k); };
  const Field rest_last = [](const Vector3d& x, const int k){ return through_x<D>(x, 1 - k); };
  const Field none = [](const Vector3d&, int){ return Vector3d::Zero(); };
  CaseDir cf("wf_flow"), cz("wf_zero_first"), cr("wf_zero_last"), cn("wf_none");
  write_case<Cell>(cf, n, u_el, no_axis, flow);
  write_case<Cell>(cz, n, u_el, no_axis, rest_first);
  write_case<Cell>(cr, n, u_el, no_axis, rest_last);
  write_case<Cell>(cn, n, u_el, no_axis, none);
  const std::vector<Vector3d> pts = box_and_sides<D>(300);
  const std::string up = "wall_upgrade=true\n";
  for (const bool hold : {false, true})
    for (const double t : {0., 0.4, 0.9}){
      INFO("t = " << t << (hold ? ", held" : ""));
      I f(with(cf, params, "up.dat", up)), z(with(cz, params, "up.dat", up)), r(with(cr, params, "up.dat", up)),
        e(with(cn, params, "up.dat", up));
      const auto ref = sample(f, pts, 0.5, hold), a = sample(z, pts, t, hold), b = sample(r, pts, t, hold);
      const auto zero = sample(e, pts, t, hold);
      double scale = 0., err = 0.;
      for (std::size_t i = 0; i < pts.size(); ++i){
        const PointValues& q = ref[i];
        scale = std::max({scale, largest(q.U), largest(q.gradU)});
        err = std::max({err, largest(a[i].U - t*q.U), largest(a[i].gradU - t*q.gradU), largest(a[i].A - q.U),
                        largest(a[i].gradA - q.gradU), largest(b[i].U - (1. - t)*q.U),
                        largest(b[i].gradU - (1. - t)*q.gradU), largest(b[i].A + q.U),
                        largest(b[i].gradA + q.gradU)});
        CHECK(zero[i].U.isZero(0.));
        CHECK(zero[i].gradU.isZero(0.));
      }
      REQUIRE(scale > 0.);
      CHECK(err <= 1e-13*scale);
    }
  // The correction is on
  I f(cf.file("up.dat")), p(with(cf, params, "plain.dat", "wall_upgrade=false\n"));
  const auto a = sample(f, pts, 0.5), b = sample(p, pts, 0.5);
  std::size_t raised = 0;
  for (std::size_t i = 0; i < pts.size(); ++i) raised += !same(a[i], b[i]);
  CHECK(raised > 0);
}

// The message of what loading params throws
template<typename I>
std::string refusal(const std::string& params){
  try { I intp(params); }
  catch (const partrac::Error& e){ return e.what(); }
  return "accepted";
}

// A schema's refusal of a parameter file
std::string schema_refusal(const partrac::Schema& s, const std::string& path){
  try { s.parse_file(path); }
  catch (const partrac::ParamError& e){
    std::string out;
    for (const auto& p : e.problems()) out += p + "\n";
    return out;
  }
  return "accepted";
}

template<typename Cell>
void refusals(const std::size_t n){
  constexpr int D = dim_of<Cell>;
  using I = SimplexInterpol<Cell>;
  CaseDir c("wf_refused");
  write_stamped<Cell>(c, n, "P2", x_axis, true, {}, at_rest<D>);
  Facets f = box_facets<Cell>(n, 3);
  write_facets(c, "facets", f, D);
  // a periodic side's facet
  std::string why = refusal<I>(with(c, "h5_params.dat", "a.dat", "wall_upgrade=true\nwall_facets=facets.xdmf\nwall_markers=1:6\n"));
  INFO(why);
  CHECK((why.find("marker 1") != std::string::npos || why.find("marker 2") != std::string::npos));
  CHECK(why.find("periodic side") != std::string::npos);
  CHECK(refusal<I>(with(c, "h5_params.dat", "b.dat", "wall_upgrade=true\nwall_facets=facets.xdmf\nwall_markers=3:6\n")) == "accepted");
  CHECK(refusal<I>(with(c, "h5_params.dat", "b2.dat", "wall_upgrade=true\nwall_facets=facets.xdmf\nwall_markers=50\n"))
        .find("selects none") != std::string::npos);
  // boundary points that are no facet's, marker 9: refused when selected alone
  Facets g = f;
  g.facets.push_back(point_at(g, Vector3d(0., 0., 0.)));
  g.facets.push_back(point_at(g, Vector3d(0., 1., 0.)));
  if (D == 3) g.facets.push_back(point_at(g, Vector3d(1., 1., 0.)));
  g.markers.push_back(9);
  write_facets(c, "extra", g, D);
  why = refusal<I>(with(c, "h5_params.dat", "c.dat", "wall_upgrade=true\nwall_facets=extra.xdmf\nwall_markers=3:9\n"));
  CHECK(why.find("marker 9") != std::string::npos);
  CHECK(why.find("matches no exterior facet") != std::string::npos);
  CHECK(refusal<I>(with(c, "h5_params.dat", "d.dat", "wall_upgrade=true\nwall_facets=extra.xdmf\nwall_markers=3:6\n")) == "accepted");
  // a point off every vertex, on a facet of side y = 0
  Facets h = f;
  // its first point inside the side, on no other
  const auto inside = [&](const std::size_t v){
    const double* x = &h.points[v*3];
    return x[0] > 0. && x[0] < 1. && (D == 2 || (x[2] > 0. && x[2] < 1.));
  };
  std::size_t j = 0;
  while (h.markers[j] != 3 || !inside(std::size_t(h.facets[j*D]))) ++j;
  const std::size_t q = std::size_t(h.facets[j*D]);
  h.points[q*3] += 0.3/double(n);
  write_facets(c, "moved", h, D);
  why = refusal<I>(with(c, "h5_params.dat", "e.dat", "wall_upgrade=true\nwall_facets=moved.xdmf\nwall_markers=3:6\n"));
  CHECK(why.find("point " + std::to_string(q)) != std::string::npos);
  CHECK(why.find("matches no vertex") != std::string::npos);
  CHECK(refusal<I>(with(c, "h5_params.dat", "e2.dat", "wall_upgrade=true\nwall_facets=moved.xdmf\nwall_markers=4:6\n")) == "accepted");
  // the schemas
  const partrac::Schema s = dolfin_h5_schema(D == 3 ? "tet" : "triangle");
  CHECK(schema_refusal(s, with(c, "h5_params.dat", "f.dat", "wall_markers=3\n")).find("wall_markers") != std::string::npos);
  // the walls are always named by marker
  CHECK(schema_refusal(s, with(c, "h5_params.dat", "f2.dat", "wall_upgrade=true\nwall_facets=facets.xdmf\n"))
          .find("wall_facets needs wall_markers") != std::string::npos);
  CHECK(schema_refusal(s, with(c, "h5_params.dat", "f3.dat", "wall_upgrade=true\nwall_facets=facets.xdmf\nwall_markers=\n"))
          .find("wall_markers is empty") != std::string::npos);
  CHECK(schema_refusal(s, with(c, "h5_params.dat", "g.dat", "wall_facets=facets.xdmf\nwall_markers=1:6\nwall_upgrade=false\n"))
        .find("wall_upgrade, which is false") != std::string::npos);
  std::ofstream(c.file("foam.dat")) << "wall_facets=facets.xdmf\n";
  CHECK(schema_refusal(openfoam_schema(), c.file("foam.dat")).find("mode openfoam") != std::string::npos);
}

std::size_t occurrences(const std::string& s, const std::string& what){
  std::size_t n = 0;
  for (std::size_t at = s.find(what); at != std::string::npos; at = s.find(what, at + 1)) ++n;
  return n;
}

const std::string velocity_walls = "Warning: wall_upgrade takes the walls from the velocity";
const std::string moving_walls = "Warning: the walls wall_facets names move";

// Without wall_facets every loader warns once that its walls come from the
// velocity; with it, never, and walls at rest do not move
template<typename Cell>
void velocity_walls_warning(const std::size_t n, const std::string& u_el){
  constexpr int D = dim_of<Cell>;
  CaseDir c("wf_warn");
  write_case<Cell>(c, n, u_el, no_axis, at_rest<D>);
  write_facets(c, "facets", box_facets<Cell>(n, D), D);
  const std::vector<Vector3d> pts = box_and_sides<D>(10);
  for (const auto& [params, run] : loaders<Cell>(u_el)){
    INFO(params);
    const std::string a = captured_cerr([&]{ run(with(c, params, "up_" + params, "wall_upgrade=true\n"), pts, 0.4); });
    const std::string b = captured_cerr([&]{
      run(with(c, params, "walls_" + params, "wall_upgrade=true\nwall_facets=facets.xdmf\nwall_markers=1:6\n"), pts, 0.4); });
    CHECK(occurrences(a, velocity_walls) == 1);
    CHECK(occurrences(b, velocity_walls) == 0);
    CHECK(occurrences(a + b, moving_walls) == 0);
  }
}

// Named walls meeting an inlet at x = 0, the first stamp at rest: clean, no
// warning; with a plug on the inlet's nodes in the second, one warning in
// every loader, naming the wall nodes on the inlet (those of the walls'
// entities: the two corners in 2D, the rim's in 3D)
template<typename Cell>
void moving_walls_warning(const std::size_t n, const std::string& u_el){
  constexpr int D = dim_of<Cell>;
  const Field plug = [](const Vector3d& x, const int k) -> Vector3d {
    return through_x<D>(x, k) + Vector3d(k == 1 && x[0] < 1e-9 ? 1. : 0., 0., 0.);
  };
  CaseDir c("wf_clean"), m("wf_moving");
  write_case<Cell>(c, n, u_el, no_axis, through_x<D>);
  write_case<Cell>(m, n, u_el, no_axis, plug);
  write_facets(c, "facets", box_facets<Cell>(n, D), D);
  write_facets(m, "facets", box_facets<Cell>(n, D), D);
  const std::size_t nodes = D == 2 ? 2 : (u_el == "P2" ? 8 : 4)*n;
  const std::string walls = "wall_upgrade=true\nwall_facets=facets.xdmf\nwall_markers=3:6\n";
  const std::vector<Vector3d> pts = box_and_sides<D>(10);
  for (const auto& [params, run] : loaders<Cell>(u_el)){
    INFO(params);
    const std::string a = captured_cerr([&]{ run(with(c, params, "walls_" + params, walls), pts, 0.4); });
    const std::string b = captured_cerr([&]{ run(with(m, params, "walls_" + params, walls), pts, 0.4); });
    CHECK(occurrences(a, moving_walls) == 0);
    CHECK(occurrences(b, moving_walls) == 1);
    CHECK(b.find("move at " + std::to_string(nodes) + " velocity nodes") != std::string::npos);
    CHECK(b.find(params == "freq_params.dat" ? "in the frequency components" : "in the first stamp with flow")
          != std::string::npos);
    CHECK(b.find("|u| up to 0.8 of the largest |u|; they are kept as no-slip walls") != std::string::npos);
  }
}

// Frequency components on named walls at rest but for a noise-level one that
// moves them: judged by their sum, no warning
template<typename Cell>
void noise_component(const std::size_t n, const std::string& u_el){
  constexpr int D = dim_of<Cell>;
  const Field f = [](const Vector3d& x, const int k) -> Vector3d {
    return k == 0 ? at_rest<D>(x, 0) : Vector3d(1e-9, 1e-9, D == 3 ? 1e-9 : 0.);
  };
  CaseDir c("wf_noise");
  write_case<Cell>(c, n, u_el, no_axis, f);
  write_facets(c, "facets", box_facets<Cell>(n, D), D);
  const std::string walls = "wall_upgrade=true\nwall_facets=facets.xdmf\nwall_markers=1:6\n";
  const std::vector<Vector3d> pts = box_and_sides<D>(10);
  const std::string a = captured_cerr([&]{
    SimplexFreqInterpol<Cell> i(with(c, "freq_params.dat", "walls_freq.dat", walls));
    sample(i, pts, 0.4);
  });
  CHECK(occurrences(a, moving_walls) == 0);
}

}  // namespace

TEST_CASE("wall_markers selects values and inclusive ranges", "[wall_facets]"){
  const wall_facets::Markers range("7:106"), list("3,5, 7:106");
  CHECK(range.selects(7));
  CHECK(range.selects(106));
  CHECK_FALSE(range.selects(6));
  CHECK_FALSE(range.selects(107));
  for (const int m : {3, 5, 7, 50, 106}) CHECK(list.selects(m));
  for (const int m : {1, 2, 4, 6, 107}) CHECK_FALSE(list.selects(m));
  CHECK(wall_facets::Markers("-2:-1").selects(-2));
  for (const std::string bad : {"", " ", "a", "5:3", "3,,4", "3,", "1:2:3", "4.5"})
    CHECK_THROWS_AS(wall_facets::Markers(bad), partrac::Error);
}

TEST_CASE("wall_facets naming the walls at rest reads as the walls found at rest", "[wall_facets]"){
  // Checkpoint P1 and P2 (located and held), XDMF P1 and the frequency loader,
  // the facet file's points with three coordinates or the mesh's two, and the
  // sides of a periodic axis left out by wall_markers
  SECTION("P1 triangles"){ same_as_rest<Triangle>(6, "P1", no_axis, 2); }
  SECTION("P2 triangles"){ same_as_rest<Triangle>(6, "P2", no_axis, 3); }
  SECTION("P2 triangles, periodic along x"){ same_as_rest<Triangle>(6, "P2", x_axis, 2); }
  SECTION("P1 tets"){ same_as_rest<Tet>(3, "P1", no_axis, 3); }
  SECTION("P2 tets"){ same_as_rest<Tet>(3, "P2", no_axis, 3); }
  SECTION("P2 tets, periodic along x"){ same_as_rest<Tet>(3, "P2", x_axis, 3); }
}

TEST_CASE("With wall_facets an inlet and an outlet are no walls and a first stamp at rest is read",
          "[wall_facets]"){
  SECTION("P1 triangles"){ inlet_outlet<Triangle>(6, "P1"); }
  SECTION("P2 triangles"){ inlet_outlet<Triangle>(6, "P2"); }
  SECTION("P1 tets"){ inlet_outlet<Tet>(4, "P1"); }
  SECTION("P2 tets"){ inlet_outlet<Tet>(4, "P2"); }
}

TEST_CASE("With wall_facets an inlet cell is no wall cell", "[wall_facets]"){
  SECTION("P1 triangles"){ inlet_cell_type<Triangle, SimplexInterpol<Triangle>>(6, "P1", "h5_params.dat"); }
  SECTION("P2 tets"){ inlet_cell_type<Tet, SimplexInterpol<Tet>>(4, "P2", "h5_params.dat"); }
  SECTION("P1 tets, XDMF"){ inlet_cell_type<Tet, XDMFInterpol<Tet>>(4, "P1", "xdmf_params.dat"); }
}

TEST_CASE("Without wall_facets the first stamp with flow fixes the walls", "[wall_facets]"){
  SECTION("P2 triangles"){ first_with_flow<Triangle, SimplexInterpol<Triangle>>(6, "P2", "h5_params.dat"); }
  SECTION("P2 tets"){ first_with_flow<Tet, SimplexInterpol<Tet>>(4, "P2", "h5_params.dat"); }
  SECTION("P1 triangles"){ first_with_flow<Triangle, SimplexInterpol<Triangle>>(6, "P1", "h5_params.dat"); }
  SECTION("P1 tets"){ first_with_flow<Tet, SimplexInterpol<Tet>>(4, "P1", "h5_params.dat"); }
  SECTION("P1 triangles, XDMF"){ first_with_flow<Triangle, XDMFInterpol<Triangle>>(6, "P1", "xdmf_params.dat"); }
  SECTION("P1 tets, XDMF"){ first_with_flow<Tet, XDMFInterpol<Tet>>(4, "P1", "xdmf_params.dat"); }
}

TEST_CASE("wall_facets refuses a facet off the mesh's walls and a point off its vertices", "[wall_facets]"){
  SECTION("triangles"){ refusals<Triangle>(6); }
  SECTION("tets"){ refusals<Tet>(3); }
}

TEST_CASE("Walls from the velocity are warned about once, and never with wall_facets", "[wall_facets]"){
  SECTION("P1 triangles"){ velocity_walls_warning<Triangle>(4, "P1"); }
  SECTION("P2 triangles"){ velocity_walls_warning<Triangle>(4, "P2"); }
  SECTION("P1 tets"){ velocity_walls_warning<Tet>(2, "P1"); }
  SECTION("P2 tets"){ velocity_walls_warning<Tet>(2, "P2"); }
}

TEST_CASE("Named walls carrying velocity are warned about once", "[wall_facets]"){
  SECTION("P1 triangles"){ moving_walls_warning<Triangle>(6, "P1"); }
  SECTION("P2 triangles"){ moving_walls_warning<Triangle>(6, "P2"); }
  SECTION("P1 tets"){ moving_walls_warning<Tet>(4, "P1"); }
  SECTION("P2 tets"){ moving_walls_warning<Tet>(4, "P2"); }
}

TEST_CASE("Frequency components are judged moving by their weighted sum", "[wall_facets]"){
  SECTION("P1 triangles"){ noise_component<Triangle>(6, "P1"); }
  SECTION("P2 triangles"){ noise_component<Triangle>(6, "P2"); }
  SECTION("P1 tets"){ noise_component<Tet>(4, "P1"); }
  SECTION("P2 tets"){ noise_component<Tet>(4, "P2"); }
}

#endif

#include <catch2/catch.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "Error.hpp"
#include "OpenFoamInterpol.hpp"
#include "cell_walk.hpp"
#include "mesh_tables.hpp"
#include "openfoam_load.hpp"
#include "wall_facets.hpp"

namespace {

// The rectangle [0, 1] x [0, 0.5], two triangles a quad of ny rows: along x,
// nf cells of width w0, then widths growing by g up to x = 1
struct Strip {
  std::vector<double> coords;
  std::vector<std::uint32_t> topo;
  std::size_t ncells = 0;
  std::vector<std::int32_t> across;
};

Strip graded_strip(const std::size_t nf, const double w0, const double g, const std::size_t ny){
  std::vector<double> xs = {0.};
  for (std::size_t i = 0; i < nf; ++i) xs.push_back(double(i + 1)*w0);
  for (double w = w0*g; xs.back() + 1.5*w < 1.; w *= g) xs.push_back(xs.back() + w);
  xs.push_back(1.);
  Strip s;
  const std::size_t nx = xs.size() - 1;
  for (std::size_t j = 0; j <= ny; ++j)
    for (const double x : xs){
      s.coords.push_back(x);
      s.coords.push_back(0.5*double(j)/double(ny));
    }
  const auto id = [&](const std::size_t i, const std::size_t j){ return std::uint32_t(j*(nx + 1) + i); };
  for (std::size_t j = 0; j < ny; ++j)
    for (std::size_t i = 0; i < nx; ++i){
      s.topo.insert(s.topo.end(), {id(i, j), id(i + 1, j), id(i + 1, j + 1)});
      s.topo.insert(s.topo.end(), {id(i, j), id(i + 1, j + 1), id(i, j + 1)});
    }
  s.ncells = s.topo.size()/3;
  mesh_tables::build_facet_neighbours<3>(s.topo, s.ncells, s.coords, 2, {false, false, false}, Vector3d::Zero(),
                                         Vector3d(1., 0.5, 0.), 1e-12, s.across);
  return s;
}

// Its exterior facets as a facet file of the mesh's points, marker 1 on y = 0
// and 2 elsewhere, each point moved by frac of the shortest boundary edge
// there in a random direction, every other one left exact
wall_facets::File strip_facets(const Strip& s, const double frac, std::vector<double>& h){
  wall_facets::File f;
  f.points = s.coords;
  const std::size_t nverts = s.coords.size()/2;
  h.assign(nverts, std::numeric_limits<double>::infinity());
  for (std::size_t q = 0; q < s.ncells*3; ++q){
    if (s.across[q] != facet_wall) continue;
    const std::uint32_t* row = &s.topo[(q/3)*3];
    const std::uint32_t a = row[(q % 3 + 1) % 3], b = row[(q % 3 + 2) % 3];
    f.facets.insert(f.facets.end(), {a, b});
    f.markers.push_back(s.coords[a*2 + 1] == 0. && s.coords[b*2 + 1] == 0. ? 1 : 2);
    const double l = std::hypot(s.coords[a*2] - s.coords[b*2], s.coords[a*2 + 1] - s.coords[b*2 + 1]);
    h[a] = std::min(h[a], l);
    h[b] = std::min(h[b], l);
  }
  std::mt19937 rng(11);
  std::uniform_real_distribution<double> angle(0., 2.*M_PI);
  for (std::size_t v = 0; v < nverts; v += 2){
    if (!std::isfinite(h[v])) continue;
    const double a = angle(rng);
    f.points[v*2] += frac*h[v]*std::cos(a);
    f.points[v*2 + 1] += frac*h[v]*std::sin(a);
  }
  return f;
}

// The exterior facets the markers select, by cell
std::vector<std::uint8_t> strip_walls(const Strip& s, const bool bottom_only){
  std::vector<std::uint8_t> w(s.ncells, 0);
  for (std::size_t q = 0; q < s.ncells*3; ++q){
    if (s.across[q] != facet_wall) continue;
    const std::uint32_t* row = &s.topo[(q/3)*3];
    const bool bottom = s.coords[row[(q % 3 + 1) % 3]*2 + 1] == 0. && s.coords[row[(q % 3 + 2) % 3]*2 + 1] == 0.;
    if (!bottom_only || bottom) w[q/3] |= std::uint8_t(1u << (q % 3));
  }
  return w;
}

}  // namespace

TEST_CASE("wall_facets matches the points of a strongly graded mesh", "[wall_facets]"){
  // Boundary edges from 1e-8 to 0.25, a ratio above 1e7; half the points
  // moved by 0.4e-3 of the shortest boundary edge at their vertex; a file
  // without markers is refused
  const Strip s = graded_strip(20000, 1e-8, 1.2, 2);
  std::vector<double> h;
  const wall_facets::File f = strip_facets(s, 0.4e-3, h);
  double lo = 1e300, hi = 0.;
  for (const double l : h) if (std::isfinite(l)){ lo = std::min(lo, l); hi = std::max(hi, l); }
  REQUIRE(hi/lo >= 1e4);
  CHECK(wall_facets::match(f, wall_facets::Markers("1:2"), s.topo, s.coords, s.ncells, 2, s.across, "graded", false)
        == strip_walls(s, false));
  CHECK(wall_facets::match(f, wall_facets::Markers("1"), s.topo, s.coords, s.ncells, 2, s.across, "graded", false)
        == strip_walls(s, true));
  // a point moved past its tolerance, in the finest cells
  wall_facets::File g = f;
  const std::size_t v = 101;
  REQUIRE(h[v] < 2e-8);
  g.points[v*2] = s.coords[v*2] + 1.5e-3*h[v];
  g.points[v*2 + 1] = s.coords[v*2 + 1];
  std::string why;
  try { wall_facets::match(g, wall_facets::Markers("1:2"), s.topo, s.coords, s.ncells, 2, s.across, "graded", false); }
  catch (const partrac::Error& e){ why = e.what(); }
  CHECK(why.find("point " + std::to_string(v) + " (") != std::string::npos);
  CHECK(why.find("matches no vertex on the mesh's boundary") != std::string::npos);
  // a file without markers
  wall_facets::File bare = f;
  bare.markers.clear();
  why.clear();
  try { wall_facets::match(bare, wall_facets::Markers("1"), s.topo, s.coords, s.ncells, 2, s.across, "graded", false); }
  catch (const partrac::Error& e){ why = e.what(); }
  CHECK(why.find("graded holds no markers") != std::string::npos);
}

TEST_CASE("Mode openfoam, whose walls wall_facets cannot name, warns of no walls from the velocity",
          "[wall_facets][openfoam]"){
  if (!openfoam_load::available()){
    WARN("no OpenFOAM reader in this build");
    return;
  }
  const std::string data = std::string(PARTRAC_SOURCE_DIR) + "/data_example/";
  const std::string said = captured_cerr([&]{
    OpenFoamInterpol<Triangle> intp(data + "openfoam_cavity/partrac_params.dat");
    intp.update(intp.get_t_max());
  });
  CHECK(said.find("wall_facets") == std::string::npos);
}
