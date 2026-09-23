// The dolfin-free loader of XDMFInterpol, on XDMF cases dolfin writes here the
// way tests/conftest.py's xdmf_dir does (one xdmf per field, the mesh in the
// velocity's h5, one grid per stamp). A P1 field's dofs are the mesh vertices,
// so what the loader has to get right is silent in a run that only looks odd:
// the values belong to the vertices the topology names, a vertex on a periodic
// max face reads its master's value, and every stamp -- not only the first --
// reaches the cells on a periodic face. The near-wall P2 scheme is checked
// physically in tests/test_wall_p2.py; here both settings have to load and to
// differ only where a wall vertex is at rest.
#ifdef USE_DOLFIN
#include <catch2/catch.hpp>

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <dolfin.h>

#include "Error.hpp"
#include "h5direct.hpp"
#include "Tet.hpp"
#include "Triangle.hpp"
#include "XDMFTetInterpol.hpp"
#include "XDMFTriangleInterpol.hpp"
#include "case_dir.hpp"
#include "dolfin_ref.hpp"
#include "taylor_hood.hpp"

namespace {

// Stamp k of the velocity, of period one in every direction, so a vertex on a
// max face and its master carry one value in a constrained space
Vector3d stamp_u(const Vector3d& x, const int k, const Uint gdim){
  Vector3d u = Vector3d::Zero();
  for (Uint d = 0; d < gdim; ++d)
    u[d] = std::sin(2.*M_PI*x[(d + Uint(k)) % gdim] + 0.4*double(d)) + 0.5 + double(k);
  return u;
}

double stamp_p(const Vector3d& x, const int k){
  return std::cos(2.*M_PI*x[k == 0 ? 0 : 1] + 0.2) - 0.75*double(k);
}

class UExpr : public dolfin::Expression {
public:
  UExpr(const Uint gdim, const int k) : dolfin::Expression(gdim), gdim_(gdim), k_(k) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    Vector3d xx = Vector3d::Zero();
    for (Uint d = 0; d < gdim_; ++d) xx[d] = x[d];
    const Vector3d u = stamp_u(xx, k_, gdim_);
    for (Uint d = 0; d < gdim_; ++d) v[d] = u[d];
  }
private:
  Uint gdim_;
  int k_;
};

class PExpr : public dolfin::Expression {
public:
  PExpr(const Uint gdim, const int k) : gdim_(gdim), k_(k) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    Vector3d xx = Vector3d::Zero();
    for (Uint d = 0; d < gdim_; ++d) xx[d] = x[d];
    v[0] = stamp_p(xx, k_);
  }
private:
  Uint gdim_;
  int k_;
};

// A velocity that vanishes on the boundary of the unit box, so every wall
// vertex is at rest and the near-wall rule applies
class WallExpr : public dolfin::Expression {
public:
  explicit WallExpr(const Uint gdim) : dolfin::Expression(gdim), gdim_(gdim) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    double b = 1.;
    for (Uint d = 0; d < gdim_; ++d) b *= x[d]*(1. - x[d]);
    for (Uint d = 0; d < gdim_; ++d) v[d] = (d == 0 ? 4. : 1.)*b;
  }
private:
  Uint gdim_;
};

// A P1 case in XDMF: one xdmf per field, the mesh in the velocity's h5, one
// grid per stamp, as dolfin's XDMFFile writes it
template<typename Cell>
void write_xdmf_case(const CaseDir& c, const std::size_t n, const bool periodic,
                     const std::string& wall_p2, const int nstamps,
                     const bool at_rest_on_walls = false){
  constexpr int nv = Cell::n_verts;
  constexpr Uint gdim = nv - 1;
  std::shared_ptr<dolfin::Mesh> mesh;
  if constexpr (nv == 3) mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitSquareMesh(n, n));
  else                   mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitCubeMesh(n, n, n));
  const std::vector<bool> per(3, periodic);
  std::shared_ptr<const dolfin::SubDomain> pbc;
  if (periodic)
    pbc = std::make_shared<PeriodicBC>(per, Vector3d::Zero(), Vector3d(1., 1., 1.), gdim);
  std::shared_ptr<dolfin::FunctionSpace> V, P;
  Uint ncoeffs_u = 0, ncoeffs_p = 0;
  taylor_hood_spaces<Cell>("P1", "P1", true, mesh, pbc, V, P, ncoeffs_u, ncoeffs_p);

  dolfin::XDMFFile xf_u(mesh->mpi_comm(), (c.path / "u.xdmf").string());
  dolfin::XDMFFile xf_p(mesh->mpi_comm(), (c.path / "p.xdmf").string());
  for (dolfin::XDMFFile* f : {&xf_u, &xf_p}){
    // one mesh in the h5, referenced by every stamp's grid
    f->parameters["functions_share_mesh"] = true;
    f->parameters["rewrite_function_mesh"] = false;
  }
  for (int k = 0; k < nstamps; ++k){
    dolfin::Function u(V), p(P);
    if (at_rest_on_walls){
      const WallExpr ue(gdim);
      u.interpolate(ue);
    }
    else {
      const UExpr ue(gdim, k);
      u.interpolate(ue);
    }
    const PExpr pe(gdim, k);
    p.interpolate(pe);
    xf_u.write(u, double(k));
    xf_p.write(p, double(k));
  }
  xf_u.close();
  xf_p.close();

  std::ofstream(c.params())
    << "u=u.xdmf\np=p.xdmf\n"
    << "periodic_x=" << (periodic ? "true" : "false")
    << "\nperiodic_y=" << (periodic ? "true" : "false")
    << "\nperiodic_z=" << (periodic && gdim == 3 ? "true" : "false")
    << "\nwall_p2=" << wall_p2 << "\n";
}

// The cell centroids of the case's mesh, and how many sit in a cell with a
// vertex on a max face: those cells read an image vertex
template<int NV>
std::vector<Vector3d> centroids(const CaseDir& c, std::size_t& touching){
  constexpr Uint gdim = NV - 1;
  const partrac::H5Id file = partrac::h5_open_read((c.path / "u.h5").string());
  std::vector<std::uint32_t> topo;
  std::vector<double> coords;
  partrac::h5_read(file, "Mesh/0/mesh/topology", topo);
  partrac::h5_read(file, "Mesh/0/mesh/geometry", coords, gdim);
  const std::size_t ncells = topo.size()/NV;
  std::vector<Vector3d> out(ncells, Vector3d::Zero());
  touching = 0;
  for (std::size_t i = 0; i < ncells; ++i){
    bool on_face = false;
    for (int k = 0; k < NV; ++k){
      const double* v = coords.data() + std::size_t(topo[i*NV + k])*gdim;
      for (Uint d = 0; d < gdim; ++d){
        out[i][d] += v[d]/double(NV);
        if (v[d] > 1. - 1e-12) on_face = true;
      }
    }
    if (on_face) ++touching;
  }
  return out;
}

// The loader's velocity at every cell centroid, against the P1 value built
// here from the stored vertex values and the cell's own basis: bit for bit,
// since both read the same numbers in the same order
template<typename Cell>
void check_against_stored(XDMFInterpol<Cell>& intp, const CaseDir& c, const int stamp){
  constexpr int nv = Cell::n_verts;
  constexpr Uint gdim = nv - 1;
  const partrac::H5Id file = partrac::h5_open_read((c.path / "u.h5").string());
  std::vector<std::uint32_t> topo;
  std::vector<double> coords, vals;
  partrac::h5_read(file, "Mesh/0/mesh/topology", topo);
  partrac::h5_read(file, "Mesh/0/mesh/geometry", coords, gdim);
  // dolfin names a stamp's values by the order it wrote them
  partrac::h5_read(file, "VisualisationVector/" + std::to_string(stamp), vals, gdim);
  const std::size_t ncells = topo.size()/nv;

  intp.update(double(stamp));
  for (std::size_t i = 0; i < ncells; ++i){
    const std::uint32_t* row = topo.data() + i*nv;
    const double* v[4] = {nullptr, nullptr, nullptr, nullptr};
    for (int k = 0; k < nv; ++k) v[k] = coords.data() + std::size_t(row[k])*gdim;
    Cell cell;
    if constexpr (nv == 4) cell = Cell(v[0], v[1], v[2], v[3]);
    else                   cell = Cell(v[0], v[1], v[2]);
    Vector3d x = Vector3d::Zero();
    for (int k = 0; k < nv; ++k)
      for (Uint d = 0; d < gdim; ++d) x[d] += v[k][d]/double(nv);
    std::array<double, 4> bary;
    cell.contains(x, bary);
    std::array<double, Cell::n_dofs_max> N;
    if constexpr (nv == 3) cell.linearbasis(bary[0], bary[1], bary[2], N.data());
    else                   cell.linearbasis(bary[0], bary[1], bary[2], bary[3], N.data());
    Vector3d want = Vector3d::Zero();
    for (Uint d = 0; d < gdim; ++d)
      for (int k = 0; k < nv; ++k)
        want[d] += N[k]*vals[std::size_t(row[k])*gdim + d];

    CellPos pos;
    REQUIRE(intp.locate(x, double(stamp), pos));
    PointValues got(1.);
    intp.evaluate(x, double(stamp), pos, got);
    for (Uint d = 0; d < gdim; ++d) REQUIRE(got.get_u()[d] == want[d]);
  }
}

// The loader's velocity and pressure at the points, against the expressions
// the stamps were written from; the field is linear in time between them
template<typename Cell>
void check_two_stamps(XDMFInterpol<Cell>& intp, const std::vector<Vector3d>& pts,
                      const double tol_u, const double tol_p){
  constexpr Uint gdim = Cell::n_verts - 1;
  for (const double t : {0., 0.4, 1.}){
    intp.update(t);
    std::size_t inside = 0;
    for (const Vector3d& x : pts){
      CellPos pos;
      if (!intp.locate(x, t, pos)) continue;
      PointValues got(1.);
      intp.evaluate(x, t, pos, got);
      const Vector3d want = (1. - t)*stamp_u(x, 0, gdim) + t*stamp_u(x, 1, gdim);
      for (Uint d = 0; d < gdim; ++d)
        REQUIRE(std::abs(got.get_u()[d] - want[d]) < tol_u);
      REQUIRE(std::abs(got.get_p() - ((1. - t)*stamp_p(x, 0) + t*stamp_p(x, 1))) < tol_p);
      ++inside;
    }
    REQUIRE(inside == pts.size());
  }
}

}  // namespace

TEST_CASE("The XDMF loader's values are the stored vertex values", "[xdmf_load]") {
  // The P1 dofs are the vertices, so a cell reads the rows of its own three or
  // four of them; wall_p2=none keeps the evaluation P1 everywhere
  SECTION("triangles"){
    CaseDir c("tri_values");
    write_xdmf_case<Triangle>(c, 8, false, "none", 2);
    XDMFTriangleInterpol intp(c.params());
    check_against_stored(intp, c, 0);
    check_against_stored(intp, c, 1);
  }
  SECTION("tets"){
    CaseDir c("tet_values");
    write_xdmf_case<Tet>(c, 4, false, "none", 2);
    XDMFTetInterpol intp(c.params());
    check_against_stored(intp, c, 0);
  }
}

TEST_CASE("Both near-wall settings load, and differ only at a wall at rest", "[xdmf_load]") {
  // wall_p2=edge replaces the P1 velocity of a cell whose wall vertices are at
  // rest by a quadratic one; away from such a cell the two settings agree
  SECTION("a field that is nowhere at rest reads the same either way"){
    CaseDir c("tri_nowall");
    write_xdmf_case<Triangle>(c, 8, false, "edge", 1);
    CaseDir n("tri_nowall_none");
    write_xdmf_case<Triangle>(n, 8, false, "none", 1);
    std::size_t touching = 0;
    const std::vector<Vector3d> pts = centroids<3>(c, touching);
    XDMFTriangleInterpol edge(c.params());
    XDMFTriangleInterpol none(n.params());
    edge.update(0.);
    none.update(0.);
    for (const Vector3d& x : pts){
      CellPos pe, pn;
      REQUIRE(edge.locate(x, 0., pe));
      REQUIRE(none.locate(x, 0., pn));
      PointValues ge(1.), gn(1.);
      edge.evaluate(x, 0., pe, ge);
      none.evaluate(x, 0., pn, gn);
      REQUIRE(ge.get_u()[0] == gn.get_u()[0]);
      REQUIRE(ge.get_u()[1] == gn.get_u()[1]);
    }
  }
  SECTION("a field at rest on the walls is quadratic there and P1 inside"){
    CaseDir c("tri_wall");
    write_xdmf_case<Triangle>(c, 8, false, "edge", 1, true);
    CaseDir n("tri_wall_none");
    write_xdmf_case<Triangle>(n, 8, false, "none", 1, true);
    XDMFTriangleInterpol edge(c.params());
    XDMFTriangleInterpol none(n.params());
    edge.update(0.);
    none.update(0.);
    const auto value = [](XDMFTriangleInterpol& intp, const Vector3d& x){
      CellPos pos;
      REQUIRE(intp.locate(x, 0., pos));
      PointValues got(1.);
      intp.evaluate(x, 0., pos, got);
      return Vector3d(got.get_u());
    };
    // a point in the first row of cells above the wall y = 0, and one in the
    // middle; off the mesh lines, since on the edge rising from a wall vertex
    // the rule leaves the tangential component linear
    const Vector3d near_wall(0.53, 0.02, 0.), inside(0.53, 0.5, 0.);
    REQUIRE(value(edge, near_wall)[0] != value(none, near_wall)[0]);
    REQUIRE(value(edge, inside)[0] == value(none, inside)[0]);
  }
}

TEST_CASE("A later XDMF stamp reaches the cells on a periodic face", "[xdmf_load]") {
  // A vertex on a max face reads its master's value, and every stamp is read
  // by vertex, so a cell touching a periodic face holds the current stamp and
  // not the first one. The tolerance is the P1 interpolation error of a sin of
  // period one on this mesh; a cell left on the wrong stamp is off by order one.
  SECTION("triangles"){
    CaseDir c("per_tri");
    write_xdmf_case<Triangle>(c, 16, true, "edge", 2);
    std::size_t touching = 0;
    const std::vector<Vector3d> pts = centroids<3>(c, touching);
    REQUIRE(touching > 0);
    XDMFTriangleInterpol intp(c.params());
    check_two_stamps(intp, pts, 5e-2, 5e-2);
  }
  SECTION("tets"){
    CaseDir c("per_tet");
    write_xdmf_case<Tet>(c, 8, true, "edge", 2);
    std::size_t touching = 0;
    const std::vector<Vector3d> pts = centroids<4>(c, touching);
    REQUIRE(touching > 0);
    XDMFTetInterpol intp(c.params());
    check_two_stamps(intp, pts, 1.5e-1, 1.5e-1);
  }
}

#endif
