// The dolfin-free loader of XDMFInterpol, on XDMF cases dolfin writes here the
// way tests/conftest.py's xdmf_dir does (one xdmf per field, the mesh in the
// velocity's h5, one grid per stamp). A P1 field's dofs are the mesh vertices,
// so what the loader has to get right is silent in a run that only looks odd:
// the values belong to the vertices the topology names, a vertex on a periodic
// max face reads its master's value, and every stamp -- not only the first --
// reaches the cells on a periodic face. The near-wall P2 scheme is checked
// physically in tests/test_wall_edge_rule.py; here both settings have to load and to
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
#include "interpolators.hpp"
#include "case_dir.hpp"
#include "dolfin_cases.hpp"
#include "dolfin_ref.hpp"
#include "taylor_hood.hpp"

namespace {

// Stamp k lifted by 1/2 + k
const Stamps stamps{0.5, 1.};

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
                     const bool wall_upgrade, const int nstamps,
                     const bool at_rest_on_walls = false){
  constexpr int nv = Cell::n_verts;
  constexpr Uint gdim = nv - 1;
  const auto mesh = unit_box<Cell>(n);
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
      const StampU ue(stamps, gdim, k);
      u.interpolate(ue);
    }
    const StampP pe(stamps, gdim, k);
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
    << "\nwall_upgrade=" << (wall_upgrade ? "true" : "false") << "\n";
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
  check_centroids<false>(intp, topo, coords, double(stamp), [&](const std::size_t i, const auto& N){
    const std::uint32_t* row = topo.data() + i*nv;
    Vector3d want = Vector3d::Zero();
    for (Uint d = 0; d < gdim; ++d)
      for (int k = 0; k < nv; ++k)
        want[d] += N[k]*vals[std::size_t(row[k])*gdim + d];
    return want;
  });
}

}  // namespace

TEST_CASE("The XDMF loader's values are the stored vertex values", "[xdmf_load]") {
  // The P1 dofs are the vertices, so a cell reads the rows of its own three or
  // four of them; wall_upgrade=false keeps the evaluation P1 everywhere
  SECTION("triangles"){
    CaseDir c("tri_values");
    write_xdmf_case<Triangle>(c, 8, false, false, 2);
    XDMFTriangleInterpol intp(c.params());
    check_against_stored(intp, c, 0);
    check_against_stored(intp, c, 1);
  }
  SECTION("tets"){
    CaseDir c("tet_values");
    write_xdmf_case<Tet>(c, 4, false, false, 2);
    XDMFTetInterpol intp(c.params());
    check_against_stored(intp, c, 0);
  }
}

TEST_CASE("Both near-wall settings load, and differ only at a wall at rest", "[xdmf_load]") {
  // wall_upgrade=true replaces the P1 velocity of a cell whose wall vertices are at
  // rest by a quadratic one; away from such a cell the two settings agree
  SECTION("a field that is nowhere at rest reads the same either way"){
    CaseDir c("tri_nowall");
    write_xdmf_case<Triangle>(c, 8, false, true, 1);
    CaseDir n("tri_nowall_none");
    write_xdmf_case<Triangle>(n, 8, false, false, 1);
    std::size_t touching = 0;
    const std::vector<Vector3d> pts = centroids<3>(c, touching, true);
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
    write_xdmf_case<Triangle>(c, 8, false, true, 1, true);
    CaseDir n("tri_wall_none");
    write_xdmf_case<Triangle>(n, 8, false, false, 1, true);
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
    write_xdmf_case<Triangle>(c, 16, true, true, 2);
    std::size_t touching = 0;
    const std::vector<Vector3d> pts = centroids<3>(c, touching, true);
    REQUIRE(touching > 0);
    XDMFTriangleInterpol intp(c.params());
    check_two_stamps(intp, stamps, pts, 5e-2, 5e-2);
  }
  SECTION("tets"){
    CaseDir c("per_tet");
    write_xdmf_case<Tet>(c, 8, true, true, 2);
    std::size_t touching = 0;
    const std::vector<Vector3d> pts = centroids<4>(c, touching, true);
    REQUIRE(touching > 0);
    XDMFTetInterpol intp(c.params());
    check_two_stamps(intp, stamps, pts, 1.5e-1, 1.5e-1);
  }
}

#endif
