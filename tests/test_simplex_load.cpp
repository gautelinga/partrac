// The dolfin-free loader of SimplexInterpol, on fixtures dolfin itself writes
// here (a mesh and a checkpointed Taylor-Hood pair, as data_example's
// generate_up.py does), so the arrays under test are the arrays dolfin stores.
// What it has to get right is silent when a run only looks slow or odd: the
// composition of the field file's rows with the mesh's cell ids, the edge
// numbering that gives a P2 cell its six midside nodes in the basis's order,
// and the refusal to read a file it does not understand. The checks against
// dolfin's own evaluation live in test_interpol_core.cpp; here the reference is
// the stored dof table read a second time, in this file, by hand.
#ifdef USE_DOLFIN
#include <catch2/catch.hpp>

#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include <dolfin.h>
#include <hdf5.h>

#include "Error.hpp"
#include "h5direct.hpp"
#include "Tet.hpp"
#include "Triangle.hpp"
#include "TetInterpol.hpp"
#include "TriangleInterpol.hpp"
#include "case_dir.hpp"
#include "dolfin_ref.hpp"
#include "taylor_hood.hpp"

namespace {

// A unit mesh with a P2 velocity and a P1 pressure, written the way the solvers
// write them. The values are arbitrary but deterministic: the loader is judged
// against the numbers in the file, not against a physical field.
template<typename Cell>
void write_case(const CaseDir& c, const std::size_t n){
  std::shared_ptr<dolfin::Mesh> mesh;
  if constexpr (Cell::n_verts == 3) mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitSquareMesh(n, n));
  else                              mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitCubeMesh(n, n, n));
  std::shared_ptr<dolfin::FunctionSpace> V, P;
  Uint ncoeffs_u = 0, ncoeffs_p = 0;
  taylor_hood_spaces<Cell>("P2", "P1", true, mesh, nullptr, V, P, ncoeffs_u, ncoeffs_p);

  dolfin::Function u(V), p(P);
  std::vector<double> uv(u.vector()->local_size()), pv(p.vector()->local_size());
  for (std::size_t i = 0; i < uv.size(); ++i) uv[i] = std::sin(0.7*double(i) + 0.3);
  for (std::size_t i = 0; i < pv.size(); ++i) pv[i] = std::cos(1.1*double(i) + 0.2);
  u.vector()->set_local(uv);
  u.vector()->apply("insert");
  p.vector()->set_local(pv);
  p.vector()->apply("insert");

  {
    dolfin::HDF5File f(MPI_COMM_WORLD, c.file("mesh.h5"), "w");
    f.write(*mesh, "mesh");
  }
  {
    dolfin::HDF5File f(MPI_COMM_WORLD, c.file("up_0.h5"), "w");
    f.write(u, "u");
    f.write(p, "p");
  }
  std::ofstream(c.path / "timestamps.dat") << "0\tup_0.h5\n1\tup_0.h5\n";
  std::ofstream(c.params())
    << "velocity_space=P2\npressure_space=P1\ntimestamps=timestamps.dat\nmesh=mesh.h5\n"
    << "periodic_x=false\nperiodic_y=false\nperiodic_z=false\n";
}

// The stored dof table of a field, one row per mesh row, as the loader composes it
std::vector<std::uint32_t> stored_rows(const CaseDir& c, const std::string& field,
                                       const std::size_t ncells, std::size_t& per_cell,
                                       std::vector<double>& vec){
  const partrac::H5Id mesh = partrac::h5_open_read(c.file("mesh.h5"));
  std::vector<std::uint64_t> gid;
  partrac::h5_read(mesh, "mesh/cell_indices", gid);
  const partrac::H5Id file = partrac::h5_open_read(c.file("up_0.h5"));
  std::vector<std::uint32_t> flat;
  std::vector<std::uint64_t> cells;
  partrac::h5_read(file, field + "/cell_dofs", flat);
  partrac::h5_read(file, field + "/cells", cells);
  partrac::h5_read(file, field + "/vector_0", vec);
  per_cell = flat.size()/ncells;
  std::vector<std::uint32_t> where(ncells, 0);
  for (std::size_t i = 0; i < ncells; ++i) where[std::size_t(gid[i])] = std::uint32_t(i);
  std::vector<std::uint32_t> rows(flat.size());
  for (std::size_t i = 0; i < ncells; ++i)
    std::copy(flat.begin() + std::ptrdiff_t(i*per_cell), flat.begin() + std::ptrdiff_t((i+1)*per_cell),
              rows.begin() + std::ptrdiff_t(std::size_t(where[std::size_t(cells[i])])*per_cell));
  return rows;
}

// Overwrites the first element of a dataset with value
template<typename T>
void poke(const std::string& path, const std::string& dataset, const T value,
          const hid_t type, const hsize_t at = 0){
  const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  REQUIRE(file >= 0);
  const hid_t dset = H5Dopen2(file, dataset.c_str(), H5P_DEFAULT);
  REQUIRE(dset >= 0);
  const hid_t space = H5Dget_space(dset);
  const hsize_t one = 1;
  H5Sselect_hyperslab(space, H5S_SELECT_SET, &at, nullptr, &one, nullptr);
  const hid_t mem = H5Screate_simple(1, &one, nullptr);
  REQUIRE(H5Dwrite(dset, type, mem, space, H5P_DEFAULT, &value) >= 0);
  H5Sclose(mem);
  H5Sclose(space);
  H5Dclose(dset);
  H5Fclose(file);
}

// Replaces a group's signature attribute
void set_signature(const std::string& path, const std::string& group, const std::string& value){
  const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  REQUIRE(file >= 0);
  const hid_t g = H5Gopen2(file, group.c_str(), H5P_DEFAULT);
  REQUIRE(g >= 0);
  H5Adelete(g, "signature");
  const hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, value.size() + 1);
  const hid_t space = H5Screate(H5S_SCALAR);
  const hid_t attr = H5Acreate2(g, "signature", type, space, H5P_DEFAULT, H5P_DEFAULT);
  REQUIRE(attr >= 0);
  REQUIRE(H5Awrite(attr, type, value.c_str()) >= 0);
  H5Aclose(attr);
  H5Sclose(space);
  H5Tclose(type);
  H5Gclose(g);
  H5Fclose(file);
}

// The loader's velocity at x, and the same value built here from the stored
// dof table and the cell's own basis: the two must agree to the last bit
template<typename Cell>
void check_against_stored(SimplexInterpol<Cell>& intp, const CaseDir& c){
  constexpr int nv = Cell::n_verts;
  constexpr Uint gdim = nv - 1;
  const partrac::H5Id mesh = partrac::h5_open_read(c.file("mesh.h5"));
  std::vector<std::uint32_t> topo;
  std::vector<double> coords;
  partrac::h5_read(mesh, "mesh/topology", topo);
  partrac::h5_read(mesh, "mesh/coordinates", coords, gdim);
  const std::size_t ncells = topo.size()/nv;
  std::size_t per_cell = 0;
  std::vector<double> vec;
  const std::vector<std::uint32_t> rows = stored_rows(c, "u", ncells, per_cell, vec);
  const std::size_t n_nodes = per_cell/gdim;

  intp.update(0.);
  std::size_t checked = 0;
  for (std::size_t i = 0; i < ncells; ++i){
    const std::uint32_t* row = topo.data() + i*nv;
    const double* v[4] = {nullptr, nullptr, nullptr, nullptr};
    for (int k = 0; k < nv; ++k) v[k] = coords.data() + std::size_t(row[k])*gdim;
    Cell cell;
    if constexpr (nv == 4) cell = Cell(v[0], v[1], v[2], v[3]);
    else                   cell = Cell(v[0], v[1], v[2]);
    // the centroid is interior, so one cell holds it and the walk cannot differ
    Vector3d x = Vector3d::Zero();
    for (int k = 0; k < nv; ++k)
      for (Uint d = 0; d < gdim; ++d) x[d] += v[k][d]/double(nv);
    std::array<double, 4> bary;
    cell.contains(x, bary);
    std::array<double, Cell::n_dofs_max> N;
    if constexpr (nv == 3) cell.quadbasis(bary[0], bary[1], bary[2], N.data());
    else                   cell.quadbasis(bary[0], bary[1], bary[2], bary[3], N.data());
    Vector3d want = Vector3d::Zero();
    for (std::size_t comp = 0; comp < gdim; ++comp)
      for (std::size_t k = 0; k < n_nodes; ++k)
        want[comp] += N[k]*vec[rows[i*per_cell + comp*n_nodes + k]];

    CellPos pos;
    REQUIRE(intp.locate(x, 0., pos));
    PointValues got(1.);
    intp.evaluate(x, 0., pos, got);
    for (Uint d = 0; d < gdim; ++d) REQUIRE(got.get_u()[d] == want[d]);
    ++checked;
  }
  REQUIRE(checked == ncells);
}

// The field of stamp k, periodic with period one in every direction of the unit
// box, so a node on a max face and its master on the min face carry one value
// and the constrained space stores one dof for the pair
Vector3d stamp_u(const Vector3d& x, const int k, const Uint gdim){
  Vector3d u = Vector3d::Zero();
  for (Uint d = 0; d < gdim; ++d)
    u[d] = std::sin(2.*M_PI*x[(d + Uint(k)) % gdim] + 0.4*double(d)) + 0.5*double(k);
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

// The parameter file of a periodic case: a reload is given its periodicity and
// its cache here, and the files it names stay as they are
void write_params(const CaseDir& c, const std::string& u_el, const std::vector<bool>& periodic,
                  const bool cache){
  std::ofstream(c.params())
    << "velocity_space=" << u_el << "\npressure_space=P1\ntimestamps=timestamps.dat\nmesh=mesh.h5\n"
    << "periodic_x=" << (periodic[0] ? "true" : "false")
    << "\nperiodic_y=" << (periodic[1] ? "true" : "false")
    << "\nperiodic_z=" << (periodic[2] ? "true" : "false")
    << "\nmesh_cache=" << (cache ? "true" : "false") << "\n";
}

// A periodic case with two stamps, written from a constrained space, as the
// solvers write a periodic run: up_0.h5 and up_1.h5 at t = 0 and t = 1
template<typename Cell>
void write_periodic_case(const CaseDir& c, const std::size_t n, const std::string& u_el,
                         const bool cache){
  constexpr int nv = Cell::n_verts;
  constexpr Uint gdim = nv - 1;
  std::shared_ptr<dolfin::Mesh> mesh;
  if constexpr (nv == 3) mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitSquareMesh(n, n));
  else                   mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitCubeMesh(n, n, n));
  const std::vector<bool> per(3, true);
  const std::shared_ptr<const dolfin::SubDomain> pbc =
    std::make_shared<PeriodicBC>(per, Vector3d::Zero(), Vector3d(1., 1., 1.), gdim);
  std::shared_ptr<dolfin::FunctionSpace> V, P;
  Uint ncoeffs_u = 0, ncoeffs_p = 0;
  taylor_hood_spaces<Cell>(u_el, "P1", true, mesh, pbc, V, P, ncoeffs_u, ncoeffs_p);

  {
    dolfin::HDF5File f(MPI_COMM_WORLD, c.file("mesh.h5"), "w");
    f.write(*mesh, "mesh");
  }
  for (int k = 0; k < 2; ++k){
    dolfin::Function u(V), p(P);
    const UExpr ue(gdim, k);
    const PExpr pe(gdim, k);
    u.interpolate(ue);
    p.interpolate(pe);
    dolfin::HDF5File f(MPI_COMM_WORLD, (c.path / ("up_" + std::to_string(k) + ".h5")).string(), "w");
    f.write(u, "u");
    f.write(p, "p");
  }
  std::ofstream(c.path / "timestamps.dat") << "0\tup_0.h5\n1\tup_1.h5\n";
  write_params(c, u_el, {true, true, true}, cache);
}

// The periodic field of stamp 0, with every node on a max face moved off its
// master by nudge: a file an unconstrained space wrote agrees across the seam
// only as closely as the values in it do
class SeamExpr : public dolfin::Expression {
public:
  SeamExpr(const Uint gdim, const double nudge)
    : dolfin::Expression(gdim), gdim_(gdim), nudge_(nudge) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    Vector3d xx = Vector3d::Zero();
    bool on_max = false;
    for (Uint d = 0; d < gdim_; ++d){
      xx[d] = x[d];
      if (x[d] > 1. - 1e-12) on_max = true;
    }
    const Vector3d u = stamp_u(xx, 0, gdim_);
    for (Uint d = 0; d < gdim_; ++d) v[d] = u[d] + (on_max ? nudge_ : 0.);
  }
private:
  Uint gdim_;
  double nudge_;
};

// One stamp of the same periodic field, written from a constrained space --
// where a node and its images share a dof -- or from an unconstrained one,
// where they are independent and nudge separates their values
template<typename Cell>
void write_seam_case(const CaseDir& c, const std::size_t n, const std::string& u_el,
                     const bool constrained, const double nudge){
  constexpr int nv = Cell::n_verts;
  constexpr Uint gdim = nv - 1;
  std::shared_ptr<dolfin::Mesh> mesh;
  if constexpr (nv == 3) mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitSquareMesh(n, n));
  else                   mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitCubeMesh(n, n, n));
  std::shared_ptr<const dolfin::SubDomain> pbc;
  if (constrained)
    pbc = std::make_shared<PeriodicBC>(std::vector<bool>(3, true), Vector3d::Zero(),
                                       Vector3d(1., 1., 1.), gdim);
  std::shared_ptr<dolfin::FunctionSpace> V, P;
  Uint ncoeffs_u = 0, ncoeffs_p = 0;
  taylor_hood_spaces<Cell>(u_el, "P1", true, mesh, pbc, V, P, ncoeffs_u, ncoeffs_p);
  {
    dolfin::HDF5File f(MPI_COMM_WORLD, c.file("mesh.h5"), "w");
    f.write(*mesh, "mesh");
  }
  dolfin::Function u(V), p(P);
  const SeamExpr ue(gdim, nudge);
  const PExpr pe(gdim, 0);
  u.interpolate(ue);
  p.interpolate(pe);
  {
    dolfin::HDF5File f(MPI_COMM_WORLD, c.file("up_0.h5"), "w");
    f.write(u, "u");
    f.write(p, "p");
  }
  std::ofstream(c.path / "timestamps.dat") << "0\tup_0.h5\n1\tup_0.h5\n";
  write_params(c, u_el, {true, true, true}, false);
}

// The centroid of every cell, and how many of them sit in a cell with a vertex
// on a periodic face: those cells read an image node, whose dof is shared with
// a master somewhere else in the box
template<int NV>
std::vector<Vector3d> centroids(const CaseDir& c, std::size_t& touching){
  constexpr Uint gdim = NV - 1;
  const partrac::H5Id mesh = partrac::h5_open_read(c.file("mesh.h5"));
  std::vector<std::uint32_t> topo;
  std::vector<double> coords;
  partrac::h5_read(mesh, "mesh/topology", topo);
  partrac::h5_read(mesh, "mesh/coordinates", coords, gdim);
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

// The loader's velocity and pressure against the expressions the two stamps
// were written from, at the first stamp, between them and at the second. The
// field is linear in time between stamps, so the reference is too. A node whose
// dof is shared with a master the loader wrote instead would hold the wrong
// stamp from the first update on, by far more than the interpolation error.
template<typename Cell>
void check_two_stamps(SimplexInterpol<Cell>& intp, const std::vector<Vector3d>& pts,
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

// The wall normal of every cell, which is the facet table seen from outside: a
// facet across a periodic seam is no wall, so these say which periodicity the
// table was built with
template<typename Cell>
std::vector<Vector3d> wall_normals(SimplexInterpol<Cell>& intp, const std::size_t ncells){
  std::vector<Vector3d> out;
  for (std::size_t i = 0; i < ncells; ++i){
    int id = int(i);
    out.push_back(intp.get_boundary_normal(Vector3d::Zero(), id));
  }
  return out;
}

// Every value the loader gives at these points and times, to compare a cached
// load with the load that wrote the cache
template<typename Cell>
std::vector<double> sample(SimplexInterpol<Cell>& intp, const std::vector<Vector3d>& pts){
  constexpr Uint gdim = Cell::n_verts - 1;
  std::vector<double> out;
  for (const double t : {0., 0.4, 1.}){
    intp.update(t);
    for (const Vector3d& x : pts){
      CellPos pos;
      REQUIRE(intp.locate(x, t, pos));
      PointValues got(1.);
      intp.evaluate(x, t, pos, got);
      for (Uint d = 0; d < gdim; ++d) out.push_back(got.get_u()[d]);
      out.push_back(got.get_p());
    }
  }
  return out;
}

}  // namespace

TEST_CASE("The loader's values are the stored dof table's", "[simplex_load]") {
  SECTION("tets"){
    CaseDir c("tet_values");
    write_case<Tet>(c, 3);
    TetInterpol intp(c.params());
    check_against_stored(intp, c);
  }
  SECTION("triangles"){
    CaseDir c("tri_values");
    write_case<Triangle>(c, 4);
    TriangleInterpol intp(c.params());
    check_against_stored(intp, c);
  }
}

TEST_CASE("A field with one stamp holds one copy of it", "[simplex_load]") {
  CaseDir c("single");
  write_case<Tet>(c, 3);
  TetInterpol intp(c.params());
  intp.update(0.);
  REQUIRE(intp.stamps_aliased());
  // and the field is then steady: no acceleration between two equal stamps
  CellPos pos;
  const Vector3d x(0.31, 0.42, 0.53);
  REQUIRE(intp.locate(x, 0., pos));
  PointValues vals(1.);
  intp.evaluate(x, 0., pos, vals);
  REQUIRE(vals.get_a() == Vector3d::Zero());
}

TEST_CASE("The loader refuses a file it cannot read as dolfin wrote it", "[simplex_load]") {
  SECTION("an element signature it does not know"){
    CaseDir c("bad_sig");
    write_case<Tet>(c, 2);
    set_signature(c.file("up_0.h5"), "u", "VectorElement(FiniteElement('Lagrange', tetrahedron, 3), dim=3)");
    REQUIRE_THROWS_AS(TetInterpol(c.params()), partrac::Error);
  }
  SECTION("an element on the wrong cell"){
    CaseDir c("bad_cell");
    write_case<Tet>(c, 2);
    set_signature(c.file("up_0.h5"), "u", "VectorElement(FiniteElement('Lagrange', triangle, 2), dim=3)");
    REQUIRE_THROWS_AS(TetInterpol(c.params()), partrac::Error);
  }
  SECTION("a ragged x_cell_dofs, which cannot be reshaped"){
    CaseDir c("ragged");
    write_case<Tet>(c, 2);
    poke<std::uint64_t>(c.file("up_0.h5"), "u/x_cell_dofs", 7, H5T_NATIVE_UINT64, 1);
    REQUIRE_THROWS_AS(TetInterpol(c.params()), partrac::Error);
  }
  SECTION("a cell the mesh file does not have"){
    CaseDir c("no_cell");
    write_case<Tet>(c, 2);
    poke<std::uint64_t>(c.file("up_0.h5"), "u/cells", 1000000, H5T_NATIVE_UINT64);
    REQUIRE_THROWS_AS(TetInterpol(c.params()), partrac::Error);
  }
  SECTION("a global cell id wider than the tables"){
    CaseDir c("wide_id");
    write_case<Tet>(c, 2);
    poke<std::int64_t>(c.file("mesh.h5"), "mesh/cell_indices", 5000000000LL, H5T_NATIVE_INT64);
    REQUIRE_THROWS_AS(TetInterpol(c.params()), partrac::Error);
  }
  SECTION("two cells that disagree about a node they share"){
    CaseDir c("clash");
    write_case<Tet>(c, 2);
    // the first cell's first vertex, pointed at another node's dof: every other
    // cell holding that vertex then reads a different value for it
    std::size_t per_cell = 0;
    std::vector<double> vec;
    const partrac::H5Id mesh = partrac::h5_open_read(c.file("mesh.h5"));
    std::vector<std::uint32_t> topo;
    partrac::h5_read(mesh, "mesh/topology", topo);
    const std::vector<std::uint32_t> rows = stored_rows(c, "u", topo.size()/4, per_cell, vec);
    std::uint32_t other = rows[0];
    for (std::uint32_t d = 0; d < std::uint32_t(vec.size()); ++d)
      if (vec[d] != vec[rows[0]]){ other = d; break; }
    poke<std::int32_t>(c.file("up_0.h5"), "u/cell_dofs", std::int32_t(other), H5T_NATIVE_INT32);
    REQUIRE_THROWS_AS(TetInterpol(c.params()), partrac::Error);
  }
}

TEST_CASE("A later stamp reaches a periodic image node", "[simplex_load]") {
  // In a periodic-reduced space one stored dof serves a node and its images, so
  // the loader keeps a dof -> nodes mapping rather than a permutation. The
  // tolerances are the interpolation error of a sin of period one on this mesh
  // (P2 cubic, P1 quadratic in the cell size); a node left on the wrong stamp
  // is off by order one, which is what this is watching for.
  SECTION("tets, P2 velocity and P1 pressure"){
    CaseDir c("per_tet_p2");
    write_periodic_case<Tet>(c, 8, "P2", false);
    std::size_t touching = 0;
    const std::vector<Vector3d> pts = centroids<4>(c, touching);
    REQUIRE(touching > 0);
    TetInterpol intp(c.params());
    check_two_stamps(intp, pts, 1e-2, 1e-1);
  }
  SECTION("triangles, P2 velocity and P1 pressure"){
    CaseDir c("per_tri_p2");
    write_periodic_case<Triangle>(c, 12, "P2", false);
    std::size_t touching = 0;
    const std::vector<Vector3d> pts = centroids<3>(c, touching);
    REQUIRE(touching > 0);
    TriangleInterpol intp(c.params());
    check_two_stamps(intp, pts, 1e-2, 6e-2);
  }
  SECTION("tets, P1 velocity: the vertices are the whole space"){
    CaseDir c("per_tet_p1");
    write_periodic_case<Tet>(c, 8, "P1", false);
    std::size_t touching = 0;
    const std::vector<Vector3d> pts = centroids<4>(c, touching);
    REQUIRE(touching > 0);
    TetInterpol intp(c.params());
    check_two_stamps(intp, pts, 1e-1, 1e-1);
  }
}

TEST_CASE("The cache file gives a periodic two-stamp field back unchanged", "[simplex_load]") {
  // The cache holds the first stamp and the dof mapping, so a cached load has
  // to read the later stamp the same way the load that wrote it did
  CaseDir c("per_cache");
  write_periodic_case<Tet>(c, 8, "P2", true);
  std::size_t touching = 0;
  const std::vector<Vector3d> pts = centroids<4>(c, touching);
  REQUIRE(touching > 0);
  std::vector<double> fresh, cached;
  {
    TetInterpol intp(c.params());
    fresh = sample(intp, pts);
  }
  REQUIRE(std::filesystem::exists(c.path / "mesh_partrac_tet.h5"));
  {
    TetInterpol intp(c.params());
    cached = sample(intp, pts);
  }
  REQUIRE(cached == fresh);
  {
    TetInterpol intp(c.params());
    check_two_stamps(intp, pts, 1e-2, 1e-1);
  }
}

TEST_CASE("A cache written with one periodicity is not read with another", "[simplex_load]") {
  // The facet table is built from the periodic flags, so a cache keyed by the
  // input files alone would hand a run the neighbours of another periodicity
  CaseDir c("cache_periodic");
  write_periodic_case<Tet>(c, 4, "P1", true);
  std::size_t touching = 0;
  const std::size_t ncells = centroids<4>(c, touching).size();
  { TetInterpol intp(c.params()); }   // periodic in every direction: writes the cache
  REQUIRE(std::filesystem::exists(c.path / "mesh_partrac_tet.h5"));

  write_params(c, "P1", {true, true, false}, false);
  std::vector<Vector3d> fresh;
  { TetInterpol intp(c.params()); fresh = wall_normals(intp, ncells); }
  write_params(c, "P1", {true, true, false}, true);
  std::vector<Vector3d> cached;
  { TetInterpol intp(c.params()); cached = wall_normals(intp, ncells); }

  REQUIRE(cached.size() == fresh.size());
  for (std::size_t i = 0; i < fresh.size(); ++i)
    REQUIRE((cached[i] - fresh[i]).norm() == 0.);
  // the walls at z = 0 and z = 1 are what tells the two tables apart
  std::size_t walls = 0;
  for (const Vector3d& n : fresh) if (n.norm() > 0.) ++walls;
  REQUIRE(walls > 0);
}

TEST_CASE("A cache whose arrays do not fit its counts is rebuilt", "[simplex_load]") {
  // A cache the key still matches but that holds fewer values than its counts
  // say would be read past its end; the load falls back on the input files
  CaseDir c("cache_short");
  write_periodic_case<Tet>(c, 4, "P1", true);
  std::size_t touching = 0;
  const std::vector<Vector3d> pts = centroids<4>(c, touching);
  std::vector<double> fresh;
  { TetInterpol intp(c.params()); fresh = sample(intp, pts); }
  const std::string cache = (c.path / "mesh_partrac_tet.h5").string();
  REQUIRE(std::filesystem::exists(cache));
  // the vertex count, which the coordinates and the values are read against
  poke<double>(cache, "scalars", 7., H5T_NATIVE_DOUBLE, 1);
  TetInterpol intp(c.params());
  REQUIRE(sample(intp, pts) == fresh);
}

TEST_CASE("A periodic field written from an unconstrained space is read through its masters",
          "[simplex_load]") {
  // Nothing in such a file pairs a node with its image, so the loader pairs
  // them by position and a cell reads the master's value, not the image's own;
  // where the two sides disagree the master still serves both and the load says so.
  SECTION("triangles, P2 velocity: the values of the constrained file"){
    CaseDir a("seam_constrained"), b("seam_free");
    write_seam_case<Triangle>(a, 12, "P2", true, 0.);
    write_seam_case<Triangle>(b, 12, "P2", false, 1e-9);
    std::size_t touching = 0;
    const std::vector<Vector3d> pts = centroids<3>(a, touching);
    REQUIRE(touching > 0);
    TriangleInterpol ia(a.params()), ib(b.params());
    const std::vector<double> va = sample(ia, pts), vb = sample(ib, pts);
    REQUIRE(va.size() == vb.size());
    // the image nodes are 1e-9 off their masters: a cell reading one shows it
    for (std::size_t i = 0; i < va.size(); ++i)
      REQUIRE(std::abs(va[i] - vb[i]) < 1e-12);
  }
  SECTION("a field whose two sides disagree loads, and reads the masters"){
    CaseDir a("seam_ref"), c("seam_broken");
    write_seam_case<Triangle>(a, 6, "P2", true, 0.);
    write_seam_case<Triangle>(c, 6, "P2", false, 1e-3);
    std::size_t touching = 0;
    const std::vector<Vector3d> pts = centroids<3>(a, touching);
    REQUIRE(touching > 0);
    TriangleInterpol ia(a.params());
    TriangleInterpol ic(c.params());   // reported, not refused
    const std::vector<double> va = sample(ia, pts), vc = sample(ic, pts);
    REQUIRE(va.size() == vc.size());
    // the images are 1e-3 off and every cell reads the master side
    for (std::size_t i = 0; i < va.size(); ++i)
      REQUIRE(std::abs(va[i] - vc[i]) < 1e-12);
  }
}

#endif
