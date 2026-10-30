// What the tests of the stamped loaders share: dolfin's unit meshes and their
// files, the fields the stamps are written from, and checks against them
#pragma once
#ifdef USE_DOLFIN

#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <catch2/catch.hpp>
#include <dolfin.h>

#include "PointValues.hpp"
#include "StampedInterpol.hpp"
#include "case_dir.hpp"
#include "h5direct.hpp"
#include "typedefs.hpp"

namespace {

// dolfin's unit square (n x n) or cube (n x n x n)
template<typename Cell>
std::shared_ptr<dolfin::Mesh> unit_box(const std::size_t n){
  if constexpr (Cell::n_verts == 3) return std::make_shared<dolfin::Mesh>(dolfin::UnitSquareMesh(n, n));
  else                              return std::make_shared<dolfin::Mesh>(dolfin::UnitCubeMesh(n, n, n));
}

// The mesh as mesh.h5
inline void write_mesh(const CaseDir& c, const dolfin::Mesh& mesh){
  dolfin::HDF5File f(MPI_COMM_WORLD, c.file("mesh.h5"), "w");
  f.write(mesh, "mesh");
}

// The cell centroids of the case's mesh, the checkpoint's mesh.h5 or XDMF's
// u.h5, and how many sit in a cell with a vertex on a max face: those cells
// read an image node, whose dof is shared with a master elsewhere in the box
template<int NV>
std::vector<Vector3d> centroids(const CaseDir& c, std::size_t& touching, const bool xdmf = false){
  constexpr Uint gdim = NV - 1;
  const partrac::H5Id file = partrac::h5_open_read(c.file(xdmf ? "u.h5" : "mesh.h5"));
  std::vector<std::uint32_t> topo;
  std::vector<double> coords;
  partrac::h5_read(file, xdmf ? "Mesh/0/mesh/topology" : "mesh/topology", topo);
  partrac::h5_read(file, xdmf ? "Mesh/0/mesh/geometry" : "mesh/coordinates", coords, gdim);
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

// The loader's velocity at every cell's centroid at time t, against want(i,
// N), N the cell's quadratic or linear basis there: bit for bit
template<bool Quad, typename Cell, typename Format, typename Want>
void check_centroids(StampedInterpol<Cell, Format>& intp, const std::vector<std::uint32_t>& topo,
                     const std::vector<double>& coords, const double t, const Want& want){
  constexpr int nv = Cell::n_verts;
  constexpr Uint gdim = nv - 1;
  const std::size_t ncells = topo.size()/nv;
  intp.update(t);
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
    if constexpr (Quad && nv == 3)  cell.quadbasis(bary[0], bary[1], bary[2], N.data());
    else if constexpr (Quad)        cell.quadbasis(bary[0], bary[1], bary[2], bary[3], N.data());
    else if constexpr (nv == 3)     cell.linearbasis(bary[0], bary[1], bary[2], N.data());
    else                            cell.linearbasis(bary[0], bary[1], bary[2], bary[3], N.data());
    const Vector3d w = want(i, N);

    CellPos pos;
    REQUIRE(intp.locate(x, t, pos));
    PointValues got(1.);
    intp.evaluate(x, t, pos, got);
    for (Uint d = 0; d < gdim; ++d) REQUIRE(got.get_u()[d] == w[d]);
  }
}

// Stamp k of the velocity, of period one in every direction of the unit box,
// so a vertex on a max face and its master carry one value in a constrained
// space; lifted by base + step k
struct Stamps {
  double base, step;
  Vector3d u(const Vector3d& x, const int k, const Uint gdim) const {
    Vector3d u = Vector3d::Zero();
    for (Uint d = 0; d < gdim; ++d)
      u[d] = std::sin(2.*M_PI*x[(d + Uint(k)) % gdim] + 0.4*double(d)) + base + step*double(k);
    return u;
  }
  double p(const Vector3d& x, const int k) const {
    return std::cos(2.*M_PI*x[k == 0 ? 0 : 1] + 0.2) - 0.75*double(k);
  }
};

class StampU : public dolfin::Expression {
public:
  StampU(const Stamps& s, const Uint gdim, const int k) : dolfin::Expression(gdim), s_(s), gdim_(gdim), k_(k) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    Vector3d xx = Vector3d::Zero();
    for (Uint d = 0; d < gdim_; ++d) xx[d] = x[d];
    const Vector3d u = s_.u(xx, k_, gdim_);
    for (Uint d = 0; d < gdim_; ++d) v[d] = u[d];
  }
private:
  Stamps s_;
  Uint gdim_;
  int k_;
};

class StampP : public dolfin::Expression {
public:
  StampP(const Stamps& s, const Uint gdim, const int k) : s_(s), gdim_(gdim), k_(k) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    Vector3d xx = Vector3d::Zero();
    for (Uint d = 0; d < gdim_; ++d) xx[d] = x[d];
    v[0] = s_.p(xx, k_);
  }
private:
  Stamps s_;
  Uint gdim_;
  int k_;
};

// The loader's velocity and pressure at the points, against the stamps they
// were written from, at the first stamp, between them and at the second; the
// field is linear in time between stamps, so the reference is too. A node
// that holds the wrong stamp is off by the step, far more than the
// interpolation error
template<typename Cell, typename Format>
void check_two_stamps(StampedInterpol<Cell, Format>& intp, const Stamps& s, const std::vector<Vector3d>& pts,
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
      const Vector3d want = (1. - t)*s.u(x, 0, gdim) + t*s.u(x, 1, gdim);
      for (Uint d = 0; d < gdim; ++d)
        REQUIRE(std::abs(got.get_u()[d] - want[d]) < tol_u);
      REQUIRE(std::abs(got.get_p() - ((1. - t)*s.p(x, 0) + t*s.p(x, 1))) < tol_p);
      ++inside;
    }
    REQUIRE(inside == pts.size());
  }
}

// A quadratic, so a P2 space holds it exactly; lifted by k/2 at stamp k
inline double phi_exact(const Vector3d& x, const Uint gdim, const int k = 0){
  return x[0]*x[0] - x[0]*x[1] + 0.5*x[1]*x[1] + (gdim == 3 ? x[2]*x[2] - 0.3*x[2] : 0.)
       + 0.5*double(k);
}

class PhiExpr : public dolfin::Expression {
public:
  explicit PhiExpr(const Uint gdim, const int k = 0) : gdim_(gdim), k_(k) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    Vector3d xx = Vector3d::Zero();
    for (Uint d = 0; d < gdim_; ++d) xx[d] = x[d];
    v[0] = phi_exact(xx, gdim_, k_);
  }
private:
  Uint gdim_;
  int k_;
};

inline void append(const std::string& params, const std::string& lines){
  std::ofstream(params, std::ios::app) << lines;
}

// Points spread over the box, a quarter of them in the first layer of cells
// above the walls
inline std::vector<Vector3d> points(const Uint gdim, const std::size_t n, const double layer){
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> any(0.01, 0.99), low(0.001, layer);
  std::vector<Vector3d> pts;
  for (std::size_t i = 0; i < n; ++i){
    Vector3d x = Vector3d::Zero();
    for (Uint d = 0; d < gdim; ++d) x[d] = any(rng);
    if (i % 4 == 0) x[i/4 % gdim] = low(rng);
    pts.push_back(x);
  }
  return pts;
}

}  // namespace

#endif
