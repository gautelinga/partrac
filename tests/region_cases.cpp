// The writers of region_cases.hpp, in a unit of their own: dolfin's headers
// and the writers are compiled once for the tests that share them.
#ifdef USE_DOLFIN
#include "region_cases.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>

#include <dolfin.h>
#include "H5Cpp.h"

#include "Tet.hpp"
#include "Triangle.hpp"
#include "divfree_poly.hpp"
#include "taylor_hood.hpp"

template<typename Cell>
std::shared_ptr<dolfin::Mesh> unit_mesh(const std::size_t n){
  std::shared_ptr<dolfin::Mesh> mesh;
  if constexpr (dim_of<Cell> == 2) mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitSquareMesh(n, n));
  else                             mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitCubeMesh(n, n, n));
  mesh->init();
  return mesh;
}

template<int D>
Vector3d periodic_u(const Vector3d& x, const int k){
  const double s = 2.*M_PI;
  Vector3d u(std::sin(s*(x[1] + 0.1*k)) + 0.3*std::cos(s*x[0]),
             std::sin(s*(x[0] + 0.2*k)) + 0.2*std::sin(s*(x[0] + x[1])), 0.);
  if constexpr (D == 3){
    u[0] += 0.2*std::cos(s*x[2]);
    u[2] = std::sin(s*(x[0] + 0.3*k)) + 0.4*std::cos(s*(x[1] - x[2]));
  }
  return u;
}

// At rest on the walls of the unit box, as the near-wall rule wants
template<int D>
Vector3d rest_u(const Vector3d& x, const int k){
  double b = 1.;
  for (int d = 0; d < D; ++d) b *= x[d]*(1. - x[d]);
  Vector3d u = Vector3d::Zero();
  for (int d = 0; d < D; ++d)
    u[d] = b*(1. + double(k))*(d == 0 ? 4.*(1. + x[1]) : 1. + x[0] + 0.5*double(d));
  return u;
}

template<int D>
class VExpr : public dolfin::Expression {
public:
  VExpr(const int k, const bool rest, const Field& field = {})
    : dolfin::Expression(D), k_(k), rest_(rest), field_(field) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    const Vector3d y(x[0], x[1], D == 3 ? x[2] : 0.);
    const Vector3d u = field_ ? field_(y, k_) : rest_ ? rest_u<D>(y, k_) : periodic_u<D>(y, k_);
    for (int d = 0; d < D; ++d) v[d] = u[d];
  }
private:
  int k_;
  bool rest_;
  Field field_;
};

class SExpr : public dolfin::Expression {
public:
  explicit SExpr(const int k) : k_(k) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    v[0] = std::cos(2.*M_PI*x[0]) + 0.1*double(k_);
  }
private:
  int k_;
};

std::shared_ptr<const dolfin::SubDomain> pbc_of(const std::vector<bool>& per, const int D){
  if (!per[0] && !per[1] && !per[2]) return nullptr;
  return std::make_shared<PeriodicBC>(per, Vector3d::Zero(), Vector3d(1., 1., D == 3 ? 1. : 0.), Uint(D));
}

std::string periodic_keys(const std::vector<bool>& per){
  return std::string("periodic_x=") + (per[0] ? "true" : "false") + "\nperiodic_y=" + (per[1] ? "true" : "false")
       + "\nperiodic_z=" + (per[2] ? "true" : "false") + "\n";
}

template<typename Cell>
void write_stamped(const CaseDir& c, const std::size_t n, const std::string& u_el,
                   const std::vector<bool>& per, const bool rest, const Warp& warp,
                   const Field& field){
  constexpr int D = dim_of<Cell>;
  auto mesh = unit_mesh<Cell>(n);
  if (warp){
    std::vector<double>& xs = mesh->coordinates();
    for (std::size_t i = 0; i < xs.size(); i += D){
      const Vector3d y = warp(Vector3d(xs[i], xs[i + 1], D == 3 ? xs[i + 2] : 0.));
      for (int d = 0; d < D; ++d) xs[i + std::size_t(d)] = y[d];
    }
  }
  std::shared_ptr<dolfin::FunctionSpace> V, P;
  Uint ncu = 0, ncp = 0;
  taylor_hood_spaces<Cell>(u_el, "P1", true, mesh, pbc_of(per, D), V, P, ncu, ncp);
  const bool xdmf = u_el == "P1";
  {
    dolfin::HDF5File f(MPI_COMM_WORLD, (c.path / "mesh.h5").string(), "w");
    f.write(*mesh, "mesh");
  }
  std::unique_ptr<dolfin::XDMFFile> xf[2];
  if (xdmf)
    for (int i = 0; i < 2; ++i){
      xf[i] = std::make_unique<dolfin::XDMFFile>(mesh->mpi_comm(), (c.path / (i ? "p.xdmf" : "u.xdmf")).string());
      xf[i]->parameters["functions_share_mesh"] = true;
      xf[i]->parameters["rewrite_function_mesh"] = false;
    }
  for (int k = 0; k < 2; ++k){
    dolfin::Function u(V), p(P);
    u.interpolate(VExpr<D>(k, rest, field));
    p.interpolate(SExpr(k));
    dolfin::HDF5File f(MPI_COMM_WORLD, (c.path / ("up_" + std::to_string(k) + ".h5")).string(), "w");
    f.write(u, "u");
    f.write(p, "p");
    if (xdmf){
      xf[0]->write(u, double(k));
      xf[1]->write(p, double(k));
    }
  }
  for (auto& f : xf) if (f) f->close();
  std::ofstream(c.path / "timestamps.dat") << "0\tup_0.h5\n1\tup_1.h5\n";
  std::ofstream(c.file("h5_params.dat")) << "velocity_space=" << u_el
    << "\npressure_space=P1\ntimestamps=timestamps.dat\nmesh=mesh.h5\n" << periodic_keys(per);
  if (xdmf)
    std::ofstream(c.file("xdmf_params.dat")) << "u=u.xdmf\np=p.xdmf\n" << periodic_keys(per);
}

template<typename Cell>
void write_freq(const CaseDir& c, const std::size_t n){
  constexpr int D = dim_of<Cell>;
  auto mesh = unit_mesh<Cell>(n);
  const std::vector<bool> per = {true, true, D == 3};
  {
    dolfin::HDF5File f(MPI_COMM_WORLD, (c.path / "mesh.h5").string(), "w");
    f.write(*mesh, "mesh");
  }
  std::ofstream stamps(c.path / "freqstamps.dat");
  for (int k = 0; k < 2; ++k){
    std::shared_ptr<dolfin::FunctionSpace> V, P;
    Uint ncu = 0, ncp = 0;
    taylor_hood_spaces<Cell>("P2", "P1", true, mesh, pbc_of(per, D), V, P, ncu, ncp);
    dolfin::Function u(V), p(P);
    u.interpolate(VExpr<D>(k, false));
    p.interpolate(SExpr(k));
    const std::string name = "up_" + std::to_string(k) + ".h5";
    dolfin::HDF5File f(MPI_COMM_WORLD, (c.path / name).string(), "w");
    f.write(u, "u");
    f.write(p, "p");
    stamps << 0.125*k << " " << 1. - 0.5*k << " " << name << "\n";
  }
  stamps.close();
  std::ofstream(c.params()) << "velocity_space=P2\npressure_space=P1\nfreqstamps=freqstamps.dat\n"
    << "mesh=mesh.h5\n" << periodic_keys(per) << "tau=1\nt_min=0\nt_max=1e8\n";
}

// Quadratic and divergence-free, so every cell of its P2 interpolant balances;
// independent of x where the mesh is periodic along x
template<int D>
Vector3d split_u(const Vector3d& x, const int k, const bool periodic){
  Vector3d u = poly_u(x, D);
  if (periodic){
    const double y = x[1], z = x[2];
    u = D == 2 ? Vector3d(0.8*y*y - 0.4*y + 0.35, 0.6, 0.)
               : Vector3d(0.3*y*y - 0.5*z*z + 0.7*y*z + 0.2*y + 0.15, 0.9*y*z + 0.25*y, -0.45*z*z - 0.25*z + 0.4);
  }
  return (1. + double(k))*u;
}

template<int D>
class SplitExpr : public dolfin::Expression {
public:
  SplitExpr(const int k, const bool periodic) : dolfin::Expression(D), k_(k), periodic_(periodic) {}
  void eval(dolfin::Array<double>& v, const dolfin::Array<double>& x) const override {
    const Vector3d u = split_u<D>(Vector3d(x[0], x[1], D == 3 ? x[2] : 0.), k_, periodic_);
    for (int d = 0; d < D; ++d) v[d] = u[d];
  }
private:
  int k_;
  bool periodic_;
};

// Every edge midpoint moved along its own edge, by an amount periodic in the
// box: no facet's flux changes, so every cell still balances, and the split's
// interior is no longer the polynomial's
template<int D>
void nudge_edges(dolfin::Function& u, const dolfin::FunctionSpace& V, const double amp){
  const std::vector<double> dc = V.tabulate_dof_coordinates();
  std::vector<double> base, vals;
  u.vector()->get_local(base);
  vals = base;
  for (dolfin::CellIterator c(*V.mesh()); !c.end(); ++c){
    std::vector<Vector3d> vx;
    for (dolfin::VertexIterator v(*c); !v.end(); ++v)
      vx.push_back(Vector3d(v->x(0), v->x(1), D == 3 ? v->x(2) : 0.));
    const auto dofs = V.dofmap()->cell_dofs(c->index());
    const int nn = int(dofs.size())/D;
    for (int k = 0; k < nn; ++k){
      Vector3d x = Vector3d::Zero();
      for (int d = 0; d < D; ++d) x[d] = dc[std::size_t(dofs[k])*D + std::size_t(d)];
      for (std::size_t a = 0; a < vx.size(); ++a)
        for (std::size_t b = a + 1; b < vx.size(); ++b){
          if ((0.5*(vx[a] + vx[b]) - x).norm() > 1e-12) continue;
          const bool a_first = std::lexicographical_compare(vx[a].data(), vx[a].data() + 3,
                                                            vx[b].data(), vx[b].data() + 3);
          const Vector3d tang = a_first ? Vector3d(vx[b] - vx[a]) : Vector3d(vx[a] - vx[b]);
          const double w = amp*std::sin(2.*M_PI*x[0] + 7.*x[1] + 3.*x[2]);
          for (int d = 0; d < D; ++d){
            const std::size_t i = std::size_t(dofs[k + d*nn]);
            vals[i] = base[i] + w*tang[d];
          }
        }
    }
  }
  u.vector()->set_local(vals);
  u.vector()->apply("insert");
}

template<typename Cell>
void write_split(const CaseDir& c, const std::size_t n, const bool periodic){
  constexpr int D = dim_of<Cell>;
  auto mesh = unit_mesh<Cell>(n);
  const std::vector<bool> per = {periodic, false, false};
  std::shared_ptr<dolfin::FunctionSpace> V, P;
  Uint ncu = 0, ncp = 0;
  taylor_hood_spaces<Cell>("P2", "P1", true, mesh, pbc_of(per, D), V, P, ncu, ncp);
  {
    dolfin::HDF5File f(MPI_COMM_WORLD, (c.path / "mesh.h5").string(), "w");
    f.write(*mesh, "mesh");
  }
  for (int k = 0; k < 2; ++k){
    dolfin::Function u(V), p(P);
    u.interpolate(SplitExpr<D>(k, periodic));
    nudge_edges<D>(u, *V, 0.3*(1. + 0.5*k));
    p.interpolate(SExpr(k));
    dolfin::HDF5File f(MPI_COMM_WORLD, (c.path / ("up_" + std::to_string(k) + ".h5")).string(), "w");
    f.write(u, "u");
    f.write(p, "p");
  }
  std::ofstream(c.path / "timestamps.dat") << "0\tup_0.h5\n1\tup_1.h5\n";
  std::ofstream(c.file("h5_params.dat")) << "velocity_space=P2\npressure_space=P1\ntimestamps=timestamps.dat\n"
    << "mesh=mesh.h5\ndivfree=true\n" << periodic_keys(per);
}

inline Vector3d lattice_u(const int n, const Vector3d& x, const int k){
  const double s = 2.*M_PI/double(n);
  return (1. + 0.5*k)*Vector3d(std::sin(s*(x[1] + 0.3*k)) + 0.3*std::cos(s*x[2]),
                               std::sin(s*(x[2] + 0.2)) + 0.4*std::cos(s*(x[0] - x[1])),
                               std::sin(s*x[0]) + 0.2*std::cos(s*(x[1] + x[2] + 0.1*k)));
}

void write_lattice(const CaseDir& c, const int n){
  const hsize_t dims[3] = {hsize_t(n), hsize_t(n), hsize_t(n)};
  const H5::DataSpace space(3, dims);
  // x fastest, as the solver writes them
  const auto at = [n](const int i, const int j, const int k){ return std::size_t((k*n + j)*n + i); };
  std::vector<int> solid(std::size_t(n)*n*n);
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) for (int k = 0; k < n; ++k)
    solid[at(i, j, k)] = lattice_solid(n, i, j, k);
  {
    H5::H5File f(c.file("output_is_solid.h5"), H5F_ACC_TRUNC);
    f.createDataSet("is_solid", H5::PredType::NATIVE_INT, space).write(solid.data(), H5::PredType::NATIVE_INT);
  }
  for (int st = 0; st < 2; ++st){
    std::vector<double> u[3], rho(solid.size()), p(solid.size());
    for (auto& v : u) v.resize(solid.size());
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) for (int k = 0; k < n; ++k){
      const std::size_t q = at(i, j, k);
      const Vector3d v = solid[q] ? Vector3d::Zero() : lattice_u(n, Vector3d(i, j, k), st);
      for (int d = 0; d < 3; ++d) u[d][q] = v[d];
      rho[q] = 1. + 0.01*std::cos(0.3*i + 0.2*j);
      p[q] = 0.1*std::sin(0.2*k + 0.1*st);
    }
    H5::H5File f(c.file("output_" + std::to_string(st) + ".h5"), H5F_ACC_TRUNC);
    const char* names[5] = {"u_x", "u_y", "u_z", "density", "pressure"};
    const std::vector<double>* data[5] = {&u[0], &u[1], &u[2], &rho, &p};
    for (int m = 0; m < 5; ++m)
      f.createDataSet(names[m], H5::PredType::NATIVE_DOUBLE, space).write(data[m]->data(), H5::PredType::NATIVE_DOUBLE);
  }
  std::ofstream(c.file("timestamps.dat")) << "0\toutput_0.h5\n1\toutput_1.h5\n";
  std::ofstream(c.file("felbm_params.dat")) << "timestamps=timestamps.dat\nis_solid_file=output_is_solid.h5\n";
}


template std::shared_ptr<dolfin::Mesh> unit_mesh<Triangle>(std::size_t);
template std::shared_ptr<dolfin::Mesh> unit_mesh<Tet>(std::size_t);
template void write_stamped<Triangle>(const CaseDir&, std::size_t, const std::string&, const std::vector<bool>&,
                                      bool, const Warp&, const Field&);
template void write_stamped<Tet>(const CaseDir&, std::size_t, const std::string&, const std::vector<bool>&,
                                 bool, const Warp&, const Field&);
template void write_freq<Triangle>(const CaseDir&, std::size_t);
template void write_freq<Tet>(const CaseDir&, std::size_t);
template void write_split<Triangle>(const CaseDir&, std::size_t, bool);
template void write_split<Tet>(const CaseDir&, std::size_t, bool);

#endif
