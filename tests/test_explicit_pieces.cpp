// A noisy explicit step whose end leaves the fluid, on a field that does not
// reflect, is taken again in pieces: halves, a failed piece halved again down
// to dt/64, each piece taking its share h/dt of the step's noise. The field is
// Taylor-Couette with the rolls and the decay off, a uniform azimuthal speed
// 1 between the cylinders, so a whole step from near the outer cylinder ends
// beyond it while pieces follow the circle. The noise the step draws is drawn
// again here from a copy of the generator, and every end, and every declined
// step, must be the rule's with that noise.
#include <catch2/catch.hpp>

#include <cmath>
#include <fstream>
#include <memory>
#include <random>
#include <vector>
#include <omp.h>

#include "AnalyticInterpol.hpp"
#include "ParticleSet.hpp"
#include "case_dir.hpp"
#include "stepping.hpp"

namespace {

// One thread for the step: the noise of particle i is then the generator's i-th triple
struct OneThread {
  const int before = omp_get_max_threads();
  OneThread(){ omp_set_num_threads(1); }
  ~OneThread(){ omp_set_num_threads(before); }
};

// The end of the step from x with this noise, by the rule of halves; false: declined
bool rule_end(AnalyticInterpol& intp, Vector3d x, const double t, const double dt, const Vector3d& noise,
              Vector3d& end, int& pieces){
  CellPos pos;
  PointValues pv(intp.get_U0());
  intp.evaluate(x, t, pos, pv);
  Vector3d dx = pv.get_u() * dt;
  dx += noise;
  pieces = 1;
  if (intp.locate(x + dx, t + dt, pos)){
    end = x + dx;
    return true;
  }
  constexpr int units = 64;
  const double unit = dt/units;
  int u = 0, size = units/2;
  pieces = 0;
  while (u < units){
    const double ta = t + u*unit, h = size*unit;
    intp.evaluate(x, ta, pos, pv);
    Vector3d d = pv.get_u() * h;
    d += noise * (h/dt);
    if (intp.locate(x + d, ta + h, pos)){
      x += d;
      u += size;
      ++pieces;
      while (size < units/2 && u % (2*size) == 0)
        size *= 2;
    }
    else if (size == 1)
      return false;
    else
      size /= 2;
  }
  end = x;
  return true;
}

}  // namespace

TEST_CASE("A noisy explicit step beyond a wall it cannot reflect off is taken in pieces sharing its noise", "[pieces]"){
  CaseDir c("explicit_pieces");
  {
    std::ofstream f(c.file("expr_params.dat"));
    f << "t_min=0.0\nt_max=1e7\nx_min=-2.5\ny_min=-2.5\nz_min=-1.5\nx_max=2.5\ny_max=2.5\nz_max=1.5\n"
         "Lx=5.0\nLy=5.0\nLz=3.0\nexpression=taylor_couette\np_inf=0.0\nrho=1.0\n"
         "R=2.5\nH=3.0\nK=0.0\na=0.4\nc1=0.0\nc2=1.2\nc3=0.0\nx0=0.0\ny0=0.0\nz0=0.0\n";
  }
  auto intp = std::make_shared<AnalyticInterpol>(c.file("expr_params.dat"));
  // Not enable_reflection: the field cannot reflect
  REQUIRE(!intp->can_reflect);
  const double t = 0., dt = 0.5, Dm = 1e-4, s0 = 2.48;
  intp->update(t);

  const Uint n = 96;
  std::vector<Vector3d> pts;
  for (Uint i = 0; i < n; ++i){
    const double a = 2*M_PI*(i + 0.25)/n;
    pts.push_back(Vector3d(s0*std::cos(a), s0*std::sin(a), 0.));
  }
  ParticleSet ps(intp, n);
  ps.carry(TransportElement::Point);
  ps.add(pts, 0);

  std::vector<std::mt19937> gens(omp_get_max_threads());
  for (std::size_t k = 0; k < gens.size(); ++k) gens[k].seed(7 + k);
  std::mt19937 g = gens[0];
  ExplicitIntegrator integrator(Dm, 1, gens);
  std::vector<Uint> out;
  {
    OneThread one;
    out = explicit_step<TransportElement::Point>(integrator, *intp, ps, t, dt);
  }

  // The step's draws, in particle order
  std::normal_distribution<double> nd(0., 1.);
  const double sigma = std::sqrt(2*Dm*dt);
  std::vector<Uint> declined;
  Uint rescued = 0;
  for (Uint i = 0; i < n; ++i){
    const Vector3d eta = {nd(g), nd(g), nd(g)};
    Vector3d end;
    int pieces = 0;
    INFO("particle " << i);
    if (!rule_end(*intp, pts[i], t, dt, sigma*eta, end, pieces)){
      declined.push_back(i);
      REQUIRE(ps.x(i) == pts[i]);
      continue;
    }
    rescued += pieces > 1;
    REQUIRE((ps.x(i) - end).norm() <= 1e-14);
    // dt itself, not the sum of the pieces
    REQUIRE(ps.t_loc(i) == dt);
  }
  INFO(rescued << " rescued, " << declined.size() << " declined");
  REQUIRE(out == declined);
  REQUIRE(rescued > n/2);
}
