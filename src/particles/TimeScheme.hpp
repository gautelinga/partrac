#ifndef __TIMESCHEME_HPP
#define __TIMESCHEME_HPP

#include <optional>
#include <random>
#include <set>
#include <vector>
#include "typedefs.hpp"
#include "Params.hpp"
#include "Interpol.hpp"
#include "ParticleSet.hpp"
#include "Integrator.hpp"
#include "ExplicitIntegrator.hpp"
#include "RKIntegrator.hpp"
#include "TransportElement.hpp"
#include "stepping.hpp"

// Time scheme: explicit or RK4, dispatched on the concrete interpolator
class TimeScheme {
public:
  TimeScheme(const partrac::Params& prm, std::vector<std::mt19937>& gens){
    // Schema allows only these
    if (prm.get<std::string>("scheme") == "RK4")
      rk4.emplace();
    else
      explicit_.emplace(prm.get<double>("Dm"), prm.get<int>("int_order"), gens);
  }
  // Accept/decline tally
  Integrator& counters(){
    return rk4 ? static_cast<Integrator&>(*rk4) : static_cast<Integrator&>(*explicit_);
  }
  template<TransportElement E = TransportElement::Point>
  std::vector<Uint> step(Interpol& intp, ParticleSet& ps, const double t, const double dt){
    return rk4 ? rk4_step<E>(*rk4, intp, ps, t, dt)
               : explicit_step<E>(*explicit_, intp, ps, t, dt);
  }
private:
  std::optional<ExplicitIntegrator> explicit_;
  std::optional<RK4Integrator> rk4;
};

#endif
