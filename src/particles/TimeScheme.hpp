#ifndef __TIMESCHEME_HPP
#define __TIMESCHEME_HPP

#include <optional>
#include <ostream>
#include <string>
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
#include "CellsIntegrator.hpp"
#include "TransportElement.hpp"
#include "stepping.hpp"

// Time scheme: explicit, RK4 or RK4cells, dispatched on the concrete interpolator
class TimeScheme {
public:
  TimeScheme(const partrac::Params& prm, std::vector<std::mt19937>& gens){
    // Schema allows only these
    const std::string scheme = prm.get<std::string>("scheme");
    if (scheme == "RK4")
      rk4.emplace();
    else if (scheme == "RK4cells")
      rk4cells.emplace();
    else
      explicit_.emplace(prm.get<double>("Dm"), prm.get<int>("int_order"), gens);
  }
  // RK4cells in cells: the gradient for its prediction; a field it cannot take refused
  void prepare(Interpol& intp) const {
    if (rk4cells && steps_in_cells(intp))
      intp.set_needs_gradient(true);
  }
  // Accept/decline tally
  Integrator& counters(){
    if (rk4cells) return *rk4cells;
    return rk4 ? static_cast<Integrator&>(*rk4) : static_cast<Integrator&>(*explicit_);
  }
  template<TransportElement E = TransportElement::Point>
  std::vector<Uint> step(Interpol& intp, ParticleSet& ps, const double t, const double dt){
    if (rk4cells) return cells_step<E>(*rk4cells, intp, ps, t, dt);
    return rk4 ? rk4_step<E>(*rk4, intp, ps, t, dt)
               : explicit_step<E>(*explicit_, intp, ps, t, dt);
  }
  // RK4cells' counts, at the end of a run
  void report(std::ostream& os) const {
    if (rk4cells) rk4cells->report(os);
  }
private:
  std::optional<ExplicitIntegrator> explicit_;
  std::optional<RK4Integrator> rk4;
  std::optional<RK4CellsIntegrator> rk4cells;
};

#endif
