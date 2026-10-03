#ifndef __CELLSINTEGRATOR_HPP
#define __CELLSINTEGRATOR_HPP

// RK4 held in the cells of a mesh field, each step cut where the path meets a facet

#include <cstdint>
#include <iostream>
#include <vector>
#include "typedefs.hpp"
#include "RKIntegrator.hpp"
#include "TransportElement.hpp"

// What the steps in cells met
struct CellsCounts {
  std::uint64_t steps = 0;         // held steps, landed or rejected
  std::uint64_t crossings = 0;     // immediate ones included
  std::uint64_t relocations = 0;
  std::uint64_t immediate = 0;     // crossings without a step
  std::uint64_t rejects = 0;
  std::uint64_t fallbacks = 0;     // major steps taken by RK4 substeps
  CellsCounts& operator+=(const CellsCounts& o){
    steps += o.steps; crossings += o.crossings; relocations += o.relocations;
    immediate += o.immediate; rejects += o.rejects; fallbacks += o.fallbacks;
    return *this;
  }
};

// Scheme RK4cells; a field without cells (analytic) takes RK4's loop
class RK4CellsIntegrator : public RK4Integrator {
public:
  RK4CellsIntegrator() : RK4Integrator("RK4cells") {}
  template<TransportElement E = TransportElement::Point, typename Interp>
  std::vector<Uint> step(Interp& intp, ParticleSet& ps, const double t, const double dt);
  const CellsCounts& counts() const { return counts_; }
  // The counts per particle-step
  void report(std::ostream& os) const {
    if (particle_steps_ == 0)
      return;
    const double n = double(particle_steps_);
    os << "RK4cells: " << particle_steps_ << " particle-steps; per particle-step " << counts_.steps/n
       << " held steps, " << counts_.crossings/n << " crossings, " << counts_.relocations/n
       << " relocations, " << counts_.immediate/n << " immediate crossings, " << counts_.rejects/n
       << " rejects, " << counts_.fallbacks/n << " fallbacks" << std::endl;
  }
private:
  CellsCounts counts_;
  std::uint64_t particle_steps_ = 0;
};

#endif
