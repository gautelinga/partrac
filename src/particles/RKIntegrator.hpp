#ifndef __RKINTEGRATOR_HPP
#define __RKINTEGRATOR_HPP

#include <algorithm>
#include <vector>
#include "typedefs.hpp"
#include "Interpol.hpp"
#include "Integrator.hpp"
#include "TransportElement.hpp"
#include <omp.h>
#include <math.h>

class RK4Integrator : public Integrator {
public:
    RK4Integrator() : RK4Integrator("Runge-Kutta 4") {};
    ~RK4Integrator() {};
  // One loop for all transport elements
  template<TransportElement E = TransportElement::Point, typename Interp>
  std::vector<Uint> step(Interp& intp, ParticleSet& ps, const double t, const double dt);
protected:
  explicit RK4Integrator(const char* name) : Integrator() {
    std::cout << "Selecting " << name << " scheme" << std::endl;
  }
};

#endif
