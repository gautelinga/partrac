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
    RK4Integrator();
    ~RK4Integrator() {};
  // One loop for all transport elements
  template<TransportElement E = TransportElement::Point, typename Interp>
  std::vector<Uint> step(Interp& intp, ParticleSet& ps, const double t, const double dt);
protected:
};

inline RK4Integrator::RK4Integrator() : Integrator() {
    std::cout << "Selecting Runge-Kutta 4 scheme" << std::endl;
}

#endif
