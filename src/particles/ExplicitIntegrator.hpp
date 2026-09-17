#ifndef __EXPLICITINTEGRATOR_HPP
#define __EXPLICITINTEGRATOR_HPP

#include <algorithm>
#include <vector>
#include "typedefs.hpp"
#include "Interpol.hpp"
#include "Integrator.hpp"
#include "TransportElement.hpp"
#include <math.h>

class ExplicitIntegrator : public Integrator {
public:
    ExplicitIntegrator(const double Dm, const int int_order, std::vector<std::mt19937>& gens);
    ~ExplicitIntegrator() {};
    // One loop for all transport elements
    template<TransportElement E = TransportElement::Point, typename Interp>
    std::vector<Uint> step(Interp& intp, ParticleSet& ps, const double t, const double dt);
protected:
  double Dm;
  int int_order;
  std::vector<std::mt19937>& gens;
  std::normal_distribution<double> rnd_normal;
};

inline ExplicitIntegrator::ExplicitIntegrator(const double Dm, const int int_order, std::vector<std::mt19937>& gens) : Dm(Dm), int_order(int_order), gens(gens), rnd_normal(0.0, 1.0) {
}

#endif
