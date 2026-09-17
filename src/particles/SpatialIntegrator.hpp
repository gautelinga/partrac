#ifndef __SPATIALINTEGRATOR_HPP
#define __SPATIALINTEGRATOR_HPP

#include <algorithm>
#include <vector>
#include "typedefs.hpp"
#include "Interpol.hpp"
#include "Integrator.hpp"
#include "ParticleSet.hpp"
#include "TransportElement.hpp"
#include <omp.h>
#include <math.h>

// Constant path-length steps along streamlines, fields frozen
class SpatialIntegrator : public Integrator {
public:
  SpatialIntegrator(const int int_order, const double u_min, const double dl_max, const double T);
  ~SpatialIntegrator() {};
  template<TransportElement E = TransportElement::Point, typename InterpolType>
  std::vector<Uint> step(InterpolType&, ParticleSet&, double t, double ds);
protected:
  double   m_u_min;
  double   m_dl_max;
  int      m_int_order;
  double   m_T;
};

inline SpatialIntegrator::SpatialIntegrator(const int int_order, const double u_min, const double dl_max, const double T)
  : Integrator(), m_u_min(u_min), m_dl_max(dl_max), m_int_order(int_order), m_T(T) {
    std::cout << "Choosing a spatial integrator of order " << int_order << "." << std::endl;
}

#endif
