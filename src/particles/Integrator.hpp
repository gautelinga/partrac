#ifndef __INTEGRATOR_HPP
#define __INTEGRATOR_HPP

#include "typedefs.hpp"
#include "Interpol.hpp"
#include "ParticleSet.hpp"

// Accept/decline tally
class Integrator {
public:
  Integrator() {};
  virtual ~Integrator() { };
  Uint get_accepted() const { return n_accepted; };
  Uint get_declined() const { return n_declined; };
protected:
  Uint n_accepted = 0;
  Uint n_declined = 0;
};

#endif