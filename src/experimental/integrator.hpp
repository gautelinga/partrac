#ifndef __EXP_INTEGRATOR_HPP
#define __EXP_INTEGRATOR_HPP

#include "typedefs.hpp"

class Integrator {
public:
  Integrator() { };
  ~Integrator() { };
  Uint get_accepted() const { return n_accepted; };
  Uint get_declined() const { return n_declined; };
  template<typename InterpolType, typename T>
  std::set<Uint> step(InterpolType&, T&, Real t, Real dt);
protected:
  Uint n_accepted = 0;
  Uint n_declined = 0;
};

#endif