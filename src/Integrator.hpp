#ifndef __INTEGRATOR_HPP
#define __INTEGRATOR_HPP

#include "typedefs.hpp"
#include "Interpol.hpp"

class Integrator {
public:
  Integrator(Interpol* intp) { this->intp = intp; };
  virtual ~Integrator() { };
  virtual Vector3d integrate(const Vector3d& x, const double t, const double dt) = 0;
  Uint get_accepted() const { return n_accepted; };
  Uint get_declined() const { return n_declined; };
  bool is_stuck() const { return stuck; };
protected:
  Interpol* intp;
  Uint n_accepted = 0;
  Uint n_declined = 0;
  bool stuck;
};

#endif