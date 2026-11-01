#ifndef __INTEGRATOR_HPP
#define __INTEGRATOR_HPP

#include "typedefs.hpp"
#include "Interpol.hpp"
#include "ParticleSet.hpp"
#include "TransportElement.hpp"

// Accept/decline tally
class Integrator {
public:
  Integrator() {};
  virtual ~Integrator() { };
  Uint get_accepted() const { return n_accepted; };
  Uint get_declined() const { return n_declined; };
protected:
  // A carried element under a second-order step reads (U.grad) gradU
  static void prepare_carried(Interpol& intp, const TransportElement E, const int int_order){
    if (E != TransportElement::Point && int_order > 1)
      intp.set_needs_convective(true);
  }
  Uint n_accepted = 0;
  Uint n_declined = 0;
};

#endif