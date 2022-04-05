#ifndef __INTEGRATOR_HPP
#define __INTEGRATOR_HPP

#include "typedefs.hpp"
#include "Interpol.hpp"
#include "ParticleSet.hpp"

class Integrator {
public:
  //Integrator(std::shared_ptr<Interpol> intp) { this->intp = intp; };
  Integrator() {};
  virtual ~Integrator() { };
  //virtual Vector3d integrate(const Vector3d& x, const double t, const double dt) = 0;
  Uint get_accepted() const { return n_accepted; };
  Uint get_declined() const { return n_declined; };
  //bool is_stuck() const { return stuck; };
  //template<typename InterpolType, typename T>
  //virtual std::set<Uint> step(InterpolType& intp, T& ps, const double t, const double dt) = 0;
  virtual std::set<Uint> step(ParticleSet& ps, const double t, const double dt) = 0;
protected:
  //std::shared_ptr<Interpol> intp;
  Uint n_accepted = 0;
  Uint n_declined = 0;
  //bool stuck;
};

#endif