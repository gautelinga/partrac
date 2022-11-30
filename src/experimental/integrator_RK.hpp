#ifndef __EXP_INTEGRATOR_RK_HPP
#define __EXP_INTEGRATOR_RK_HPP

#include "typedefs.hpp"
#include "integrator.hpp"
#include <math.h>

class Integrator_RK4 : public Integrator {
public:
  Integrator_RK4();
  ~Integrator_RK4() {};
  template<typename InterpolType, typename T>
  std::set<Uint> step(InterpolType&, T&, Real t, Real dt);
protected:
};

Integrator_RK4::Integrator_RK4() : Integrator() {
    std::cout << "Selecting Runge-Kutta 4 scheme" << std::endl;
}

template<typename InterpolType, typename T>
std::set<Uint> Integrator_RK4::step(InterpolType& intp, T& ps, const Real t, const Real dt) {
    std::set<Uint> outside_nodes;
    Uint i = 0;
    for (auto & particle : ps.particles() ){
        Vector x = particle.x();

        intp.probe(x, t);
        Vector k1 = intp.get_u();
        intp.probe(x + k1 * dt/2, t + dt/2);
        Vector k2 = intp.get_u();
        intp.probe(x + k2 * dt/2, t + dt/2);
        Vector k3 = intp.get_u();
        intp.probe(x + k3 * dt, t + dt);
        Vector k4 = intp.get_u();

        Vector dx = (k1 + 2*k2 + 2*k3 + k4) * dt/6;
        intp.probe(x + dx, t+dt);
        if (intp.inside_domain()){
            ++n_accepted;
            particle.x() = x + dx;
        }
        else {
            //++n_stuck;
            outside_nodes.insert(i);
            ++n_declined;
        }
        ++i;
    }
    return outside_nodes;
}

#endif