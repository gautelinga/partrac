#ifndef __EXP_INTEGRATOR_EXPLICIT_H
#define __EXP_INTEGRATOR_EXPLICIT_H

#include "typedefs.hpp"
#include "integrator.hpp"
#include <math.h>

class Integrator_Explicit : public Integrator {
public:
  Integrator_Explicit(const Real Dm, const int int_order, std::mt19937& gen);
  ~Integrator_Explicit() {};
  template<typename InterpolType, typename T>
  std::set<Uint> step(InterpolType&, T&, Real t, Real dt);
protected:
  Real Dm;
  int int_order;
  std::mt19937& gen;
  std::normal_distribution<Real> rnd_normal;
};

Integrator_Explicit::Integrator_Explicit(const Real Dm, const int int_order, std::mt19937& gen)
  : Integrator(), gen(gen), Dm(Dm), int_order(int_order), rnd_normal(0.0, 1.0) {
    std::cout << "Choosing an explicit integrator of order " << int_order << " with diffusivity " << Dm << std::endl;
}

template<typename InterpolType, typename T>
std::set<Uint> Integrator_Explicit::step(InterpolType& intp, T& ps, const Real t, const Real dt) {
    std::set<Uint> outside_nodes;
    Real sqrt2Dmdt = sqrt(2*Dm*dt);
    Uint i = 0;
    for (auto & particle : ps.particles() ){
        Vector x = particle.x();
        intp.probe(x, t);
        Vector dx = intp.get_u() * dt;

        // Second-order terms
        if (int_order >= 2){
            //dx_rw += 0.5*a_rw[irw]*dt2;
            dx += 0.5 * (intp.get_a() + intp.get_Ju()) * dt * dt;
        }
        if (Dm > 0.0){
            Vector eta = {rnd_normal(gen),
                            rnd_normal(gen),
                            rnd_normal(gen)};
            dx += sqrt2Dmdt * eta;
        }
        intp.probe(x+dx, t+dt);
        if (intp.inside_domain()){
            ++n_accepted;
            particle.x() = x + dx;
        }
        else {
            outside_nodes.insert(i);
            ++n_declined;
        }
        ++i;
    }
    return outside_nodes;
}

#endif