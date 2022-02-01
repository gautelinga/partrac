#ifndef __EXPLICITINTEGRATOR_HPP
#define __EXPLICITINTEGRATOR_HPP

#include "typedefs.hpp"
#include "Interpol.hpp"
#include "Integrator.hpp"
#include <math.h>

class ExplicitIntegrator : public Integrator {
public:
  ExplicitIntegrator(Interpol* intp, const double Dm, const int int_order, std::mt19937& gen);
  ~ExplicitIntegrator() {};
  Vector3d integrate(const Vector3d& x, const double t, const double dt);
protected:
  double Dm;
  int int_order;
  std::mt19937& gen;
  std::normal_distribution<double> rnd_normal;
};

ExplicitIntegrator::ExplicitIntegrator(Interpol* intp, const double Dm, const int int_order, std::mt19937& gen) : Integrator(intp), gen(gen), Dm(Dm), int_order(int_order), rnd_normal(0.0, 1.0) {
}

Vector3d ExplicitIntegrator::integrate(const Vector3d& x, const double t, const double dt) {
    intp->probe(x, t);
    Vector3d dx_rw = intp->get_u() * dt;

    // Set elongation
    //ps.e_rw[irw] = 0.;

    // Second-order terms
    if (int_order >= 2){
        //dx_rw += 0.5*a_rw[irw]*dt2;
        dx_rw += 0.5 * (intp->get_a() + intp->get_Ju()) * dt * dt;
    }
    if (Dm > 0.0){
        Vector3d eta = {rnd_normal(gen),
                        rnd_normal(gen),
                        rnd_normal(gen)};
        double sqrt2Dmdt = sqrt(2*Dm*dt);
        dx_rw += sqrt2Dmdt * eta;
    }
    intp->probe(x+dx_rw, t+dt);
    if (intp->inside_domain()){
        stuck = false;
        ++n_accepted;
        return dx_rw;
    }
    else {
        ++n_declined;
        stuck = true;
        return {0., 0., 0.};
    }
}

#endif