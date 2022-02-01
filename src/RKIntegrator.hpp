#ifndef __RKINTEGRATOR_HPP
#define __RKINTEGRATOR_HPP

#include "typedefs.hpp"
#include "Interpol.hpp"
#include "Integrator.hpp"
#include <math.h>

class RK4Integrator : public Integrator {
public:
  RK4Integrator(Interpol* intp);
  ~RK4Integrator() {};
  Vector3d integrate(const Vector3d& x, const double t, const double dt);
protected:
};

RK4Integrator::RK4Integrator(Interpol* intp) : Integrator(intp) {
    std::cout << "Selecting Runge-Kutta 4 scheme" << std::endl;
}

Vector3d RK4Integrator::integrate(const Vector3d& x, const double t, const double dt) {
    /*
    k1 = f(t_n, y_n)
    k2 = f(t_n + h/2, y_n + (h/2)*k1)
    k3 = f(t_n + h/2, y_n + (h/2)*k2)
    k4 = f(t_n + h,   y_n + h * k3)
    return y_n + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    */
    intp->probe(x, t);
    Vector3d k1 = intp->get_u();
    intp->probe(x + k1 * dt/2, t + dt/2);
    Vector3d k2 = intp->get_u();
    intp->probe(x + k2 * dt/2, t + dt/2);
    Vector3d k3 = intp->get_u();
    intp->probe(x + k3 * dt, t + dt);
    Vector3d k4 = intp->get_u();

    Vector3d dx = (k1 + 2*k2 + 2*k3 + k4) * dt/6;

    intp->probe(x+dx, t+dt);
    if (intp->inside_domain()){
        stuck = false;
        ++n_accepted;
        return dx;
    }
    else {
        stuck = true;
        ++n_declined;
        return {0., 0., 0.};
    }
}

#endif