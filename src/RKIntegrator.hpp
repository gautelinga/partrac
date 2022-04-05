#ifndef __RKINTEGRATOR_HPP
#define __RKINTEGRATOR_HPP

#include "typedefs.hpp"
#include "Interpol.hpp"
#include "Integrator.hpp"
#include <math.h>

class RK4Integrator : public Integrator {
public:
    //RK4Integrator(std::shared_ptr<Interpol> intp);
    RK4Integrator();
    ~RK4Integrator() {};
    //Vector3d integrate(const Vector3d& x, const double t, const double dt);
    //template<typename InterpolType, typename T>
    //std::set<Uint> step(InterpolType& intp, T& ps, const double t, const double dt);
    std::set<Uint> step(ParticleSet& ps, const double t, const double dt);
protected:
};

//RK4Integrator::RK4Integrator(std::shared_ptr<Interpol> intp) : Integrator(intp) {
RK4Integrator::RK4Integrator() : Integrator() {
    std::cout << "Selecting Runge-Kutta 4 scheme" << std::endl;
}

/*
Vector3d RK4Integrator::integrate(const Vector3d& x, const double t, const double dt) {
    // k1 = f(t_n, y_n)
    // k2 = f(t_n + h/2, y_n + (h/2)*k1)
    // k3 = f(t_n + h/2, y_n + (h/2)*k2)
    // k4 = f(t_n + h,   y_n + h * k3)
    // return y_n + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
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
}*/
//template<typename InterpolType, typename T>
//std::set<Uint> RK4Integrator::step(InterpolType& intp, T& ps, const double t, const double dt) {
std::set<Uint> RK4Integrator::step(ParticleSet& ps, const double t, const double dt) {
    std::set<Uint> outside_nodes;
    Vector3d dx, k1, k2, k3, k4;

    auto & intp = *ps.interpolator();

    for (Uint i=0; i < ps.N(); ++i){
        Vector3d x = ps.x(i);

        intp.probe(x, t);
        k1 = intp.get_u();
        intp.probe(x + k1 * dt/2, t + dt/2);
        k2 = intp.get_u();
        intp.probe(x + k2 * dt/2, t + dt/2);
        k3 = intp.get_u();
        intp.probe(x + k3 * dt, t + dt);
        k4 = intp.get_u();

        dx = (k1 + 2*k2 + 2*k3 + k4) * dt/6;

        intp.probe(x+dx, t+dt);
        if (intp.inside_domain()){
            ps.set_x(i, x + dx);
            ps.set_tau(i, ps.tau(i) + dt);
            ++n_accepted;
        }
        else {
            outside_nodes.insert(i);
            ++n_declined;
        }
    }
    return outside_nodes;
}

#endif