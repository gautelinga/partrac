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
  template<typename InterpolType, typename T>
  std::set<Uint> step_vec(InterpolType&, T&, Real t, Real dt);
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
        int cell_id = particle.cell_id(); // to accelerate search

        intp.probe(x, t, cell_id);
        Vector k1 = intp.get_u();
        intp.probe(x + k1 * dt/2, t + dt/2, cell_id);
        Vector k2 = intp.get_u();
        intp.probe(x + k2 * dt/2, t + dt/2, cell_id);
        Vector k3 = intp.get_u();
        intp.probe(x + k3 * dt, t + dt, cell_id);
        Vector k4 = intp.get_u();

        Vector dx = (k1 + 2*k2 + 2*k3 + k4) * dt/6;
        intp.probe(x + dx, t+dt, cell_id);
        if (intp.inside_domain()){
            ++n_accepted;
            particle.x() = x + dx;
            particle.cell_id() = cell_id;
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

template<typename InterpolType, typename T>
std::set<Uint> Integrator_RK4::step_vec(InterpolType& intp, T& ps, const Real t, const Real dt) {
    std::set<Uint> outside_nodes;
    Uint i = 0;
    for (auto & particle : ps.particles() ){
        Vector x = particle.x();
        Vector n = particle.n();
        int cell_id = particle.cell_id(); // to accelerate search

        intp.probe(x, t, cell_id);
        Vector k1 = intp.get_u();
        Matrix J1 = intp.get_J();
        Vector F1 = J1 * n;

        intp.probe(x + k1 * dt/2, t + dt/2, cell_id);
        Vector k2 = intp.get_u();
        Vector n2 = n + F1 * dt/2;

        Matrix J2 = intp.get_J();
        Vector F2 = J2 * n2;

        intp.probe(x + k2 * dt/2, t + dt/2, cell_id);
        Vector k3 = intp.get_u();
        Vector n3 = n + F2 * dt/2;
        
        Matrix J3 = intp.get_J();
        Vector F3 = J3 * n3;

        intp.probe(x + k3 * dt, t + dt, cell_id);
        Vector k4 = intp.get_u();
        Vector n4 = n + F3 * dt;

        Matrix J4 = intp.get_J();
        Vector F4 = J4 * n4;

        Vector dx = (k1 + 2*k2 + 2*k3 + k4) * dt/6;
        Vector el = n + (F1 + 2*F2 + 2*F3 + F4) * dt/6;

        intp.probe(x + dx, t+dt, cell_id);
        if (intp.inside_domain()){
            ++n_accepted;
            particle.x() = x + dx;
            particle.n() = el/el.norm();
            particle.w() += log(el.norm());

            Matrix J = intp.get_J();
            particle.S() = particle.n().transpose() * J * particle.n();
            particle.cell_id() = cell_id;
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