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
  template<typename InterpolType, typename T>
  std::set<Uint> step_tensor(InterpolType&, T&, Real t, Real dt);
protected:
};

Integrator_RK4::Integrator_RK4() : Integrator() {
    std::cout << "Selecting Runge-Kutta 4 scheme" << std::endl;
}

template<typename InterpolType, typename T>
std::set<Uint> Integrator_RK4::step(InterpolType& intp, T& ps, const Real t, const Real dt) {
    std::set<Uint> outside_nodes;
    // Particles are independent; only the tally is shared
    #pragma omp parallel
    {
    std::set<Uint> outside_nodes_loc;
    Uint n_accepted_loc = 0;
    Uint n_declined_loc = 0;
    #pragma omp for
    for (Uint i=0; i < ps.particles().size(); ++i){
        auto & particle = ps.particles()[i];
        Vector x = particle.x();
        int cell_id = particle.cell_id(); // to accelerate search
        PointValues ptvals(intp.get_U0());

        intp.locate(x, t, cell_id);
        intp.evaluate(x, t, cell_id, ptvals);
        Vector k1 = ptvals.get_u();
        intp.locate(x + k1 * dt/2, t + dt/2, cell_id);
        intp.evaluate(x + k1 * dt/2, t + dt/2, cell_id, ptvals);
        Vector k2 = ptvals.get_u();
        intp.locate(x + k2 * dt/2, t + dt/2, cell_id);
        intp.evaluate(x + k2 * dt/2, t + dt/2, cell_id, ptvals);
        Vector k3 = ptvals.get_u();
        intp.locate(x + k3 * dt, t + dt, cell_id);
        intp.evaluate(x + k3 * dt, t + dt, cell_id, ptvals);
        Vector k4 = ptvals.get_u();

        Vector dx = (k1 + 2*k2 + 2*k3 + k4) * dt/6;
        if (intp.locate(x + dx, t+dt, cell_id)){
            ++n_accepted_loc;
            particle.x() = x + dx;
            particle.cell_id() = cell_id;
        }
        else {
            //++n_stuck;
            outside_nodes_loc.insert(i);
            ++n_declined_loc;
        }
    }
    #pragma omp critical
    {
        outside_nodes.insert(outside_nodes_loc.begin(), outside_nodes_loc.end());
        n_accepted += n_accepted_loc;
        n_declined += n_declined_loc;
    }
    }
    return outside_nodes;
}

template<typename InterpolType, typename T>
std::set<Uint> Integrator_RK4::step_vec(InterpolType& intp, T& ps, const Real t, const Real dt) {
    std::set<Uint> outside_nodes;
    // Uint i = 0;
    #pragma omp parallel
    {
        std::set<Uint> outside_nodes_loc;
        Uint n_accepted_loc = 0;
        Uint n_declined_loc = 0;

        //for (auto & particle : ps.particles() ){
        #pragma omp for
        for (Uint i=0; i < ps.particles().size(); ++i)
        {
            auto & particle = ps.particles()[i];
            Vector x = particle.x();
            Vector n = particle.n();
            int cell_id = particle.cell_id(); // to accelerate search

            PointValues ptvals(intp.get_U0());
            
            Vector3d k1 = {0., 0., 0.};
            Vector3d k2 = {0., 0., 0.};
            Vector3d k3 = {0., 0., 0.};
            Vector3d k4 = {0., 0., 0.};
            
            Vector3d F1 = {0., 0., 0.};
            Vector3d F2 = {0., 0., 0.};
            Vector3d F3 = {0., 0., 0.};
            Vector3d F4 = {0., 0., 0.};

            bool is_inside = intp.locate(x, t, cell_id);
            if (is_inside)
            {
                intp.evaluate(x, t, cell_id, ptvals);
                k1 = ptvals.get_u();
                Matrix J1 = ptvals.get_J();
                F1 = J1 * n;
            }
            is_inside = intp.locate(x + k1 * dt/2, t + dt/2, cell_id);
            if (is_inside)
            {
                intp.evaluate(x + k1 * dt/2, t + dt/2, cell_id, ptvals);
                k2 = ptvals.get_u();
                Vector n2 = n + F1 * dt/2;
                Matrix J2 = ptvals.get_J();
                F2 = J2 * n2;
            }
            is_inside = intp.locate(x + k2 * dt/2, t + dt/2, cell_id);
            if (is_inside)
            {
                intp.evaluate(x + k2 * dt/2, t + dt/2, cell_id, ptvals);
                k3 = ptvals.get_u();
                Vector n3 = n + F2 * dt/2;
                Matrix J3 = ptvals.get_J();
                F3 = J3 * n3;
            }
            is_inside = intp.locate(x + k3 * dt, t + dt, cell_id);
            if (is_inside)
            {
                intp.evaluate(x + k3 * dt, t + dt, cell_id, ptvals);
                k4 = ptvals.get_u();
                Vector n4 = n + F3 * dt;
                Matrix J4 = ptvals.get_J();
                F4 = J4 * n4;
            }
            Vector dx = (k1 + 2*k2 + 2*k3 + k4) * dt/6;
            Vector el = n + (F1 + 2*F2 + 2*F3 + F4) * dt/6;

            is_inside = intp.locate(x + dx, t+dt, cell_id);
            if (is_inside){
                ++n_accepted_loc;
                particle.x() = x + dx;
                particle.n() = el/el.norm();
                particle.w() += log(el.norm());

                intp.evaluate(x + dx, t+dt, cell_id, ptvals);
                Matrix J = ptvals.get_J();
                particle.S() = particle.n().transpose() * J * particle.n();
                particle.cell_id() = cell_id;
            }
            else {
                //++n_stuck;
                outside_nodes_loc.insert(i);
                ++n_declined_loc;
            }
            //++i;
        }
        #pragma omp critical
        {
            outside_nodes.insert(outside_nodes_loc.begin(), outside_nodes_loc.end());
            n_accepted += n_accepted_loc;
            n_declined += n_declined_loc;
        }
    }
    return outside_nodes;
}

template<typename InterpolType, typename T>
std::set<Uint> Integrator_RK4::step_tensor(InterpolType& intp, T& ps, const Real t, const Real dt) {
    std::set<Uint> outside_nodes;
    #pragma omp parallel
    {
    std::set<Uint> outside_nodes_loc;
    Uint n_accepted_loc = 0;
    Uint n_declined_loc = 0;
    #pragma omp for
    for (Uint i=0; i < ps.particles().size(); ++i){
        auto & particle = ps.particles()[i];
        Vector x = particle.x();
        Matrix F = particle.F();
        int cell_id = particle.cell_id(); // to accelerate search
        PointValues ptvals(intp.get_U0());

        intp.locate(x, t, cell_id);
        intp.evaluate(x, t, cell_id, ptvals);
        Vector k1 = ptvals.get_u();
        
        Matrix J1 = ptvals.get_J();
        Matrix dFdt1 = J1 * F;

        intp.locate(x + k1 * dt/2, t + dt/2, cell_id);
        intp.evaluate(x + k1 * dt/2, t + dt/2, cell_id, ptvals);
        Vector k2 = ptvals.get_u();
        Matrix F2 = F + dFdt1 * dt/2;
        Matrix J2 = ptvals.get_J();
        Matrix dFdt2 = J2 * F2;

        intp.locate(x + k2 * dt/2, t + dt/2, cell_id);
        intp.evaluate(x + k2 * dt/2, t + dt/2, cell_id, ptvals);
        Vector k3 = ptvals.get_u();
        Matrix F3 = F + dFdt2 * dt/2;
        Matrix J3 = ptvals.get_J();
        Matrix dFdt3 = J3 * F3;

        intp.locate(x + k3 * dt, t + dt, cell_id);
        intp.evaluate(x + k3 * dt, t + dt, cell_id, ptvals);
        Vector k4 = ptvals.get_u();
        Matrix F4 = F + dFdt3 * dt;
        Matrix J4 = ptvals.get_J();
        Matrix dFdt4 = J4 * F4;

        Vector dx = (k1 + 2*k2 + 2*k3 + k4) * dt/6;
        Matrix dF = (dFdt1 + 2*dFdt2 + 2*dFdt3 + dFdt4) * dt/6;

        if (intp.locate(x + dx, t+dt, cell_id)){
            ++n_accepted_loc;
            particle.x() = x + dx;
            particle.F() = F + dF;
            //particle.n() = el/el.norm();
            //particle.w() += log(el.norm());

            //Matrix J = intp.get_J();
            //particle.S() = particle.n().transpose() * J * particle.n();
            particle.cell_id() = cell_id;
        }
        else {
            //++n_stuck;
            outside_nodes_loc.insert(i);
            ++n_declined_loc;
        }
    }
    #pragma omp critical
    {
        outside_nodes.insert(outside_nodes_loc.begin(), outside_nodes_loc.end());
        n_accepted += n_accepted_loc;
        n_declined += n_declined_loc;
    }
    }
    return outside_nodes;
}

#endif