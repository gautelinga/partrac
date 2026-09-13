#ifndef __RKINTEGRATOR_HPP
#define __RKINTEGRATOR_HPP

#include "typedefs.hpp"
#include "Interpol.hpp"
#include "Integrator.hpp"
#include <omp.h>
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
  template<typename Interp>
  std::set<Uint> step(Interp& intp, ParticleSet& ps, const double t, const double dt);
protected:
};

//RK4Integrator::RK4Integrator(std::shared_ptr<Interpol> intp) : Integrator(intp) {
RK4Integrator::RK4Integrator() : Integrator() {
    std::cout << "Selecting Runge-Kutta 4 scheme" << std::endl;
}

//template<typename InterpolType, typename T>
//std::set<Uint> RK4Integrator::step(InterpolType& intp, T& ps, const double t, const double dt) {
std::set<Uint> RK4Integrator::step(ParticleSet& ps, const double t, const double dt) {
    return step(*ps.interpolator(), ps, t, dt);
}

template<typename Interp>
std::set<Uint> RK4Integrator::step(Interp& intp, ParticleSet& ps, const double t, const double dt) {
    std::set<Uint> outside_nodes;
    // Particles are independent; only the tally is shared
    #pragma omp parallel
    {
    std::set<Uint> outside_nodes_loc;
    Uint n_accepted_loc = 0;
    Uint n_declined_loc = 0;
    Vector3d dx, k1, k2, k3, k4;

    #pragma omp for
    for (Uint i=0; i < ps.N(); ++i){
        Vector3d x = ps.x(i);
        int cell_id = ps.get_cell_id(i);

        PointValues ptvals(intp.get_U0());

        intp.locate(x, t, cell_id);
        intp.evaluate(x, t, cell_id, ptvals);
        k1 = ptvals.get_u();
        intp.locate(x + k1 * dt/2, t + dt/2, cell_id);
        intp.evaluate(x + k1 * dt/2, t + dt/2, cell_id, ptvals);
        k2 = ptvals.get_u();
        intp.locate(x + k2 * dt/2, t + dt/2, cell_id);
        intp.evaluate(x + k2 * dt/2, t + dt/2, cell_id, ptvals);
        k3 = ptvals.get_u();
        intp.locate(x + k3 * dt, t + dt, cell_id);
        intp.evaluate(x + k3 * dt, t + dt, cell_id, ptvals);
        k4 = ptvals.get_u();

        dx = (k1 + 2*k2 + 2*k3 + k4) * dt/6;

        // only containment is needed here, so no evaluate
        if (intp.locate(x+dx, t+dt, cell_id)){
            ps.set_x(i, x + dx);
            ps.set_t_loc(i, ps.t_loc(i) + dt);
            ps.set_cell_id(i, cell_id);
            ++n_accepted_loc;
        }
        else {
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