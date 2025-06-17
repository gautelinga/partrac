#ifndef __EXP_INTEGRATOR_SPATIAL_HPP
#define __EXP_INTEGRATOR_SPATIAL_HPP

#include "typedefs.hpp"
#include "integrator.hpp"
#include <math.h>

class Integrator_Spatial : public Integrator {
public:
  Integrator_Spatial(const int int_order, const Real un_min, const Real dl_max);
  ~Integrator_Spatial() {};
  template<typename InterpolType, typename T>
  std::set<Uint> step_vec(InterpolType&, T&, Real t, Real ds);
protected:
  Real   m_un_min;
  Real   m_dl_max;
  int    m_int_order;
};

Integrator_Spatial::Integrator_Spatial(const int int_order, const Real un_min, const Real dl_max)
  : Integrator(), m_int_order(int_order), m_un_min(un_min), m_dl_max(dl_max) {
    std::cout << "Choosing a spatial integrator of order " << int_order << "." << std::endl;
}

template<typename InterpolType, typename T>
std::set<Uint> Integrator_Spatial::step_vec(InterpolType& intp, T& ps, const Real t, const Real ds) {
    std::set<Uint> outside_nodes;
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

            Vector dx;
            Vector el = n;
            double dt;

            PointValues ptvals(intp.get_U0());

            bool is_inside = intp.probe_light(x, t, cell_id);
            if (is_inside)
            {
                intp.probe_heavy(x, t, cell_id, ptvals);
                Vector u1 = ptvals.get_u();
                Matrix J1 = ptvals.get_J();

                double u0 = std::max(u1.norm(), m_un_min);

                dt = ds/u0;
                dx = u1 * dt;
                el += J1 * n * dt;

                if (m_int_order >= 2){
                    dx += 0.5*(J1*u1 + ptvals.get_a()) * dt * dt;
                    el += 0.5*(J1*(J1*n) + ptvals.get_grada()*n) * dt * dt;
                }
                is_inside = intp.probe_light(x + dx, t+dt, cell_id);
                if (is_inside){
                    ++n_accepted_loc;
                    particle.x() = x + dx;
                    particle.n() = el/el.norm();
                    particle.w() += log(el.norm());

                    intp.probe_heavy(x + dx, t+dt, cell_id, ptvals);
                    Matrix J = ptvals.get_J();
                    particle.S() = particle.n().transpose() * J * particle.n();
                    particle.cell_id() = cell_id;
                    particle.tau() += dt;
                }
                else {
                    //++n_stuck;
                    outside_nodes_loc.insert(i);
                    ++n_declined_loc;
                }
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

#endif