#ifndef __EXP_INTEGRATOR_DIRECTIONAL_HPP
#define __EXP_INTEGRATOR_DIRECTIONAL_HPP

#include "typedefs.hpp"
#include "integrator.hpp"
#include <math.h>

class Integrator_Directional : public Integrator {
public:
  Integrator_Directional(const Vector& direction, const int int_order, const Real un_min, const Real dl_max);
  ~Integrator_Directional() {};
  template<typename InterpolType, typename T>
  std::set<Uint> step(InterpolType&, T&, Real t, Real s);
protected:
  Vector m_direction;
  Real   m_un_min;
  Real   m_dl_max;
  int    m_int_order;
};

Integrator_Directional::Integrator_Directional(const Vector& direction, const int int_order, const Real un_min, const Real dl_max)
  : Integrator(), m_direction(direction), m_int_order(int_order), m_un_min(un_min), m_dl_max(dl_max) {
    std::cout << "Choosing a directional integrator of order " << int_order << "." << std::endl;
}

template<typename InterpolType, typename T>
std::set<Uint> Integrator_Directional::step(InterpolType& intp, T& ps, const Real t, const Real s) {
    std::set<Uint> outside_nodes;
    Uint i = 0;
    bool is_inside;
    double s_prev, un_est, dt;
    Vector dx;

    for (auto & particle : ps.particles() ){
        Vector x = particle.x();

        int cell_id = -1;
        PointValues ptvals(intp.get_U0());
        intp.locate(x, t, cell_id);
        intp.evaluate(x, t, cell_id, ptvals);

        Vector u_1 = ptvals.get_u();

        un_est = u_1.dot(m_direction);

        is_inside = false;
        if (un_est > m_un_min){
            s_prev = x.dot(m_direction);
            dt = (s - s_prev) / un_est;

            dx = u_1 * dt;

            // Second-order terms
            if (m_int_order >= 2){
                dx += 0.5 * (ptvals.get_a() + ptvals.get_Ju()) * dt * dt;
            }

            if (dx.norm() < m_dl_max){
                // Frozen time, otherwise: locate(x+dx, t+dt, cell_id)
                is_inside = intp.locate(x + dx, t, cell_id);
            }
            else {
                std::cout << "Step too long (dl=" << dx.norm() << "), consider doing something smart!" << std::endl;
            }
        }
        // count things
        if (is_inside){
            ++n_accepted;
            particle.x() = x + dx;
            particle.tau() += dt;
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