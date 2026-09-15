#ifndef __RKINTEGRATOR_HPP
#define __RKINTEGRATOR_HPP

#include <algorithm>
#include <vector>
#include "typedefs.hpp"
#include "Interpol.hpp"
#include "Integrator.hpp"
#include "TransportElement.hpp"
#include <omp.h>
#include <math.h>

class RK4Integrator : public Integrator {
public:
    RK4Integrator();
    ~RK4Integrator() {};
  // One loop for all transport elements
  template<TransportElement E = TransportElement::Point, typename Interp>
  std::vector<Uint> step(Interp& intp, ParticleSet& ps, const double t, const double dt);
protected:
};

inline RK4Integrator::RK4Integrator() : Integrator() {
    std::cout << "Selecting Runge-Kutta 4 scheme" << std::endl;
}

template<TransportElement E, typename Interp>
PARTRAC_HOT_LOOP
std::vector<Uint> RK4Integrator::step(Interp& intp, ParticleSet& ps, const double t, const double dt) {
    std::vector<Uint> outside_nodes;
    // Parallel over particles
    #pragma omp parallel
    {
    std::vector<Uint> outside_nodes_loc;
    Uint n_accepted_loc = 0;
    Uint n_declined_loc = 0;
    Vector3d dx, k1, k2, k3, k4;

    #pragma omp for
    for (Uint i=0; i < ps.N(); ++i){
        Vector3d x = ps.x(i);
        CellPos pos;
        pos.id = ps.get_cell_id(i);

        PointValues ptvals(intp.get_U0());

        if constexpr (E == TransportElement::Point){
            intp.locate(x, t, pos);
            intp.evaluate(x, t, pos, ptvals);
            k1 = ptvals.get_u();
            intp.locate(x + k1 * dt/2, t + dt/2, pos);
            intp.evaluate(x + k1 * dt/2, t + dt/2, pos, ptvals);
            k2 = ptvals.get_u();
            intp.locate(x + k2 * dt/2, t + dt/2, pos);
            intp.evaluate(x + k2 * dt/2, t + dt/2, pos, ptvals);
            k3 = ptvals.get_u();
            intp.locate(x + k3 * dt, t + dt, pos);
            intp.evaluate(x + k3 * dt, t + dt, pos, ptvals);
            k4 = ptvals.get_u();

            dx = (k1 + 2*k2 + 2*k3 + k4) * dt/6;

            // Containment only
            if (intp.locate(x+dx, t+dt, pos)){
                ps.set_x(i, x + dx);
                ps.set_t_loc(i, ps.t_loc(i) + dt);
                ps.set_cell_id(i, pos.id);
                ++n_accepted_loc;
            }
            else {
                outside_nodes_loc.push_back(i);
                ++n_declined_loc;
            }
        }
        if constexpr (E == TransportElement::Vector){
            // Line element: d(rhohat)/dt = J rhohat; outside stages contribute zero
            const Vector3d n = ps.rhohat(i);
            k1 = k2 = k3 = k4 = Vector3d::Zero();
            Vector3d F1 = Vector3d::Zero(), F2 = Vector3d::Zero(), F3 = Vector3d::Zero(), F4 = Vector3d::Zero();

            bool is_inside = intp.locate(x, t, pos);
            if (is_inside){
                intp.evaluate(x, t, pos, ptvals);
                k1 = ptvals.get_u();
                Matrix3d J1 = ptvals.get_J();
                F1 = J1 * n;
            }
            is_inside = intp.locate(x + k1 * dt/2, t + dt/2, pos);
            if (is_inside){
                intp.evaluate(x + k1 * dt/2, t + dt/2, pos, ptvals);
                k2 = ptvals.get_u();
                Vector3d n2 = n + F1 * dt/2;
                Matrix3d J2 = ptvals.get_J();
                F2 = J2 * n2;
            }
            is_inside = intp.locate(x + k2 * dt/2, t + dt/2, pos);
            if (is_inside){
                intp.evaluate(x + k2 * dt/2, t + dt/2, pos, ptvals);
                k3 = ptvals.get_u();
                Vector3d n3 = n + F2 * dt/2;
                Matrix3d J3 = ptvals.get_J();
                F3 = J3 * n3;
            }
            is_inside = intp.locate(x + k3 * dt, t + dt, pos);
            if (is_inside){
                intp.evaluate(x + k3 * dt, t + dt, pos, ptvals);
                k4 = ptvals.get_u();
                Vector3d n4 = n + F3 * dt;
                Matrix3d J4 = ptvals.get_J();
                F4 = J4 * n4;
            }
            dx = (k1 + 2*k2 + 2*k3 + k4) * dt/6;
            Vector3d el = n + (F1 + 2*F2 + 2*F3 + F4) * dt/6;

            if (intp.locate(x + dx, t+dt, pos)){
                ps.set_x(i, x + dx);
                ps.set_t_loc(i, ps.t_loc(i) + dt);
                const double len = el.norm();
                ps.set_rhohat(i, el/len);
                ps.set_w(i, ps.w(i) + log(len));
                ps.set_cell_id(i, pos.id);
                ++n_accepted_loc;
            }
            else {
                outside_nodes_loc.push_back(i);
                ++n_declined_loc;
            }
        }
        if constexpr (E == TransportElement::Tensor){
            // Deformation gradient: dF/dt = J F
            const Matrix3d F = ps.F(i);

            intp.locate(x, t, pos);
            intp.evaluate(x, t, pos, ptvals);
            k1 = ptvals.get_u();
            Matrix3d J1 = ptvals.get_J();
            Matrix3d dFdt1 = J1 * F;

            intp.locate(x + k1 * dt/2, t + dt/2, pos);
            intp.evaluate(x + k1 * dt/2, t + dt/2, pos, ptvals);
            k2 = ptvals.get_u();
            Matrix3d F2 = F + dFdt1 * dt/2;
            Matrix3d J2 = ptvals.get_J();
            Matrix3d dFdt2 = J2 * F2;

            intp.locate(x + k2 * dt/2, t + dt/2, pos);
            intp.evaluate(x + k2 * dt/2, t + dt/2, pos, ptvals);
            k3 = ptvals.get_u();
            Matrix3d F3 = F + dFdt2 * dt/2;
            Matrix3d J3 = ptvals.get_J();
            Matrix3d dFdt3 = J3 * F3;

            intp.locate(x + k3 * dt, t + dt, pos);
            intp.evaluate(x + k3 * dt, t + dt, pos, ptvals);
            k4 = ptvals.get_u();
            Matrix3d F4 = F + dFdt3 * dt;
            Matrix3d J4 = ptvals.get_J();
            Matrix3d dFdt4 = J4 * F4;

            dx = (k1 + 2*k2 + 2*k3 + k4) * dt/6;
            Matrix3d dF = (dFdt1 + 2*dFdt2 + 2*dFdt3 + dFdt4) * dt/6;

            if (intp.locate(x + dx, t+dt, pos)){
                ps.set_x(i, x + dx);
                ps.set_t_loc(i, ps.t_loc(i) + dt);
                ps.set_F(i, F + dF);
                ps.set_cell_id(i, pos.id);
                ++n_accepted_loc;
            }
            else {
                outside_nodes_loc.push_back(i);
                ++n_declined_loc;
            }
        }
    }
    #pragma omp critical
    {
        outside_nodes.insert(outside_nodes.end(), outside_nodes_loc.begin(), outside_nodes_loc.end());
        n_accepted += n_accepted_loc;
        n_declined += n_declined_loc;
    }
    }
    // Sort: threads merge out of order
    std::sort(outside_nodes.begin(), outside_nodes.end());
    return outside_nodes;
}

#endif
