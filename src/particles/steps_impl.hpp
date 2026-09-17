#ifndef __STEPS_IMPL_HPP
#define __STEPS_IMPL_HPP

// Loop bodies of the steps; included only by the per-interpolator step sources

#include <algorithm>
#include <cmath>
#include <vector>
#include <omp.h>
#include "ParticleSet.hpp"
#include "RKIntegrator.hpp"
#include "ExplicitIntegrator.hpp"
#include "SpatialIntegrator.hpp"

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
            // Outside stages contribute zero
            k1 = k2 = k3 = k4 = Vector3d::Zero();
            if (intp.locate(x, t, pos)){
                intp.evaluate(x, t, pos, ptvals);
                k1 = ptvals.get_u();
            }
            if (intp.locate(x + k1 * dt/2, t + dt/2, pos)){
                intp.evaluate(x + k1 * dt/2, t + dt/2, pos, ptvals);
                k2 = ptvals.get_u();
            }
            if (intp.locate(x + k2 * dt/2, t + dt/2, pos)){
                intp.evaluate(x + k2 * dt/2, t + dt/2, pos, ptvals);
                k3 = ptvals.get_u();
            }
            if (intp.locate(x + k3 * dt, t + dt, pos)){
                intp.evaluate(x + k3 * dt, t + dt, pos, ptvals);
                k4 = ptvals.get_u();
            }

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
            // Deformation gradient: dF/dt = J F; outside stages contribute zero
            const Matrix3d F = ps.F(i);
            k1 = k2 = k3 = k4 = Vector3d::Zero();
            Matrix3d dFdt1 = Matrix3d::Zero(), dFdt2 = Matrix3d::Zero(), dFdt3 = Matrix3d::Zero(), dFdt4 = Matrix3d::Zero();

            if (intp.locate(x, t, pos)){
                intp.evaluate(x, t, pos, ptvals);
                k1 = ptvals.get_u();
                Matrix3d J1 = ptvals.get_J();
                dFdt1 = J1 * F;
            }
            if (intp.locate(x + k1 * dt/2, t + dt/2, pos)){
                intp.evaluate(x + k1 * dt/2, t + dt/2, pos, ptvals);
                k2 = ptvals.get_u();
                Matrix3d F2 = F + dFdt1 * dt/2;
                Matrix3d J2 = ptvals.get_J();
                dFdt2 = J2 * F2;
            }
            if (intp.locate(x + k2 * dt/2, t + dt/2, pos)){
                intp.evaluate(x + k2 * dt/2, t + dt/2, pos, ptvals);
                k3 = ptvals.get_u();
                Matrix3d F3 = F + dFdt2 * dt/2;
                Matrix3d J3 = ptvals.get_J();
                dFdt3 = J3 * F3;
            }
            if (intp.locate(x + k3 * dt, t + dt, pos)){
                intp.evaluate(x + k3 * dt, t + dt, pos, ptvals);
                k4 = ptvals.get_u();
                Matrix3d F4 = F + dFdt3 * dt;
                Matrix3d J4 = ptvals.get_J();
                dFdt4 = J4 * F4;
            }

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

template<TransportElement E, typename Interp>
PARTRAC_HOT_LOOP
std::vector<Uint> ExplicitIntegrator::step(Interp& intp, ParticleSet& ps, const double t, const double dt) {
    std::vector<Uint> outside_nodes;

    #pragma omp parallel
    {
        std::normal_distribution<double> _rnd_normal(0., 1.0);
        std::mt19937& gen = gens[omp_get_thread_num()];
        std::vector<Uint> outside_nodes_loc;
        Uint n_accepted_loc = 0;
        Uint n_declined_loc = 0;

        double sqrt2Dmdt = sqrt(2*Dm*dt);
        const bool reflects = Dm > 0.0 && intp.can_reflect;

        # pragma omp for
        for (Uint i=0; i < ps.N(); ++i){
            Vector3d x = ps.x(i);
            CellPos pos;
            pos.id = ps.get_cell_id(i);

            PointValues ptvals(intp.get_U0());

            // Carried element
            [[maybe_unused]] Vector3d n0, el;
            [[maybe_unused]] Matrix3d F0, Fel;
            if constexpr (E == TransportElement::Vector){ n0 = ps.rhohat(i); el = n0; }
            if constexpr (E == TransportElement::Tensor){ F0 = ps.F(i); Fel = F0; }

            const bool x_inside = intp.locate(x, t, pos);
            Vector3d dx_rw = Vector3d::Zero();
            if (x_inside){
                // Outside: zero velocity
                intp.evaluate(x, t, pos, ptvals);
                dx_rw = ptvals.get_u() * dt;
                if (int_order >= 2){
                    dx_rw += 0.5 * (ptvals.get_a() + ptvals.get_Ju()) * dt * dt;
                }
                // The gradient, for what carries it
                if constexpr (E != TransportElement::Point){
                    const Matrix3d J1 = ptvals.get_J();
                    if constexpr (E == TransportElement::Vector){
                        el += J1 * el * dt;
                        if (int_order >= 2)
                            el += 0.5*(J1*(J1*n0) + ptvals.get_grada()*n0) * dt * dt;
                    }
                    if constexpr (E == TransportElement::Tensor){
                        Fel += J1 * Fel * dt;
                        if (int_order >= 2)
                            Fel += 0.5*(J1*(J1*F0) + ptvals.get_grada()*F0) * dt * dt;
                    }
                }
            }
            if (Dm > 0.0){
                // TODO: Consider trying multiple times
                Vector3d eta = {_rnd_normal(gen),
                                _rnd_normal(gen),
                                _rnd_normal(gen)};
                dx_rw += sqrt2Dmdt * eta;
            }
            // Diffusive steps walk off the walls
            const bool is_inside = (reflects && x_inside)
                ? intp.reflect(x, dx_rw, pos)
                : intp.locate(x+dx_rw, t+dt, pos);
            if (is_inside){
                ps.set_x(i, x + dx_rw);
                ps.set_t_loc(i, ps.t_loc(i) + dt);
                ps.set_cell_id(i, pos.id);
                if constexpr (E == TransportElement::Vector){
                    const double len = el.norm();
                    ps.set_rhohat(i, el/len);
                    ps.set_w(i, ps.w(i) + log(len));
                }
                if constexpr (E == TransportElement::Tensor)
                    ps.set_F(i, Fel);
                ++n_accepted_loc;
            }
            else {
                outside_nodes_loc.push_back(i);
                ++n_declined_loc;
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

template<TransportElement E, typename InterpolType>
PARTRAC_HOT_LOOP
std::vector<Uint> SpatialIntegrator::step(InterpolType& intp, ParticleSet& ps, const double t, const double ds) {
    std::vector<Uint> outside_nodes;
    // Nodes are independent; only the tally is shared
    #pragma omp parallel
    {
    std::vector<Uint> outside_nodes_loc;
    Uint n_accepted_loc = 0;
    Uint n_declined_loc = 0;
    bool is_inside;
    double uabs_est, dt;
    Vector3d dx;
    Vector3d el;
    Matrix3d F;

    #pragma omp for
    for (Uint i=0; i < ps.N(); ++i){
        // Done
        if (ps.t_loc(i) >= m_T){
            outside_nodes_loc.push_back(i);
            ++n_declined_loc;
            continue;
        }
        Vector3d x = ps.x(i);
        CellPos pos;
        pos.id = ps.get_cell_id(i);

        PointValues ptvals(intp.get_U0());
        // Outside: zero velocity, so declined below
        if (intp.locate(x, t, pos))
            intp.evaluate(x, t, pos, ptvals);

        Vector3d u_1 = ptvals.get_u();

        uabs_est = u_1.norm();

        is_inside = false;
        if (uabs_est > m_u_min){
            dt = ds / uabs_est;

            dx = u_1 * dt;

            // Second-order terms
            if (m_int_order >= 2){
                dx += 0.5 * (ptvals.get_a() + ptvals.get_Ju()) * dt * dt;
            }
            if constexpr (E == TransportElement::Vector){
                const Vector3d n = ps.rhohat(i);
                const Matrix3d J1 = ptvals.get_J();
                el = n + J1 * n * dt;
                if (m_int_order >= 2)
                    el += 0.5 * (J1 * (J1 * n) + ptvals.get_grada() * n) * dt * dt;
            }
            if constexpr (E == TransportElement::Tensor){
                const Matrix3d F0 = ps.F(i);
                const Matrix3d J1 = ptvals.get_J();
                F = F0 + J1 * F0 * dt;
                if (m_int_order >= 2)
                    F += 0.5 * (J1 * (J1 * F0) + ptvals.get_grada() * F0) * dt * dt;
            }

            if (dx.norm() < m_dl_max){
                // Frozen time, otherwise: locate(x+dx, t+dt, pos)
                is_inside = intp.locate(x + dx, t, pos);
            }
            else {
                #pragma omp critical
                std::cout << "Step too long (dl=" << dx.norm() << "), consider doing something smart!" << std::endl;
            }
        }
        // count things
        if (is_inside){
            ++n_accepted_loc;
            ps.set_x(i, x + dx);
            ps.set_t_loc(i, ps.t_loc(i) + dt);
            ps.set_cell_id(i, pos.id);
            if constexpr (E == TransportElement::Vector){
                const double len = el.norm();
                ps.set_rhohat(i, el/len);
                ps.set_w(i, ps.w(i) + log(len));
            }
            if constexpr (E == TransportElement::Tensor)
                ps.set_F(i, F);
        }
        else {
            outside_nodes_loc.push_back(i);
            ++n_declined_loc;
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

template<typename Interp>
PARTRAC_HOT_LOOP
void ParticleSet::update_fields(Interp& intp, const double t, std::map<std::string, bool> &output_fields){
  // Read before the threads: operator[] inserts
  const bool do_rho = output_fields["rho"];
  const bool do_p = output_fields["p"];
  const bool do_J = has_J && output_fields["J"];
  // phi always updated (vector statistics)
  const bool do_phi = has_phi;
  const bool do_cell_type = has_cell_type;
  const bool do_S = element == TransportElement::Vector;

  #pragma omp parallel for
  for (Uint irw=0; irw < N(); ++irw){
    CellPos pos;
    pos.id = get_cell_id(irw);
    PointValues ptvals(intp.get_U0());
    // Outside: zero fields
    if (intp.locate(x_rw[irw], t, pos))
      intp.evaluate(x_rw[irw], t, pos, ptvals);
    // Stretching rate
    if (do_S){
      const Matrix3d J = ptvals.get_J();
      S_rw[irw] = rhohat_rw[irw].transpose() * J * rhohat_rw[irw];
    }
    // Always, for the statistics
    u_rw[irw] = ptvals.get_u();
    if (do_rho){
      //rho_rw[irw] = intp->get_rho();
      rho_rw[irw] = ptvals.get_rho();
    }
    if (do_p){
      //p_rw[irw] = intp->get_p();
      p_rw[irw] = ptvals.get_p();
    }
    if (do_J)
      J_rw[irw] = ptvals.get_J();
    if (do_phi)
      phi_rw[irw] = ptvals.get_phi();
    if (do_cell_type)
      cell_type_rw[irw] = ptvals.get_cell_type();
  }
}

#endif
