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
#include "tableau.hpp"
#include <type_traits>
#include <utility>

// A step reads velocity, acceleration and their gradients: an interpolator
// with an evaluate_motion leaves the rest out; called on the concrete type,
// so never through the vtable
template<typename Interp, typename = void>
struct has_evaluate_motion : std::false_type {};
template<typename Interp>
struct has_evaluate_motion<Interp, std::void_t<decltype(std::declval<Interp&>().evaluate_motion(
  std::declval<const Vector3d&>(), 0., std::declval<const CellPos&>(), std::declval<PointValues&>()))>> : std::true_type {};

template<typename Interp>
inline void evaluate_motion(Interp& intp, const Vector3d& x, const double t, const CellPos& pos, PointValues& ptvals){
  if constexpr (has_evaluate_motion<Interp>::value) intp.evaluate_motion(x, t, pos, ptvals);
  else intp.evaluate(x, t, pos, ptvals);
}

// A stage located from scratch: inside if a cell holds it
template<typename Interp>
struct LocatedEval {
  Interp& intp;
  PointValues& ptvals;
  CellPos& pos;
  bool operator()(const Vector3d& x, const double t){
    if (!intp.locate(x, t, pos))
      return false;
    evaluate_motion(intp, x, t, pos, ptvals);
    return true;
  }
  Vector3d u(){ return ptvals.get_u(); }
  Matrix3d J(){ return ptvals.get_J(); }
  bool end(const Vector3d& x, const double t){ return intp.locate(x, t, pos); }
};

// A step with a stage or its end outside, again in 2, 4, then 8 substeps with
// every stage inside; false: the particle leaves the fluid. Out of the loop:
// rare and bulky
template<TransportElement E, typename Interp>
__attribute__((noinline))
bool rk4_substeps(Interp& intp, ParticleSet& ps, const Uint i, const double t, const double dt){
  for (int m = 2; m <= 8; m *= 2){
    const double h = dt/m;
    Vector3d x = ps.x(i);
    CellPos pos;
    pos.id = ps.get_cell_id(i);
    Vector3d n = Vector3d::Zero();
    double w = 0.;
    Matrix3d F = Matrix3d::Identity();
    if constexpr (E == TransportElement::Vector){ n = ps.rhohat(i); w = ps.w(i); }
    if constexpr (E == TransportElement::Tensor) F = ps.frame(i);
    PointValues ptvals(intp.get_U0());
    LocatedEval<Interp> ev{intp, ptvals, pos};
    bool ok = true;
    for (int s = 0; s < m && ok; ++s){
      Vector3d dx, el;
      Matrix3d dF;
      ok = rk_stages<RK4Tableau, E, true>(ev, x, n, F, t + s*h, h, dx, el, dF);
      if (!ok)
        break;
      if constexpr (E == TransportElement::Vector){
        const double len = el.norm();
        n = el/len;
        w += log(len);
      }
      if constexpr (E == TransportElement::Tensor) F += dF;
      x += dx;
    }
    if (!ok)
      continue;
    ps.set_x(i, x);
    ps.set_t_loc(i, ps.t_loc(i) + dt);
    ps.set_cell_id(i, pos.id);
    if constexpr (E == TransportElement::Vector){ ps.set_rhohat(i, n); ps.set_w(i, w); }
    if constexpr (E == TransportElement::Tensor) ps.advance_frame(i, F);
    return true;
  }
  return false;
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

    #pragma omp for
    for (Uint i=0; i < ps.N(); ++i){
        Vector3d x = ps.x(i);
        CellPos pos;
        pos.id = ps.get_cell_id(i);
        PointValues ptvals(intp.get_U0());
        // Carried element
        [[maybe_unused]] const Vector3d n = E == TransportElement::Vector ? ps.rhohat(i) : Vector3d::Zero();
        [[maybe_unused]] const Matrix3d F = E == TransportElement::Tensor ? ps.frame(i) : Matrix3d::Zero();
        Vector3d dx;
        [[maybe_unused]] Vector3d el;
        [[maybe_unused]] Matrix3d dF;

        LocatedEval<Interp> ev{intp, ptvals, pos};
        // An outside stage or end: substeps
        if (rk_stages<RK4Tableau, E, false>(ev, x, n, F, t, dt, dx, el, dF)){
            ps.set_x(i, x + dx);
            ps.set_t_loc(i, ps.t_loc(i) + dt);
            if constexpr (E == TransportElement::Vector){
                const double len = el.norm();
                ps.set_rhohat(i, el/len);
                ps.set_w(i, ps.w(i) + log(len));
            }
            if constexpr (E == TransportElement::Tensor)
                ps.advance_frame(i, F + dF);
            ps.set_cell_id(i, pos.id);
            ++n_accepted_loc;
        }
        else if (rk4_substeps<E>(intp, ps, i, t, dt))
            ++n_accepted_loc;
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
            if constexpr (E == TransportElement::Tensor){ F0 = ps.frame(i); Fel = F0; }

            const bool x_inside = intp.locate(x, t, pos);
            Vector3d dx_rw = Vector3d::Zero();
            if (x_inside){
                // Outside: zero velocity
                evaluate_motion(intp, x, t, pos, ptvals);
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
                    ps.advance_frame(i, Fel);
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
            evaluate_motion(intp, x, t, pos, ptvals);

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
                const Matrix3d F0 = ps.frame(i);
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
                ps.advance_frame(i, F);
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
void ParticleSet::update_fields(Interp& intp, const double t, const OutputFields& output_fields){
  const bool do_rho = output_fields.rho;
  const bool do_p = output_fields.p;
  const bool do_J = has_J && output_fields.J;
  // phi always updated (vector statistics)
  const bool do_phi = has_phi;
  const bool do_cell_type = has_cell_type;
  const bool do_S = element == TransportElement::Vector;
  const bool do_S3 = element == TransportElement::Tensor;

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
    // Stretching rates along the frame
    if (do_S3){
      Matrix3d Q;
      Vector3d s, U;
      settled(irw, Q, s, U);
      S3_rw[irw] = (Q.transpose() * ptvals.get_J() * Q).diagonal();
    }
    // Always, for the statistics
    u_rw[irw] = ptvals.get_u();
    if (do_rho){
      rho_rw[irw] = ptvals.get_rho();
    }
    if (do_p){
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
