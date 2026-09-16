#ifndef __SPATIALINTEGRATOR_HPP
#define __SPATIALINTEGRATOR_HPP

#include <algorithm>
#include <vector>
#include "typedefs.hpp"
#include "Interpol.hpp"
#include "Integrator.hpp"
#include "ParticleSet.hpp"
#include "TransportElement.hpp"
#include <omp.h>
#include <math.h>

// Constant path-length steps along streamlines, fields frozen
class SpatialIntegrator : public Integrator {
public:
  SpatialIntegrator(const int int_order, const double u_min, const double dl_max, const double T);
  ~SpatialIntegrator() {};
  template<TransportElement E = TransportElement::Point, typename InterpolType>
  std::vector<Uint> step(InterpolType&, ParticleSet&, double t, double ds);
protected:
  double   m_u_min;
  double   m_dl_max;
  int      m_int_order;
  double   m_T;
};

inline SpatialIntegrator::SpatialIntegrator(const int int_order, const double u_min, const double dl_max, const double T)
  : Integrator(), m_u_min(u_min), m_dl_max(dl_max), m_int_order(int_order), m_T(T) {
    std::cout << "Choosing a spatial integrator of order " << int_order << "." << std::endl;
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

#endif
