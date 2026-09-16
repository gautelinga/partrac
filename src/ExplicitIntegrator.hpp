#ifndef __EXPLICITINTEGRATOR_HPP
#define __EXPLICITINTEGRATOR_HPP

#include <algorithm>
#include <vector>
#include "typedefs.hpp"
#include "Interpol.hpp"
#include "Integrator.hpp"
#include "TransportElement.hpp"
#include <math.h>

class ExplicitIntegrator : public Integrator {
public:
    ExplicitIntegrator(const double Dm, const int int_order, std::vector<std::mt19937>& gens);
    ~ExplicitIntegrator() {};
    // One loop for all transport elements
    template<TransportElement E = TransportElement::Point, typename Interp>
    std::vector<Uint> step(Interp& intp, ParticleSet& ps, const double t, const double dt);
protected:
  double Dm;
  int int_order;
  std::vector<std::mt19937>& gens;
  std::normal_distribution<double> rnd_normal;
};

inline ExplicitIntegrator::ExplicitIntegrator(const double Dm, const int int_order, std::vector<std::mt19937>& gens) : Dm(Dm), int_order(int_order), gens(gens), rnd_normal(0.0, 1.0) {
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

            bool is_inside = intp.locate(x, t, pos);
            Vector3d dx_rw = Vector3d::Zero();
            if (is_inside){
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
            is_inside = intp.locate(x+dx_rw, t+dt, pos);
            if (!is_inside && intp.can_reflect){
                pos.id = ps.get_cell_id(i);
                intp.reflect(x, dx_rw, t, dt, pos.id);
                is_inside = intp.locate(x+dx_rw, t+dt, pos);
            }
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

#endif
