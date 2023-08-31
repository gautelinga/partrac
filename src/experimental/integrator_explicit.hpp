#ifndef __EXP_INTEGRATOR_EXPLICIT_H
#define __EXP_INTEGRATOR_EXPLICIT_H

#include "typedefs.hpp"
#include "integrator.hpp"
#include <math.h>

class Integrator_Explicit : public Integrator {
public:
  Integrator_Explicit(const Real Dm, const int int_order, std::vector<std::mt19937>& gens);
  ~Integrator_Explicit() {};
  template<typename InterpolType, typename T>
  std::set<Uint> step(InterpolType&, T&, Real t, Real dt);
  template<typename InterpolType, typename T>
  void step_parallel(InterpolType&, T&, Real t, Real dt);
protected:
  Real Dm;
  int int_order;
  std::vector<std::mt19937>& gens;
  std::normal_distribution<Real> rnd_normal;
};

Integrator_Explicit::Integrator_Explicit(const Real Dm, const int int_order, std::vector<std::mt19937>& gens)
  : Integrator(), gens(gens), Dm(Dm), int_order(int_order), rnd_normal(0.0, 1.0) {
    std::cout << "Choosing an explicit integrator of order " << int_order << " with diffusivity " << Dm << std::endl;
}

template<typename InterpolType, typename T>
std::set<Uint> Integrator_Explicit::step(InterpolType& intp, T& ps, const Real t, const Real dt) {
    std::set<Uint> outside_nodes;
    Real sqrt2Dmdt = sqrt(2*Dm*dt);
    Uint i = 0;
    for (auto & particle : ps.particles() ){
        Vector x = particle.x();
        int cell_id = particle.cell_id();

        intp.probe(x, t, cell_id);
        Vector dx = intp.get_u() * dt;

        // Second-order terms
        if (int_order >= 2){
            //dx_rw += 0.5*a_rw[irw]*dt2;
            dx += 0.5 * (intp.get_a() + intp.get_Ju()) * dt * dt;
        }
        if (Dm > 0.0){
            std::mt19937& gen = gens[0]; // Not parallel yet
            Vector eta = {rnd_normal(gen),
                          rnd_normal(gen),
                          rnd_normal(gen)};
            dx += sqrt2Dmdt * eta;
        }
        intp.probe(x+dx, t+dt, cell_id);
        if (intp.inside_domain()){
            ++n_accepted;
            particle.x() = x + dx;
            particle.cell_id() = cell_id;
        }
        else {
            outside_nodes.insert(i);
            ++n_declined;
        }
        ++i;
    }
    return outside_nodes;
}

template<typename InterpolType, typename T>
void Integrator_Explicit::step_parallel(InterpolType& intp, T& ps, const Real t, const Real dt) {
    std::set<Uint> outside_nodes;
    Real sqrt2Dmdt = sqrt(2*Dm*dt);
    // Uint i = 0;
    double U0 = intp.get_U0();

    #pragma omp parallel 
    {
        std::normal_distribution<Real> _rnd_normal(0., 1.0);
        std::mt19937& gen = gens[omp_get_thread_num()];

        # pragma omp for
        for (auto & particle : ps.particles() ){
            Vector3d x = particle.x();
            int cell_id = particle.cell_id();
            PointValues ptvals(U0);

            bool is_inside = intp.probe_light(x, t, cell_id);
            intp.probe_heavy(x, t, cell_id, ptvals);

            Vector3d dx = ptvals.get_u() * dt;

            // Second-order terms
            if (int_order >= 2) dx += 0.5*(ptvals.get_Ju() + ptvals.get_a()) * dt * dt;
            if (Dm > 0.0){
                Vector eta = {_rnd_normal(gen), _rnd_normal(gen), _rnd_normal(gen)};
                dx += sqrt2Dmdt * eta;
            }
            is_inside = intp.probe_light(x+dx, t+dt, cell_id);
            if (is_inside){
                //++n_accepted;
                particle.x() = x + dx;
                particle.cell_id() = cell_id;
            }
        }
    }
}

#endif