#ifndef __EXPLICITINTEGRATOR_HPP
#define __EXPLICITINTEGRATOR_HPP

#include "typedefs.hpp"
#include "Interpol.hpp"
#include "Integrator.hpp"
#include <math.h>

class ExplicitIntegrator : public Integrator {
public:
    //ExplicitIntegrator(std::shared_ptr<Interpol> intp, const double Dm, const int int_order, std::mt19937& gen);
    ExplicitIntegrator(const double Dm, const int int_order, std::vector<std::mt19937>& gens);
    ~ExplicitIntegrator() {};
    // Vector3d integrate(const Vector3d& x, const double t, const double dt);
    //template<typename InterpolType, typename T>
    //std::set<Uint> step(InterpolType& intp, T& ps, const double t, const double dt);
    std::set<Uint> step(ParticleSet& ps, const double t, const double dt);
    // the same loop with the interpolator known by type; see interpol_dispatch.hpp
    template<typename Interp>
    std::set<Uint> step(Interp& intp, ParticleSet& ps, const double t, const double dt);
protected:
  double Dm;
  int int_order;
  //std::mt19937& gen;
  std::vector<std::mt19937>& gens;
  std::normal_distribution<double> rnd_normal;
};

/*
ExplicitIntegrator::ExplicitIntegrator(std::shared_ptr<Interpol> intp, const double Dm, const int int_order, std::mt19937& gen) : Integrator(intp), gen(gen), Dm(Dm), int_order(int_order), rnd_normal(0.0, 1.0) {
}*/
ExplicitIntegrator::ExplicitIntegrator(const double Dm, const int int_order, std::vector<std::mt19937>& gens) : Dm(Dm), int_order(int_order), gens(gens), rnd_normal(0.0, 1.0) {
}


//template<typename InterpolType, typename T>
//std::set<Uint> ExplicitIntegrator::step(InterpolType& intp, T& ps, const double t, const double dt) {
std::set<Uint> ExplicitIntegrator::step(ParticleSet& ps, const double t, const double dt) {
    return step(*ps.interpolator(), ps, t, dt);
}

template<typename Interp>
std::set<Uint> ExplicitIntegrator::step(Interp& intp, ParticleSet& ps, const double t, const double dt) {
    std::set<Uint> outside_nodes;

    #pragma omp parallel 
    {
        std::normal_distribution<double> _rnd_normal(0., 1.0);
        std::mt19937& gen = gens[omp_get_thread_num()];
        std::set<Uint> outside_nodes_loc;
        Uint n_accepted_loc = 0;
        Uint n_declined_loc = 0;

        double sqrt2Dmdt = sqrt(2*Dm*dt);

        # pragma omp for
        for (Uint i=0; i < ps.N(); ++i){
            Vector3d x = ps.x(i);
            int cell_id = ps.get_cell_id(i);

            PointValues ptvals(intp.get_U0());

            bool is_inside = intp.locate(x, t, cell_id);
            intp.evaluate(x, t, cell_id, ptvals);

            //Vector3d dx_rw = intp.get_u() * dt;
            Vector3d dx_rw = ptvals.get_u() * dt;

            // Second-order terms
            if (int_order >= 2){
                //dx_rw += 0.5*a_rw[irw]*dt2;
                //dx_rw += 0.5 * (intp.get_a() + intp.get_Ju()) * dt * dt;
                dx_rw += 0.5 * (ptvals.get_a() + ptvals.get_Ju()) * dt * dt;
            }
            if (Dm > 0.0){
                // TODO: Consider trying multiple times
                Vector3d eta = {_rnd_normal(gen),
                                _rnd_normal(gen),
                                _rnd_normal(gen)};
                dx_rw += sqrt2Dmdt * eta;
            }
            is_inside = intp.locate(x+dx_rw, t+dt, cell_id);
            if (!is_inside && intp.can_reflect){
                cell_id = ps.get_cell_id(i);
                intp.reflect(x, dx_rw, t, dt, cell_id);
                is_inside = intp.locate(x+dx_rw, t+dt, cell_id);
            }
            if (is_inside){
                ps.set_x(i, x + dx_rw);
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