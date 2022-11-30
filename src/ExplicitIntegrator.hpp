#ifndef __EXPLICITINTEGRATOR_HPP
#define __EXPLICITINTEGRATOR_HPP

#include "typedefs.hpp"
#include "Interpol.hpp"
#include "Integrator.hpp"
#include <math.h>

class ExplicitIntegrator : public Integrator {
public:
    //ExplicitIntegrator(std::shared_ptr<Interpol> intp, const double Dm, const int int_order, std::mt19937& gen);
    ExplicitIntegrator(const double Dm, const int int_order, std::mt19937& gen);
    ~ExplicitIntegrator() {};
    // Vector3d integrate(const Vector3d& x, const double t, const double dt);
    //template<typename InterpolType, typename T>
    //std::set<Uint> step(InterpolType& intp, T& ps, const double t, const double dt);
    std::set<Uint> step(ParticleSet& ps, const double t, const double dt);
protected:
  double Dm;
  int int_order;
  std::mt19937& gen;
  std::normal_distribution<double> rnd_normal;
};

/*
ExplicitIntegrator::ExplicitIntegrator(std::shared_ptr<Interpol> intp, const double Dm, const int int_order, std::mt19937& gen) : Integrator(intp), gen(gen), Dm(Dm), int_order(int_order), rnd_normal(0.0, 1.0) {
}*/
ExplicitIntegrator::ExplicitIntegrator(const double Dm, const int int_order, std::mt19937& gen) : gen(gen), Dm(Dm), int_order(int_order), rnd_normal(0.0, 1.0) {
}

/*
Vector3d ExplicitIntegrator::integrate(const Vector3d& x, const double t, const double dt) {
    intp->probe(x, t);
    Vector3d dx_rw = intp->get_u() * dt;

    // Set elongation
    //ps.e_rw[irw] = 0.;

    // Second-order terms
    if (int_order >= 2){
        //dx_rw += 0.5*a_rw[irw]*dt2;
        dx_rw += 0.5 * (intp->get_a() + intp->get_Ju()) * dt * dt;
    }
    if (Dm > 0.0){
        Vector3d eta = {rnd_normal(gen),
                        rnd_normal(gen),
                        rnd_normal(gen)};
        double sqrt2Dmdt = sqrt(2*Dm*dt);
        dx_rw += sqrt2Dmdt * eta;
    }
    intp->probe(x+dx_rw, t+dt);
    if (intp->inside_domain()){
        stuck = false;
        ++n_accepted;
        return dx_rw;
    }
    else {
        ++n_declined;
        stuck = true;
        return {0., 0., 0.};
    }
}*/

//template<typename InterpolType, typename T>
//std::set<Uint> ExplicitIntegrator::step(InterpolType& intp, T& ps, const double t, const double dt) {
std::set<Uint> ExplicitIntegrator::step(ParticleSet& ps, const double t, const double dt) {
    std::set<Uint> outside_nodes;
    double sqrt2Dmdt = sqrt(2*Dm*dt);

    auto & intp = *ps.interpolator();

    for (Uint i=0; i < ps.N(); ++i){
        Vector3d x = ps.x(i);
        int cell_id = ps.get_cell_id(i);

        intp.probe(x, t, cell_id);
        
        Vector3d dx_rw = intp.get_u() * dt;

        // Second-order terms
        if (int_order >= 2){
            //dx_rw += 0.5*a_rw[irw]*dt2;
            dx_rw += 0.5 * (intp.get_a() + intp.get_Ju()) * dt * dt;
        }
        if (Dm > 0.0){
            Vector3d eta = {rnd_normal(gen),
                            rnd_normal(gen),
                            rnd_normal(gen)};      
            dx_rw += sqrt2Dmdt * eta;
        }
        intp.probe(x+dx_rw, t+dt, cell_id);
        if (intp.inside_domain()){
            ps.set_x(i, x + dx_rw);
            ps.set_t_loc(i, ps.t_loc(i) + dt);
            ps.set_cell_id(i, cell_id);
            ++n_accepted;
        }
        else {
            outside_nodes.insert(i);
            ++n_declined;
        }
    }
    return outside_nodes;
}

#endif