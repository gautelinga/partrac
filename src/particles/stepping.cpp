#include "stepping.hpp"
#include "interpol_dispatch.hpp"

template<TransportElement E>
std::vector<Uint> rk4_step(RK4Integrator& integrator, Interpol& intp, ParticleSet& ps, const double t, const double dt){
  return with_concrete(intp, [&](auto& ip){ return integrator.template step<E>(ip, ps, t, dt); });
}

template<TransportElement E>
std::vector<Uint> explicit_step(ExplicitIntegrator& integrator, Interpol& intp, ParticleSet& ps, const double t, const double dt){
  return with_concrete(intp, [&](auto& ip){ return integrator.template step<E>(ip, ps, t, dt); });
}

template<TransportElement E>
std::vector<Uint> spatial_step(SpatialIntegrator& integrator, Interpol& intp, ParticleSet& ps, const double t, const double ds){
  return with_concrete(intp, [&](auto& ip){ return integrator.template step<E>(ip, ps, t, ds); });
}

void update_fields(ParticleSet& ps, Interpol& intp, const double t, std::map<std::string, bool>& output_fields){
  with_concrete(intp, [&](auto& ip){ ps.update_fields(ip, t, output_fields); });
}

// The elements each step carries
template std::vector<Uint> rk4_step<TransportElement::Point>(RK4Integrator&, Interpol&, ParticleSet&, const double, const double);
template std::vector<Uint> rk4_step<TransportElement::Vector>(RK4Integrator&, Interpol&, ParticleSet&, const double, const double);
template std::vector<Uint> rk4_step<TransportElement::Tensor>(RK4Integrator&, Interpol&, ParticleSet&, const double, const double);
template std::vector<Uint> explicit_step<TransportElement::Point>(ExplicitIntegrator&, Interpol&, ParticleSet&, const double, const double);
template std::vector<Uint> explicit_step<TransportElement::Vector>(ExplicitIntegrator&, Interpol&, ParticleSet&, const double, const double);
template std::vector<Uint> explicit_step<TransportElement::Tensor>(ExplicitIntegrator&, Interpol&, ParticleSet&, const double, const double);
template std::vector<Uint> spatial_step<TransportElement::Point>(SpatialIntegrator&, Interpol&, ParticleSet&, const double, const double);
template std::vector<Uint> spatial_step<TransportElement::Vector>(SpatialIntegrator&, Interpol&, ParticleSet&, const double, const double);
