#ifndef __STEPPING_HPP
#define __STEPPING_HPP

// Steps on the abstract interpolator, dispatched to the loops compiled for each concrete one

#include <map>
#include <string>
#include <vector>
#include "typedefs.hpp"
#include "Interpol.hpp"
#include "ParticleSet.hpp"
#include "RKIntegrator.hpp"
#include "ExplicitIntegrator.hpp"
#include "SpatialIntegrator.hpp"
#include "TransportElement.hpp"

template<TransportElement E>
std::vector<Uint> rk4_step(RK4Integrator& integrator, Interpol& intp, ParticleSet& ps, const double t, const double dt);

template<TransportElement E>
std::vector<Uint> explicit_step(ExplicitIntegrator& integrator, Interpol& intp, ParticleSet& ps, const double t, const double dt);

// Fields frozen at t, path length ds
template<TransportElement E>
std::vector<Uint> spatial_step(SpatialIntegrator& integrator, Interpol& intp, ParticleSet& ps, const double t, const double ds);

// Fields at the particles
void update_fields(ParticleSet& ps, Interpol& intp, const double t, std::map<std::string, bool>& output_fields);

#endif
