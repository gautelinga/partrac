#ifndef __STEPPING_HPP
#define __STEPPING_HPP

// Steps on the abstract interpolator, dispatched to the loops compiled for each concrete one

#include <string>
#include <vector>
#include "typedefs.hpp"
#include "Interpol.hpp"
#include "ParticleSet.hpp"
#include "RKIntegrator.hpp"
#include "CellsIntegrator.hpp"
#include "ExplicitIntegrator.hpp"
#include "SpatialIntegrator.hpp"
#include "TransportElement.hpp"

template<TransportElement E>
std::vector<Uint> rk4_step(RK4Integrator& integrator, Interpol& intp, ParticleSet& ps, const double t, const double dt);

// In cells on a mesh with regions; RK4's loop on an analytic field
template<TransportElement E>
std::vector<Uint> cells_step(RK4CellsIntegrator& integrator, Interpol& intp, ParticleSet& ps, const double t, const double dt);

// Whether RK4cells steps in intp's cells, a mesh with regions; fails unless
// that or an analytic field
bool steps_in_cells(Interpol& intp);

template<TransportElement E>
std::vector<Uint> explicit_step(ExplicitIntegrator& integrator, Interpol& intp, ParticleSet& ps, const double t, const double dt);

// Fields frozen at t, path length ds
template<TransportElement E>
std::vector<Uint> spatial_step(SpatialIntegrator& integrator, Interpol& intp, ParticleSet& ps, const double t, const double ds);

// Fields at the particles
void update_fields(ParticleSet& ps, Interpol& intp, const double t, const OutputFields& output_fields);

#endif
