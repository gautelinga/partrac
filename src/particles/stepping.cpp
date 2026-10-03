#include <type_traits>
#include "stepping.hpp"
#include "held_eval.hpp"
#include "interpol_dispatch.hpp"

template<TransportElement E>
std::vector<Uint> rk4_step(RK4Integrator& integrator, Interpol& intp, ParticleSet& ps, const double t, const double dt){
  return with_concrete(intp, [&](auto& ip){ return integrator.template step<E>(ip, ps, t, dt); });
}

template<TransportElement E>
std::vector<Uint> cells_step(RK4CellsIntegrator& integrator, Interpol& intp, ParticleSet& ps, const double t, const double dt){
  return with_concrete(intp, [&](auto& ip){
    using I = std::decay_t<decltype(ip)>;
    if constexpr (has_regions<I>::value)
      return integrator.template step<E>(ip, ps, t, dt);
    else
      return integrator.RK4Integrator::template step<E>(ip, ps, t, dt);
  });
}

// The mesh loaders' cells; 0 for the rest
template<typename Interp, typename = void>
struct has_cell_count : std::false_type {};
template<typename Interp>
struct has_cell_count<Interp, std::void_t<decltype(std::declval<const Interp&>().cell_count())>> : std::true_type {};

Uint cell_count(Interpol& intp){
  return with_concrete(intp, [](auto& ip) -> Uint {
    if constexpr (has_cell_count<std::decay_t<decltype(ip)>>::value) return ip.cell_count();
    else return 0;
  });
}

bool steps_in_cells(Interpol& intp){
  return with_concrete(intp, [](auto& ip){
    using I = std::decay_t<decltype(ip)>;
    if (!has_regions<I>::value && !std::is_same_v<I, AnalyticInterpol>)
      partrac::fail("scheme=RK4cells steps in the cells of a mesh (not mode=fenics) or of felbm "
                    "(interpolation=linear), or on an analytic field; use RK4 for this one");
    return has_regions<I>::value;
  });
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
template std::vector<Uint> cells_step<TransportElement::Point>(RK4CellsIntegrator&, Interpol&, ParticleSet&, const double, const double);
template std::vector<Uint> cells_step<TransportElement::Vector>(RK4CellsIntegrator&, Interpol&, ParticleSet&, const double, const double);
template std::vector<Uint> cells_step<TransportElement::Tensor>(RK4CellsIntegrator&, Interpol&, ParticleSet&, const double, const double);
template std::vector<Uint> explicit_step<TransportElement::Point>(ExplicitIntegrator&, Interpol&, ParticleSet&, const double, const double);
template std::vector<Uint> explicit_step<TransportElement::Vector>(ExplicitIntegrator&, Interpol&, ParticleSet&, const double, const double);
template std::vector<Uint> explicit_step<TransportElement::Tensor>(ExplicitIntegrator&, Interpol&, ParticleSet&, const double, const double);
template std::vector<Uint> spatial_step<TransportElement::Point>(SpatialIntegrator&, Interpol&, ParticleSet&, const double, const double);
template std::vector<Uint> spatial_step<TransportElement::Vector>(SpatialIntegrator&, Interpol&, ParticleSet&, const double, const double);
