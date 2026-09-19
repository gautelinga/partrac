#ifdef USE_DOLFIN
#include "Error.hpp"
#include "dolfin_spaces.hpp"
#include "Triangle.hpp"
#include "Tet.hpp"
#include "dolfin_elements/P1_2.h"
#include "dolfin_elements/P2_2.h"
#include "dolfin_elements/vP1_2.h"
#include "dolfin_elements/vP2_2.h"
#include "dolfin_elements/P1_3.h"
#include "dolfin_elements/P2_3.h"
#include "dolfin_elements/vP1_3.h"
#include "dolfin_elements/vP2_3.h"
#include <iostream>

namespace {

// P1 or P2 of this dimension; ncoeffs is the element's dofs per cell
template<typename Cell, bool Vector>
std::shared_ptr<dolfin::FunctionSpace> space_of(const std::string& el,
                                                std::shared_ptr<dolfin::Mesh> mesh,
                                                std::shared_ptr<const dolfin::SubDomain> constrained_domain,
                                                Uint& ncoeffs, const char* what){
  constexpr bool two_d = Cell::n_verts == 3;
  if (el == "P1"){
    ncoeffs = Cell::n_verts;
    if constexpr (two_d && Vector)   return std::make_shared<vP1_2::FunctionSpace>(mesh, constrained_domain);
    else if constexpr (two_d)        return std::make_shared<P1_2::FunctionSpace>(mesh, constrained_domain);
    else if constexpr (Vector)       return std::make_shared<vP1_3::FunctionSpace>(mesh, constrained_domain);
    else                             return std::make_shared<P1_3::FunctionSpace>(mesh, constrained_domain);
  }
  if (el == "P2"){
    ncoeffs = Cell::n_dofs_max;
    if constexpr (two_d && Vector)   return std::make_shared<vP2_2::FunctionSpace>(mesh, constrained_domain);
    else if constexpr (two_d)        return std::make_shared<P2_2::FunctionSpace>(mesh, constrained_domain);
    else if constexpr (Vector)       return std::make_shared<vP2_3::FunctionSpace>(mesh, constrained_domain);
    else                             return std::make_shared<P2_3::FunctionSpace>(mesh, constrained_domain);
  }
  partrac::fail("unrecognized ", what, " element: ", el);
}

}  // namespace

template<typename Cell>
void taylor_hood_spaces(const std::string& u_el, const std::string& p_el,
                        const bool include_pressure,
                        std::shared_ptr<dolfin::Mesh> mesh,
                        std::shared_ptr<const dolfin::SubDomain> constrained_domain,
                        std::shared_ptr<dolfin::FunctionSpace>& u_space,
                        std::shared_ptr<dolfin::FunctionSpace>& p_space,
                        Uint& ncoeffs_u, Uint& ncoeffs_p){
  // Velocity
  u_space = space_of<Cell, true>(u_el, mesh, constrained_domain, ncoeffs_u, "velocity");

  // Pressure
  if (include_pressure){
    p_space = space_of<Cell, false>(p_el, mesh, constrained_domain, ncoeffs_p, "pressure");
  }
  else {
    std::cout << "Note: Ignoring pressure." << std::endl;
  }
}

template void taylor_hood_spaces<Triangle>(const std::string&, const std::string&, const bool,
                                           std::shared_ptr<dolfin::Mesh>,
                                           std::shared_ptr<const dolfin::SubDomain>,
                                           std::shared_ptr<dolfin::FunctionSpace>&,
                                           std::shared_ptr<dolfin::FunctionSpace>&, Uint&, Uint&);
template void taylor_hood_spaces<Tet>(const std::string&, const std::string&, const bool,
                                      std::shared_ptr<dolfin::Mesh>,
                                      std::shared_ptr<const dolfin::SubDomain>,
                                      std::shared_ptr<dolfin::FunctionSpace>&,
                                      std::shared_ptr<dolfin::FunctionSpace>&, Uint&, Uint&);

#endif
