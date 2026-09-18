#ifdef USE_DOLFIN
#ifndef __DOLFIN_SPACES_HPP
#define __DOLFIN_SPACES_HPP

// The Taylor-Hood element spaces of a mesh loader, by the names in its
// parameter file; the element classes differ by dimension

#include <memory>
#include <string>
#include <dolfin.h>
#include "typedefs.hpp"

template<typename Cell>
void taylor_hood_spaces(const std::string& u_el, const std::string& p_el,
                        const bool include_pressure,
                        std::shared_ptr<dolfin::Mesh> mesh,
                        std::shared_ptr<const dolfin::SubDomain> constrained_domain,
                        std::shared_ptr<dolfin::FunctionSpace>& u_space,
                        std::shared_ptr<dolfin::FunctionSpace>& p_space,
                        Uint& ncoeffs_u, Uint& ncoeffs_p);

#endif
#endif
