#ifndef __TESTS_TAYLOR_HOOD_HPP
#define __TESTS_TAYLOR_HOOD_HPP
#ifdef USE_DOLFIN

// The P1/P2 velocity-pressure pair the fixtures here are written in, and the
// dofs per cell of each element. The loaders read the element off the file, so
// this lives with the tests that write the files.

#include <dolfin.h>
#include <memory>
#include <string>

#include "Error.hpp"
#include "dolfin_ref.hpp"
#include "typedefs.hpp"

template<typename Cell>
void taylor_hood_spaces(const std::string& u_el, const std::string& p_el,
                        const bool include_pressure,
                        std::shared_ptr<dolfin::Mesh> mesh,
                        std::shared_ptr<const dolfin::SubDomain> constrained_domain,
                        std::shared_ptr<dolfin::FunctionSpace>& u_space,
                        std::shared_ptr<dolfin::FunctionSpace>& p_space,
                        Uint& ncoeffs_u, Uint& ncoeffs_p){
  constexpr int D = Cell::n_verts - 1;
  const auto ncoeffs_of = [](const std::string& el, const char* what) -> Uint {
    if (el == "P1") return Uint(Cell::n_verts);
    if (el == "P2") return Uint(Cell::n_dofs_max);
    partrac::fail("unrecognized ", what, " element: ", el);
  };
  ncoeffs_u = ncoeffs_of(u_el, "velocity");
  u_space = lagrange_space<D, true>(u_el, mesh, constrained_domain, "velocity");
  if (include_pressure){
    ncoeffs_p = ncoeffs_of(p_el, "pressure");
    p_space = lagrange_space<D, false>(p_el, mesh, constrained_domain, "pressure");
  }
}

#endif
#endif
