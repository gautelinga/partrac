#ifndef __INTERPOL_FACTORY_HPP
#define __INTERPOL_FACTORY_HPP

#include <memory>
#include <string>
#include <vector>
#include "Interpol.hpp"

// The modes set_interpolate_mode reads, for the schemas
inline std::vector<std::string> interpol_modes(){
  return {"analytic", "structured", "lbm", "felbm", "fenics",
          "tet", "triangle", "trianglefreq", "tetfreq", "xdmftriangle", "xdmftet", "openfoam"};
}

// Interpolator for mode, reading infilename
void set_interpolate_mode(std::shared_ptr<Interpol>& intp, const std::string& mode, const std::string& infilename);

#endif
