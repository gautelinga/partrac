#ifndef __INTERPOL_FACTORY_HPP
#define __INTERPOL_FACTORY_HPP

#include <memory>
#include <string>
#include "Interpol.hpp"

// Interpolator for mode, reading infilename
void set_interpolate_mode(std::shared_ptr<Interpol>& intp, const std::string& mode, const std::string& infilename);

#endif
