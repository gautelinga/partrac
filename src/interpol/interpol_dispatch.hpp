#ifndef __INTERPOL_DISPATCH_HPP
#define __INTERPOL_DISPATCH_HPP

#include <iostream>
#include <typeinfo>
#include "Error.hpp"
#include "Interpol.hpp"
#include "interpolators.hpp"

// Call f with the concrete interpolator, tried in the table's order (src/CMakeLists.txt)
template<typename F>
inline auto with_concrete(Interpol& ip, F&& f){
#define PARTRAC_TRY(T) if (auto* p = dynamic_cast<T*>(&ip)) return f(*p);
  PARTRAC_INTERPOLATORS(PARTRAC_TRY)
#undef PARTRAC_TRY
  partrac::fail("with_concrete: interpolator type ", typeid(ip).name(), " is not in the dispatch list");
  return decltype(f(ip))();   // for the return type only; not reached
}

#endif
