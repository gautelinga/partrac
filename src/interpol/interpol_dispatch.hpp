#ifndef __INTERPOL_DISPATCH_HPP
#define __INTERPOL_DISPATCH_HPP

#include <iostream>
#include <typeinfo>
#include "Interpol.hpp"
#include "AnalyticInterpol.hpp"
#include "StructuredInterpol.hpp"
#ifdef USE_DOLFIN
#include "TriangleInterpol.hpp"
#include "TetInterpol.hpp"
#include "TriangleFreqInterpol.hpp"
#include "XDMFTriangleInterpol.hpp"
#include "XDMFTetInterpol.hpp"
#include "DolfInterpol.hpp"
#endif

// Call f with the concrete interpolator; the same list as set_interpolate_mode
// and PARTRAC_STEP_INTERPOLATORS in src/CMakeLists.txt
template<typename F>
inline auto with_concrete(Interpol& ip, F&& f){
  if (auto* p = dynamic_cast<AnalyticInterpol*>(&ip)) return f(*p);
  if (auto* p = dynamic_cast<StructuredInterpol*>(&ip)) return f(*p);
#ifdef USE_DOLFIN
  if (auto* p = dynamic_cast<TriangleInterpol*>(&ip)) return f(*p);
  if (auto* p = dynamic_cast<TetInterpol*>(&ip)) return f(*p);
  if (auto* p = dynamic_cast<TriangleFreqInterpol*>(&ip)) return f(*p);
  if (auto* p = dynamic_cast<XDMFTriangleInterpol*>(&ip)) return f(*p);
  if (auto* p = dynamic_cast<XDMFTetInterpol*>(&ip)) return f(*p);
  if (auto* p = dynamic_cast<DolfInterpol*>(&ip)) return f(*p);
#endif
  std::cerr << "with_concrete: interpolator type " << typeid(ip).name()
            << " is not in the dispatch list" << std::endl;
  exit(1);
  return decltype(f(ip))();   // for the return type only; not reached
}

#endif
