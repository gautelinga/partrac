#ifndef __INTERPOL_DISPATCH_HPP
#define __INTERPOL_DISPATCH_HPP

#include <iostream>
#include <typeinfo>
#include "Error.hpp"
#include "Interpol.hpp"
#include "AnalyticInterpol.hpp"
#include "StructuredInterpol.hpp"
#include "TriangleInterpol.hpp"
#include "TetInterpol.hpp"
#include "TriangleFreqInterpol.hpp"
#include "XDMFTriangleInterpol.hpp"
#include "XDMFTetInterpol.hpp"
#ifdef USE_DOLFIN
#include "DolfTriangleInterpol.hpp"
#include "DolfTetInterpol.hpp"
#endif

// Call f with the concrete interpolator; the same list as set_interpolate_mode
// and PARTRAC_STEP_INTERPOLATORS in src/CMakeLists.txt
template<typename F>
inline auto with_concrete(Interpol& ip, F&& f){
  if (auto* p = dynamic_cast<AnalyticInterpol*>(&ip)) return f(*p);
  if (auto* p = dynamic_cast<StructuredInterpol*>(&ip)) return f(*p);
  if (auto* p = dynamic_cast<StructuredConstInterpol*>(&ip)) return f(*p);
  if (auto* p = dynamic_cast<TriangleInterpol*>(&ip)) return f(*p);
  if (auto* p = dynamic_cast<TetInterpol*>(&ip)) return f(*p);
  if (auto* p = dynamic_cast<TriangleFreqInterpol*>(&ip)) return f(*p);
  if (auto* p = dynamic_cast<XDMFTriangleInterpol*>(&ip)) return f(*p);
  if (auto* p = dynamic_cast<XDMFTetInterpol*>(&ip)) return f(*p);
#ifdef USE_DOLFIN
  if (auto* p = dynamic_cast<DolfTriangleInterpol*>(&ip)) return f(*p);
  if (auto* p = dynamic_cast<DolfTetInterpol*>(&ip)) return f(*p);
#endif
  partrac::fail("with_concrete: interpolator type ", typeid(ip).name(), " is not in the dispatch list");
  return decltype(f(ip))();   // for the return type only; not reached
}

#endif
