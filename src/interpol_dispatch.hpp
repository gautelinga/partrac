#ifndef __INTERPOL_DISPATCH_HPP
#define __INTERPOL_DISPATCH_HPP

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

// Call f with the interpolator as its concrete type, so that the particle loop
// inside f is compiled against locate and evaluate directly: no virtual call
// per particle per stage, and the loop-invariant loads hoisted out of it. The
// interpolator is still chosen at run time from the mode string, once here
// rather than once per particle. Anything not listed falls back to the base
// class and the virtual calls it always had.
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
  return f(ip);
}

#endif
