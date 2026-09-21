#ifndef __PARTRAC_PREDICATES_HPP
#define __PARTRAC_PREDICATES_HPP

// Shewchuk's adaptive exact geometric predicates (predicates.c, public domain,
// kept verbatim). exactinit() must run once before any predicate call.

extern "C" {
  void exactinit(void);
  double orient2d(double* pa, double* pb, double* pc);
  double orient3d(double* pa, double* pb, double* pc, double* pd);
}

#endif
