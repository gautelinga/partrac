#ifndef __PHASE_TIMING_HPP
#define __PHASE_TIMING_HPP

// Wall-clock phases of a loader, one stderr line each, when PARTRAC_TIMING=1.
// A loader is built once from one thread, so one clock for the process is
// enough and no phase nests inside another.

namespace partrac {

// Starts a run of phases
void phase_begin(const char* what);
// Seconds since the previous phase
void phase(const char* name);
// Seconds since phase_begin
void phase_total();

}

#endif
