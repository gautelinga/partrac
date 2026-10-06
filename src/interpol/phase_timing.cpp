#include "phase_timing.hpp"

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iostream>

namespace partrac {

namespace {
using phase_clock = std::chrono::steady_clock;
phase_clock::time_point start_, last_;
const char* what_ = "";

bool timing_on(){
  static const bool on = [](){
    const char* e = std::getenv("PARTRAC_TIMING");
    return e && std::strcmp(e, "1") == 0;
  }();
  return on;
}
}

void phase_begin(const char* what){
  if (!timing_on()) return;
  what_ = what;
  start_ = last_ = phase_clock::now();
}

void phase(const char* name){
  if (!timing_on()) return;
  const auto now = phase_clock::now();
  const double dt = std::chrono::duration<double>(now - last_).count();
  last_ = now;
  std::cerr << "[timing] " << what_ << ": " << name << " " << dt << " s" << std::endl;
}

void phase_total(){
  if (!timing_on()) return;
  const double dt = std::chrono::duration<double>(phase_clock::now() - start_).count();
  std::cerr << "[timing] " << what_ << ": total " << dt << " s" << std::endl;
}

}
