#ifndef __SCHEMA_CHECKS_HPP
#define __SCHEMA_CHECKS_HPP

// Shared assertions for the app schemas.

#include <catch2/catch.hpp>
#include <string>
#include "Params.hpp"

inline void check_schema(partrac::Schema& s){
  INFO(s.app());
  REQUIRE(!s.app().empty());
  REQUIRE_NOTHROW(s.validate_self());

  const std::string h = s.help();
  REQUIRE(h.find("Usage: " + s.app()) != std::string::npos);
  REQUIRE(h.find("Required:") != std::string::npos);

  // dump_intv and stat_intv must not round down to zero steps
  const std::vector<std::string> args{s.app(), "mesh.dat", "dt=0.4",
                                      "dump_intv=0.1", "stat_intv=0.1"};
  try {
    partrac::Params p = s.parse(args);
    if (p.has("dump_intv")) REQUIRE(p.get<double>("dump_intv") >= 0.4);
    if (p.has("stat_intv")) REQUIRE(p.get<double>("stat_intv") >= 0.4);
  } catch (const partrac::ParamError&) {
    // needs more required parameters than this probe supplies
  }
}

#endif
