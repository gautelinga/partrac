#include <iostream>

#include "Error.hpp"
#include "wait_policy.hpp"
#include "TracerApp.hpp"

#include "tracervectors_spatial_schema.hpp"

static int run(int argc, char* argv[])
{
  std::cout << "Initialized tracervectors_spatial." << std::endl;

  partrac::Params prm = partrac::parse_or_exit(tracervectors_spatial_schema(), argc, argv);
  return run_spatial_tracers<TransportElement::Vector>(prm, "SpatialTracerVectors");
}

int main(int argc, char* argv[])
{
  partrac::wait_passively(argv);
  return partrac::report_errors([&]{ return run(argc, argv); });
}
