#include <iostream>

#include "TracerApp.hpp"

#include "tracervectors_spatial_schema.hpp"

int main(int argc, char* argv[])
{
  std::cout << "Initialized tracervectors_spatial." << std::endl;

  if (argc < 2){
    std::cout << "Please specify an input file." << std::endl;
    return 1;
  }
  partrac::Params prm = partrac::parse_or_exit(tracervectors_spatial_schema(), argc, argv);
  return run_spatial_tracers<TransportElement::Vector>(prm, "SpatialTracerVectors");
}
