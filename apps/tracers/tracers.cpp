#include <iostream>

#include "TracerApp.hpp"

#include "tracers_schema.hpp"

int main(int argc, char* argv[])
{
  std::cout << "Initialized tracers." << std::endl;

  if (argc < 2){
    std::cout << "Please specify an input file." << std::endl;
    return 1;
  }
  partrac::Params prm = partrac::parse_or_exit(tracers_schema(), argc, argv);
  return run_tracers<TransportElement::Point>(prm, "Tracers");
}
