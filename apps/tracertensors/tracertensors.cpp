#include <iostream>

#include "Error.hpp"
#include "wait_policy.hpp"
#include "TracerApp.hpp"

#include "tracertensors_schema.hpp"

static int run(int argc, char* argv[])
{
  std::cout << "Initialized tracertensors." << std::endl;

  partrac::Params prm = partrac::parse_or_exit(tracertensors_schema(), argc, argv);
  return run_tracers<TransportElement::Tensor>(prm, "TracerTensors");
}

int main(int argc, char* argv[])
{
  partrac::wait_passively(argv);
  return partrac::report_errors([&]{ return run(argc, argv); });
}
