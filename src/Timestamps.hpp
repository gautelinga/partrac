#ifndef __TIMESTAMPS_HPP
#define __TIMESTAMPS_HPP

#include <fstream>
#include <map>
#include <string>
#include "utils.hpp"

class Timestamps {
public:
  Timestamps() { }
  Timestamps(const std::string&);
  void initialize(const std::string&);
  void update(const double);
  std::string get_folder() const { return folder; };
  std::string get_first() { return folder + "/" + stamps[0]; };
  //  std::string get_is_solid(){ return folder + "/output_is_solid.h5"; };
  StampPair get(const double);
  double get_t_min() const { return t_min; };
  double get_t_max() const { return t_max; };
private:
  std::map<double, std::string> stamps;
  double t_min = 1e14;
  double t_max = -1e14;
  std::string folder;
  std::string filename;
};

#endif