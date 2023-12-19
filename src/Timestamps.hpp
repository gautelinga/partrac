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
  void initialize(std::vector<std::pair<double, std::string>>& );
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

// Special classes to handle timeseries which can combine different fields and fields at different times

class MultiStamp {
public:
  MultiStamp(const double t_in, const Uint it_in) : t(t_in), it(it_in) {};
  ~MultiStamp() {};
  double t;
  Uint it;
};

class MultiStampPair {
public:
  MultiStampPair(const double t_prev, const Uint it_prev,
	               const double t_next, const Uint it_next) :
    prev(t_prev, it_prev), next(t_next, it_next) {};
  MultiStamp prev;
  MultiStamp next;
  double weight_next(const double t){ return (t-this->prev.t)/(this->next.t-this->prev.t); };
  double weight_prev(const double t){ return 1.-this->weight_next(t); };
};

class MultiTimestamps {
public:
  MultiTimestamps() { }
  void initialize(const std::vector<std::pair<double, std::vector<std::string>>>& );
  void add(const std::string&, const std::vector<std::pair<double, std::vector<std::string>>>& );
  // void update(const double);
  // std::string get_folder() const { return folder; };
  // std::string get_first() { return folder + "/" + stamps[0]; };
  // std::string get_is_solid(){ return folder + "/output_is_solid.h5"; };
  MultiStampPair get(const double);
  double get_t_min() const { return t_min; };
  double get_t_max() const { return t_max; };
  std::vector<std::string> get_path(const std::string& field, const Uint it){ return stamps[field][it]; }
private:
  std::vector<double> t_;
  std::map<std::string, std::vector<std::vector<std::string>>> stamps;
  double t_min = 1e14;
  double t_max = -1e14;
  std::string folder;
  //std::string filename;
};


#endif