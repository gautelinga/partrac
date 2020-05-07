#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "io.hpp"
#include "utils.hpp"

#ifndef __TIMESTAMPS_HPP
#define __TIMESTAMPS_HPP

using namespace std;

class Timestamps {
public:
  Timestamps(string);
  void update(const double);
  string get_folder(){ return folder; };
  string get_first(){ return folder + "/" + stamps[0]; };
  string get_is_solid(){ return folder + "/output_is_solid.h5"; };
  StampPair get(const double);
  double get_t_min() { return t_min; };
  double get_t_max() { return t_max; };
private:
  map<double, string> stamps;
  double t_min = 1e14;
  double t_max = -1e14;
  string folder;
  string filename;
};

Timestamps::Timestamps(string infilename){
  verify_file_exists(infilename);
  ifstream input(infilename);
  string fname;
  double key;
  while (input >> key >> fname){
    stamps[key] = fname;
    if (key < t_min){
      t_min = key;
    }
    else if (key > t_max){
      t_max = key;
    }
  }
  std::size_t botDirPos = infilename.find_last_of("/");
  std::size_t extPos = infilename.find_last_of(".");

  folder = infilename.substr(0, botDirPos);
  filename = infilename.substr(botDirPos+1, extPos-botDirPos-1);
}

void Timestamps::update(const double){
  exit(0);
}

StampPair Timestamps::get(const double t){
  double t_prev, t_next;
  string filename_prev, filename_next;
  filename_next = stamps.begin()->second;
  t_next = stamps.begin()->first;
  for (map<double, string>::iterator it=stamps.begin()++;
       it!=stamps.end(); ++it){
    filename_prev = filename_next;
    filename_next = it->second;
    t_prev = t_next;
    t_next = it->first;
    if (t_prev <= t && t_next > t){
      StampPair Pair(t_prev, filename_prev,
		     t_next, filename_next);
      return Pair;
    }
  }
  double t_last = (--stamps.end())->first;
  string filename_last = (--stamps.end())->second;
  if (t >= t_last){
    StampPair Pair(t_last, filename_last, t_last, filename_last);
    return Pair;
  }
  double t_first = stamps.begin()->first;
  string filename_first = stamps.begin()->second;
  StampPair Pair(t_first, filename_first, t_first, filename_first);
  return Pair;
}

#endif
