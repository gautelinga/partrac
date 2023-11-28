#include <iostream>
#include "Timestamps.hpp"
#include "io.hpp"


Timestamps::Timestamps(const std::string& infilename){
  initialize(infilename);
}

void Timestamps::initialize(const std::string& infilename){
  verify_file_exists(infilename);
  std::ifstream input(infilename);
  std::string fname;
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
  std::cout << "Timestamps::update not implemented" <<std::endl;
  exit(0);
}

void Timestamps::initialize(std::vector<std::pair<double, std::string>>& items){
  for ( auto & item : items ){
    double tkey = item.first;
    std::string val = item.second;
    stamps[tkey] = val;
    if (tkey < t_min){
      t_min = tkey;
    }
    else if (tkey > t_max){
      t_max = tkey;
    }
  }
  folder = "";
  filename = "";
}

StampPair Timestamps::get(const double t){
  double t_prev, t_next;
  std::string filename_prev, filename_next;
  filename_next = stamps.begin()->second;
  t_next = stamps.begin()->first;
  for (std::map<double, std::string>::iterator it=stamps.begin()++;
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
  std::string filename_last = (--stamps.end())->second;
  if (t >= t_last){
    StampPair Pair(t_last, filename_last, t_last, filename_last);
    return Pair;
  }
  double t_first = stamps.begin()->first;
  std::string filename_first = stamps.begin()->second;
  StampPair Pair(t_first, filename_first, t_first, filename_first);
  return Pair;
}
