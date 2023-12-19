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
  // TODO: Make more efficient
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

//
void MultiTimestamps::initialize(const std::vector<std::pair<double, std::vector<std::string>>>& items){
  t_.resize(items.size());
  stamps["u"].resize(items.size()); // initialize vector
  for (Uint i=0; i < items.size(); ++i){
    auto & item = items[i];
    double tkey = item.first;
    t_[i] = tkey;
    stamps["u"][i] = item.second;
    if (tkey < t_min){
      t_min = tkey;
    }
    else if (tkey > t_max){
      t_max = tkey;
    }
  }
  //folder = "";
  //filename = "";
}

void MultiTimestamps::add(const std::string& field, const std::vector<std::pair<double, std::vector<std::string>>>& items){
  stamps[field].resize(items.size()); // initialize vector
  for (Uint i=0; i < items.size(); ++i){
    auto tkey = items[i].first;
    if (abs(t_[i] - tkey) > 1e-10)
    {
      std::cout << "ERROR: XDMF tkey for field " << field << " is not matching!" << std::endl;
      exit(0);
    }
    stamps[field][i] = items[i].second;
  }
}

MultiStampPair MultiTimestamps::get(const double t){
  if (t_.size() > 0)
  {
    for (Uint _it=1; _it < t_.size(); ++_it)
    {
      if (t_[_it-1] <= t && t_[_it] > t)
      {
        MultiStampPair Pair(t_[_it-1], _it-1, t_[_it], _it);
        return Pair;
      }
    }
  }
  Uint it_last = t_.size()-1;
  double t_last = t_[it_last];
  if (t >= t_last){
    MultiStampPair Pair(t_last, it_last, t_last, it_last);
    return Pair;
  }
  double t_first = t_[0];
  MultiStampPair Pair(t_first, 0, t_first, 0);
  return Pair;
}