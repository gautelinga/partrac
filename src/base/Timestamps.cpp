#include <iostream>
#include <iterator>
#include "Error.hpp"
#include "Timestamps.hpp"
#include "files.hpp"


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
    // Unsorted input: separate tests
    if (key < t_min){
      t_min = key;
    }
    if (key > t_max){
      t_max = key;
    }
  }
  std::size_t botDirPos = infilename.find_last_of("/");
  std::size_t extPos = infilename.find_last_of(".");

  folder = infilename.substr(0, botDirPos);
  filename = infilename.substr(botDirPos+1, extPos-botDirPos-1);
}

void Timestamps::update(const double){
  partrac::fail("Timestamps::update is not implemented");
}

void Timestamps::initialize(std::vector<std::pair<double, std::string>>& items){
  for ( auto & item : items ){
    double tkey = item.first;
    std::string val = item.second;
    stamps[tkey] = val;
    if (tkey < t_min){
      t_min = tkey;
    }
    if (tkey > t_max){
      t_max = tkey;
    }
  }
  folder = "";
  filename = "";
}

StampPair Timestamps::get(const double t){
  if (stamps.empty())
    partrac::fail("Timestamps: no time stamps");
  // First stamp after t
  const auto next = stamps.upper_bound(t);
  // Before the first or past the last: that stamp twice
  if (next == stamps.begin())
    return StampPair(next->first, next->second, next->first, next->second);
  const auto prev = std::prev(next);
  if (next == stamps.end())
    return StampPair(prev->first, prev->second, prev->first, prev->second);
  return StampPair(prev->first, prev->second, next->first, next->second);
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
    if (tkey > t_max){
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
      partrac::fail("XDMF: the time keys of field ", field, " do not match");
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