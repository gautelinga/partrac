#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
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

double Timestamps::next_after(const double t) const {
  const auto next = stamps.upper_bound(t);
  return next == stamps.end() ? std::numeric_limits<double>::infinity() : next->first;
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
  // Searched by bisection
  if (!std::is_sorted(t_.begin(), t_.end()))
    partrac::fail("XDMF: the time keys are not in increasing order");
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
  if (t_.empty())
    partrac::fail("XDMF: no time stamps");
  // First stamp after t
  const Uint next = std::upper_bound(t_.begin(), t_.end(), t) - t_.begin();
  // Before the first or past the last: that stamp twice
  if (next == 0)
    return MultiStampPair(t_[0], 0, t_[0], 0);
  if (next == t_.size())
    return MultiStampPair(t_[next-1], next-1, t_[next-1], next-1);
  return MultiStampPair(t_[next-1], next-1, t_[next], next);
}

double MultiTimestamps::next_after(const double t) const {
  const auto next = std::upper_bound(t_.begin(), t_.end(), t);
  return next == t_.end() ? std::numeric_limits<double>::infinity() : *next;
}
