#ifndef __STAMP_BUFFER_HPP
#define __STAMP_BUFFER_HPP

// The two stamps a loader blends a field between: one buffer per stamp, keyed
// by the stamp it holds, so a stamp a buffer already holds is never read again.

#include <string>
#include <utility>

namespace partrac {

// The two times of a bracket
struct StampTimes {
  double prev = 0.;
  double next = 0.;
};

// Always load once; keep the last bracket past t_max
inline bool stamp_reload(const bool held, const double t, const double t_max,
                         const StampTimes& have, const StampTimes& want){
  return !held || ((have.prev != want.prev || have.next != want.next) && t < t_max);
}

// What holding a wanted stamp took
enum class StampFill { Held, Swapped, Read, Aliased };

// What an update reports for one stamp
inline std::string stamp_note(const StampFill fill, const std::string& name){
  switch (fill){
  case StampFill::Swapped: return "swapping...";
  case StampFill::Read:    return "file = " + name;
  case StampFill::Aliased: return "the previous stamp";
  default:                 return "held";
  }
}

// Stamp carries one stamp's values and whatever was computed from them; Key is
// what a loader calls a stamp
template<typename Stamp, typename Key>
class StampBuffer {
public:
  // Buffer a, which a loader's setup may fill with the first stamp itself
  Stamp& a() { return a_; }
  void hold_a(const Key& key){ key_a_ = key; a_held_ = true; }
  const Key& key_a() const { return key_a_; }

  // The wanted bracket, whatever is missing read through read(key, stamp)
  template<typename Read>
  std::pair<StampFill, StampFill> load(const Key& prev_key, const Key& next_key, Read&& read);

  const Stamp& prev() const { return *prev_; }
  const Stamp& next() const { return *next_; }
  // The two stamps are the same values, not a copy of them
  bool aliased() const { return prev_ == next_; }

private:
  static bool holds(const bool held, const Key& key, const Key& want){
    return held && key == want;
  }

  Stamp a_, b_;
  Key key_a_{}, key_b_{};
  bool a_held_ = false, b_held_ = false;
  const Stamp* prev_ = nullptr;
  const Stamp* next_ = nullptr;
};

template<typename Stamp, typename Key>
template<typename Read>
std::pair<StampFill, StampFill> StampBuffer<Stamp, Key>::load(const Key& prev_key,
                                                              const Key& next_key, Read&& read)
{
  StampFill prev_fill = StampFill::Held;
  // The next stamp's buffer may hold the one wanted as previous, or the
  // previous stamp's the one wanted as next
  const bool prev_in_b = holds(b_held_, key_b_, prev_key) && !holds(a_held_, key_a_, prev_key);
  const bool next_in_a = next_key != prev_key && holds(a_held_, key_a_, next_key)
                      && !holds(b_held_, key_b_, next_key);
  if (prev_in_b || next_in_a){
    std::swap(a_, b_);
    std::swap(key_a_, key_b_);
    std::swap(a_held_, b_held_);
    if (prev_in_b) prev_fill = StampFill::Swapped;
  }
  if (!holds(a_held_, key_a_, prev_key)){
    read(prev_key, a_);
    key_a_ = prev_key;
    a_held_ = true;
    prev_fill = StampFill::Read;
  }
  prev_ = &a_;

  // Single stamp: the same values, not a copy of them
  if (next_key == prev_key){
    next_ = &a_;
    return {prev_fill, StampFill::Aliased};
  }

  StampFill next_fill = StampFill::Held;
  if (!holds(b_held_, key_b_, next_key)){
    read(next_key, b_);
    key_b_ = next_key;
    b_held_ = true;
    next_fill = StampFill::Read;
  }
  next_ = &b_;
  return {prev_fill, next_fill};
}

}  // namespace partrac

#endif
