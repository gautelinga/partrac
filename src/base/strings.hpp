#ifndef __STRINGS_HPP
#define __STRINGS_HPP

#include <map>
#include <set>
#include <string>
#include <vector>

template<typename T>
bool contains(const std::set<T>& container, const T &elem){
  return container.find(elem) != container.end();
}

template<typename T1, typename T2>
bool contains(const std::map<T1, T2>& container, const T1 &elem){
  return container.find(elem) != container.end();
}

inline std::vector<std::string> split_string(const std::string s, const std::string delim){
  std::vector<std::string> s_;
  auto start = 0U;
  auto end = s.find(delim);
  while (end != std::string::npos){
    s_.push_back(s.substr(start, end-start));
    start = end + delim.length();
    end = s.find(delim, start);
  }
  s_.push_back(s.substr(start, end));
  return s_;
}

inline bool contains(const std::string s, const std::string c){
  return (s.find(c) != std::string::npos);
}

// Directions after the shape are non-empty combinations of x, y and z.
inline bool init_mode_dirs_ok(const std::vector<std::string>& key){
  for (std::size_t i = 1; i < key.size(); ++i)
    if (key[i].empty() || key[i].find_first_not_of("xyz") != std::string::npos)
      return false;
  return true;
}

#endif
