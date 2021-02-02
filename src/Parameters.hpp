#ifndef __PARAMETERS_HPP
#define __PARAMETERS_HPP

#include <string>
#include <variant>
#include <map>
#include <iostream>
#include "typedefs.hpp"

class Parameters {
public:
  Parameters() {}
  Parameters(int argc, char* argv[]);
  void set_traj_defaults();
  void parse_cmd(int argc, char* argv[]);
  void parse_file(std::string infile);

  template<typename T>
  void set(std::string key, T val) { dict_[key] = val; }

  template<typename T>
  T get(std::string key) const { return std::get<T>(dict_.at(key)); }
  std::string getstr(std::string key) const { return std::get<std::string>(dict_.at(key)); }
  double getd(std::string key) const { return std::get<double>(dict_.at(key)); }
  int geti(std::string key) const { return std::get<int>(dict_.at(key)); }
  bool getb(std::string key) const { return std::get<bool>(dict_.at(key)); }
  Uint getui(std::string key) const { return std::get<Uint>(dict_.at(key)); }

  void print(std::ostream& out = std::cout) const;
  void dump(std::string folder) const;
  void dump(std::string folder, const double t) const;

  bool contains(std::string k) const;

  static void trim(std::string& s);

private:

  typedef std::variant<std::string, double, bool, int, Uint> V;
  std::map<std::string, V> dict_;

  void write_params_to_file(std::string filename) const;

  static void ltrim(std::string& s);
  static void rtrim(std::string& s);


};

#endif
