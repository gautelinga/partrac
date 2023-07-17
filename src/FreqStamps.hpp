#ifndef __FREQSTAMPS_HPP
#define __FREQSTAMPS_HPP

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include "utils.hpp"
#include "io.hpp"

class FreqStamp {
public:
  FreqStamp( const double t_in
           , const double a_in
           , const std::string filename_in
           ) : t(t_in), a(a_in), filename(filename_in) {};
  ~FreqStamp() {};
  double t; // time shift
  double a; // amplitude
  std::string filename;
};

class FreqStamps {
public:
  FreqStamps() { }
  //FreqStamps(const std::string&);
  void initialize(const std::string& infilename){
    verify_file_exists(infilename);
    std::ifstream input(infilename);
    std::string fname;
    double t, a;
    while (input >> t >> a >> fname){
      FreqStamp f(t, a, fname);
      ordered_stamps.push_back(f);
    }
    std::size_t botDirPos = infilename.find_last_of("/");
    std::size_t extPos = infilename.find_last_of(".");

    folder = infilename.substr(0, botDirPos);
    filename = infilename.substr(botDirPos+1, extPos-botDirPos-1);
  };
  std::string get_folder() const { return folder; };
  FreqStamp& get(const int i) { return ordered_stamps[i]; };
  int size() const { return ordered_stamps.size(); };
  std::string get_filename(const int i) { return ordered_stamps[i].filename; };
  double get_t(const int i) { return ordered_stamps[i].t; };
  double get_a(const int i) { return ordered_stamps[i].a; };
private:
  std::vector<FreqStamp> ordered_stamps;  // ordered by line number
  std::string folder;
  std::string filename;
};

#endif