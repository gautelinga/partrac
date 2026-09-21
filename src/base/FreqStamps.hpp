#ifndef __FREQSTAMPS_HPP
#define __FREQSTAMPS_HPP

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include "Error.hpp"
#include "files.hpp"

class FreqStamp {
public:
  FreqStamp( const double t_in
           , const double omega_in
           , const double a_in
           , const std::string filename_in
           ) : t(t_in), omega(omega_in), a(a_in), filename(filename_in) {};
  ~FreqStamp() {};
  double t;       // `t a file`: phase in time, cos(omega0*(k*t_run + t)); `omega phi a file`: phi
  double omega;   // `omega phi a file`: angular frequency, cos(omega*t_run + phi)
  double a;       // amplitude
  std::string filename;
};

class FreqStamps {
public:
  FreqStamps() { }
  //FreqStamps(const std::string&);
  // A line is `t a file`, a harmonic of the base frequency with k its line
  // number, or `omega phi a file`; one form a file. A frequency may appear on
  // several lines, as the cosine and sine parts of a Fourier mode (phi = 0 and
  // -pi/2).
  void initialize(const std::string& infilename){
    verify_file_exists(infilename);
    std::ifstream input(infilename);
    std::string line;
    std::size_t lineno = 0, width = 0;
    while (std::getline(input, line)){
      ++lineno;
      if (!line.empty() && line.back() == '\r') line.pop_back();
      std::istringstream ls(line);
      std::vector<std::string> tok;
      for (std::string w; ls >> w; ) tok.push_back(w);
      // Blank lines and comments
      if (tok.empty() || tok[0][0] == '#') continue;
      if (tok.size() != 3 && tok.size() != 4)
        partrac::fail(infilename, ":", lineno, ": ", tok.size(), " columns; a line is `t a file` or `omega phi a file`");
      if (width != 0 && tok.size() != width)
        partrac::fail(infilename, ":", lineno, ": ", tok.size(), " columns after lines of ", width, "; one form a file");
      width = tok.size();
      std::vector<double> v(width - 1);
      for (std::size_t i = 0; i + 1 < width; ++i){
        std::size_t used = 0;
        try { v[i] = std::stod(tok[i], &used); }
        catch (const std::exception&){ used = 0; }
        if (used == 0 || used != tok[i].size())
          partrac::fail(infilename, ":", lineno, ": '", tok[i], "' is not a number");
      }
      const int k = int(ordered_stamps.size());
      const bool mean = width == 3 ? k == 0 : v[0] == 0.;
      const double phase = width == 3 ? v[0] : v[1];
      // A mean has no phase: its cosine would only rescale it
      if (mean && phase != 0.)
        partrac::fail(infilename, ":", lineno, ": a component of frequency 0 has phase ", phase,
                      "; it must be 0");
      if (width == 3) ordered_stamps.push_back(FreqStamp(v[0], 0., v[1], tok[2]));
      else            ordered_stamps.push_back(FreqStamp(v[1], v[0], v[2], tok[3]));
    }
    omega_given_ = width == 4;
  };
  FreqStamp& get(const int i) { return ordered_stamps[i]; };
  int size() const { return ordered_stamps.size(); };
  // Lines of `omega phi a file`, which need no base frequency
  bool omega_given() const { return omega_given_; };
private:
  std::vector<FreqStamp> ordered_stamps;  // in the file's order
  bool omega_given_ = false;
};

#endif