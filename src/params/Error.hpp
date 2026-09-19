#ifndef __PARTRAC_ERROR_HPP
#define __PARTRAC_ERROR_HPP

// A run that cannot go on: bad input, a broken file, a broken invariant
// (its message then begins "internal:"). The libraries throw it; each main
// reports it through report_errors.

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace partrac {

class Error : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

// Throws Error, the arguments streamed into its message
template<typename... Args>
[[noreturn]] void fail(const Args&... args){
  std::ostringstream s;
  (s << ... << args);
  throw Error(s.str());
}

// A main's body; an Error goes to stderr and the run exits with 2, as a parameter
// error; any other exception, from dolfin or the standard library, with 1
template<typename Body>
int report_errors(Body&& body){
  try {
    return body();
  } catch (const Error& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 2;
  } catch (const std::exception& e) {
    std::cerr << "Error (unexpected): " << e.what() << std::endl;
    return 1;
  }
}

}  // namespace partrac

#endif
