#ifndef __PARAMS_HPP
#define __PARAMS_HPP

// Each app declares the parameters it accepts in a Schema, then parses the
// command line against it. Unknown keys and missing required ones are errors.
// Parsing throws ParamError rather than calling exit(), so it can be tested.

#include <cstddef>
#include <functional>
#include <iosfwd>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

namespace partrac {

// Every parameter is stored as one of four kinds.  Narrowing to the type the
// caller asks for happens in get<T>(), which is where the range checks live.
enum class Kind { Bool, Int, Real, Str };

// Where a resolved value came from.  Used by print() and by warn() predicates
// that need to distinguish "user asked for this" from "this is the default".
enum class Source { Unset, Default, File, Cmdline, Runtime };

using Value = std::variant<bool, long long, double, std::string>;

std::string to_string(Kind k);
std::string to_string(Source s);
std::string value_to_string(const Value& v);

// Carries every problem found during one parse, not just the first.
class ParamError : public std::runtime_error {
public:
  ParamError(std::string app, std::vector<std::string> problems);
  const std::vector<std::string>& problems() const { return m_problems; }
  const std::string& app() const { return m_app; }
private:
  std::string m_app;
  std::vector<std::string> m_problems;
};

class Params;
class Schema;

using Pred = std::function<bool(const Params&)>;

namespace detail {

template <typename> inline constexpr bool always_false = false;

template <typename T>
constexpr Kind kind_of() {
  if constexpr (std::is_same_v<T, bool>)             return Kind::Bool;
  else if constexpr (std::is_same_v<T, std::string>) return Kind::Str;
  else if constexpr (std::is_floating_point_v<T>)    return Kind::Real;
  else if constexpr (std::is_integral_v<T>)          return Kind::Int;
  else static_assert(always_false<T>,
      "parameter type must be bool, an integer, a float, or std::string");
}

template <typename T>
Value to_value(T v) {
  if constexpr (std::is_same_v<T, bool>)             return Value(v);
  else if constexpr (std::is_same_v<T, std::string>) return Value(std::move(v));
  else if constexpr (std::is_floating_point_v<T>)    return Value(static_cast<double>(v));
  else                                               return Value(static_cast<long long>(v));
}

struct Entry {
  std::string key;
  Kind kind = Kind::Real;
  std::string doc;
  bool required = false;      // require<T>
  bool is_runtime = false;    // set by the program; rejected on the command line
  std::optional<Value> def;   // absent for require<T> and require_if<T>
  Pred required_when;         // require_if<T>
  std::string required_why;
  std::vector<std::string> choices;
  std::string choices_sep;      // if set, choices apply to the first token
};

struct Constraint {
  Pred pred;
  std::string message;
};

struct SchemaImpl {
  std::string app;
  std::string positional_doc;
  std::vector<Entry> entries;              // registration order, drives dump/print
  std::map<std::string, std::size_t> index;
  std::vector<Constraint> checks;
  std::vector<Constraint> warns;
  std::vector<std::function<void(Params&)>> finalizers;
  bool strict_file = true;

  const Entry* find(const std::string& key) const;
  Entry& add(Entry e);
};

}  // namespace detail

// Resolved parameter values.  Copyable; carries a pointer to the schema it was
// produced from so it can validate set() and drive print()/dump().
class Params {
public:
  Params() = default;

  template <typename T> T get(const std::string& key) const;
  template <typename T> T get_or(const std::string& key, T fallback) const;

  bool has(const std::string& key) const;      // resolved to a value
  bool was_set(const std::string& key) const;  // given explicitly by the user
  Source source(const std::string& key) const;

  template <typename T> void set(const std::string& key, T v);

  // the positional argument: the interpolator input file
  const std::string& input_file() const { return m_positional; }
  bool check_only() const { return m_check_only; }
  bool help_requested() const { return m_help; }

  void print(std::ostream& out) const;
  void print() const;
  void dump(const std::string& folder) const;              // folder/params.dat
  void dump(const std::string& folder, double t) const;    // folder/params_from_t<t>.dat

  const detail::SchemaImpl* schema() const { return m_schema.get(); }

private:
  friend class Schema;

  const Value& raw(const std::string& key) const;
  void store(const std::string& key, Value v, Source src);
  void set_checked(const std::string& key, Value v, Kind k);

  bool         as_bool(const std::string& key, const Value& v) const;
  long long    as_int (const std::string& key, const Value& v) const;
  double       as_real(const std::string& key, const Value& v) const;
  std::string  as_str (const std::string& key, const Value& v) const;

  [[noreturn]] void fail_range(const std::string& key, long long x,
                               const std::string& why) const;

  template <typename T> T narrow(const std::string& key, long long x) const;

  void write_to(const std::string& filename) const;

  std::shared_ptr<const detail::SchemaImpl> m_schema;
  std::map<std::string, Value> m_values;
  std::map<std::string, Source> m_sources;
  std::string m_positional;
  bool m_check_only = false;
  bool m_help = false;
};

class Schema {
public:
  explicit Schema(std::string app_name,
                  std::string positional_doc = "<interpolator_file>");

  // No default; an error if absent from both the checkpoint file and the CLI.
  template <typename T> Schema& require(std::string key, std::string doc);

  // Default used when the parameter is not given.
  template <typename T> Schema& opt(std::string key, T def, std::string doc);

  // Required only when `when` holds, evaluated after parsing so it can depend
  // on another parameter's value.  `why` is quoted back in the error message.
  template <typename T> Schema& require_if(std::string key, Pred when,
                                           std::string why, std::string doc);

  // As above, but with a default so the key can still be read when the
  // condition does not hold.  It must then be given explicitly when it does.
  template <typename T> Schema& require_if(std::string key, T def, Pred when,
                                           std::string why, std::string doc);

  // Computed by the program: rejected on the command line, read back from a
  // checkpoint file, written to dumps.  This is what Lx/Ly/Lz are.
  template <typename T> Schema& runtime(std::string key, T init, std::string doc);

  Schema& choices(const std::string& key, std::vector<std::string> allowed);
  Schema& token_choices(const std::string& key, const std::string& sep,
                        std::vector<std::string> allowed);
  Schema& check(Pred ok, std::string message);
  Schema& warn(Pred trigger, std::string message);
  Schema& finalize(std::function<void(Params&)> f);

  // Unknown keys in a checkpoint file: error (true, default) or warn (false).
  Schema& strict_file(bool on);

  Params parse(int argc, char* argv[]) const;
  Params parse(const std::vector<std::string>& args) const;  // args[0] is the program

  std::string help() const;
  const std::string& app() const;

  // Checks the schema itself: duplicate keys, defaults matching their declared
  // kind, and that every predicate only reads unconditionally-present keys.
  // Intended for unit tests -- throws ParamError on a malformed schema.
  void validate_self() const;

private:
  detail::Entry& entry(const std::string& key);
  std::shared_ptr<detail::SchemaImpl> m_impl;
};

// Convenience wrappers for apps: print to stderr and exit non-zero rather than
// propagating the exception.  parse_or_abort (MPI-safe) lives in ParamsMPI.hpp.
Params parse_or_exit(const Schema& s, int argc, char* argv[]);
void report(const ParamError& e, std::ostream& out);

// Template definitions

template <typename T>
Schema& Schema::require(std::string key, std::string doc) {
  detail::Entry e;
  e.key = std::move(key);
  e.kind = detail::kind_of<T>();
  e.doc = std::move(doc);
  e.required = true;
  m_impl->add(std::move(e));
  return *this;
}

template <typename T>
Schema& Schema::opt(std::string key, T def, std::string doc) {
  detail::Entry e;
  e.key = std::move(key);
  e.kind = detail::kind_of<T>();
  e.doc = std::move(doc);
  e.def = detail::to_value<T>(std::move(def));
  m_impl->add(std::move(e));
  return *this;
}

template <typename T>
Schema& Schema::require_if(std::string key, Pred when, std::string why,
                           std::string doc) {
  detail::Entry e;
  e.key = std::move(key);
  e.kind = detail::kind_of<T>();
  e.doc = std::move(doc);
  e.required_when = std::move(when);
  e.required_why = std::move(why);
  m_impl->add(std::move(e));
  return *this;
}

template <typename T>
Schema& Schema::require_if(std::string key, T def, Pred when, std::string why,
                           std::string doc) {
  detail::Entry e;
  e.key = std::move(key);
  e.kind = detail::kind_of<T>();
  e.doc = std::move(doc);
  e.def = detail::to_value<T>(std::move(def));
  e.required_when = std::move(when);
  e.required_why = std::move(why);
  m_impl->add(std::move(e));
  return *this;
}

template <typename T>
Schema& Schema::runtime(std::string key, T init, std::string doc) {
  detail::Entry e;
  e.key = std::move(key);
  e.kind = detail::kind_of<T>();
  e.doc = std::move(doc);
  e.is_runtime = true;
  e.def = detail::to_value<T>(std::move(init));
  m_impl->add(std::move(e));
  return *this;
}

template <typename T>
T Params::narrow(const std::string& key, long long x) const {
  if constexpr (std::is_unsigned_v<T>) {
    if (x < 0)
      fail_range(key, x, "must not be negative");
    const auto ux = static_cast<unsigned long long>(x);
    if (ux > static_cast<unsigned long long>(std::numeric_limits<T>::max()))
      fail_range(key, x, "is too large");
    return static_cast<T>(ux);
  } else {
    if (x < static_cast<long long>(std::numeric_limits<T>::min()) ||
        x > static_cast<long long>(std::numeric_limits<T>::max()))
      fail_range(key, x, "is out of range");
    return static_cast<T>(x);
  }
}

template <typename T>
T Params::get(const std::string& key) const {
  const Value& v = raw(key);
  if constexpr (std::is_same_v<T, bool>)             return as_bool(key, v);
  else if constexpr (std::is_same_v<T, std::string>) return as_str(key, v);
  else if constexpr (std::is_floating_point_v<T>)    return static_cast<T>(as_real(key, v));
  else if constexpr (std::is_integral_v<T>)          return narrow<T>(key, as_int(key, v));
  else static_assert(detail::always_false<T>, "unsupported parameter type");
}

template <typename T>
T Params::get_or(const std::string& key, T fallback) const {
  return has(key) ? get<T>(key) : fallback;
}

template <typename T>
void Params::set(const std::string& key, T v) {
  set_checked(key, detail::to_value<T>(std::move(v)), detail::kind_of<T>());
}

}  // namespace partrac

#endif
