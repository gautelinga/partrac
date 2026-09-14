#include "Params.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>

namespace partrac {

// Small helpers

namespace {

std::string trim(const std::string& s) {
  const char* ws = " \t\n\r\f\v";
  const auto b = s.find_first_not_of(ws);
  if (b == std::string::npos) return "";
  const auto e = s.find_last_not_of(ws);
  return s.substr(b, e - b + 1);
}

// Parse a real, requiring the whole string to be consumed
bool try_parse_real(const std::string& s, double& out) {
  if (s.empty()) return false;
  errno = 0;
  char* end = nullptr;
  const double d = std::strtod(s.c_str(), &end);
  if (end == s.c_str()) return false;
  while (*end != '\0' && std::isspace(static_cast<unsigned char>(*end))) ++end;
  if (*end != '\0') return false;
  if (!std::isfinite(d)) return false;
  out = d;
  return true;
}

// Via double, to allow scientific notation
bool try_parse_int(const std::string& s, long long& out) {
  double d;
  if (!try_parse_real(s, d)) return false;
  if (std::trunc(d) != d) return false;
  if (std::fabs(d) > 9007199254740992.0) return false;  // 2^53
  out = static_cast<long long>(d);
  return true;
}

bool try_parse_bool(const std::string& s, bool& out) {
  if (s == "true"  || s == "True"  || s == "TRUE"  || s == "1" ||
      s == "yes"   || s == "on")  { out = true;  return true; }
  if (s == "false" || s == "False" || s == "FALSE" || s == "0" ||
      s == "no"    || s == "off") { out = false; return true; }
  return false;
}

// Damerau-Levenshtein, so a transposition such as Nwr for Nrw costs 1
std::size_t edit_distance(const std::string& a, const std::string& b) {
  const std::size_t n = a.size(), m = b.size();
  std::vector<std::vector<std::size_t>> d(n + 1, std::vector<std::size_t>(m + 1));
  for (std::size_t i = 0; i <= n; ++i) d[i][0] = i;
  for (std::size_t j = 0; j <= m; ++j) d[0][j] = j;
  for (std::size_t i = 1; i <= n; ++i) {
    for (std::size_t j = 1; j <= m; ++j) {
      const std::size_t cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
      d[i][j] = std::min({d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + cost});
      if (i > 1 && j > 1 && a[i - 1] == b[j - 2] && a[i - 2] == b[j - 1])
        d[i][j] = std::min(d[i][j], d[i - 2][j - 2] + 1);
    }
  }
  return d[n][m];
}

// Nearest registered key, if it is close enough to be worth suggesting.
std::string suggest(const std::string& key,
                    const std::vector<detail::Entry>& entries) {
  std::string best;
  std::size_t best_d = std::string::npos;
  for (const auto& e : entries) {
    const std::size_t d = edit_distance(key, e.key);
    if (d < best_d) { best_d = d; best = e.key; }
  }
  // short keys need a tight limit, or 'nx' suggests 'Dm'
  const std::size_t limit = std::max<std::size_t>(1, std::min<std::size_t>(3, key.size() / 3));
  return (best_d <= limit) ? best : "";
}

bool convert(const detail::Entry& e, const std::string& raw, Value& out,
             std::string& err) {
  switch (e.kind) {
    case Kind::Bool: {
      bool b;
      if (!try_parse_bool(raw, b)) {
        err = "'" + raw + "' is not a boolean (use true/false)";
        return false;
      }
      out = Value(b);
      return true;
    }
    case Kind::Int: {
      long long i;
      if (!try_parse_int(raw, i)) {
        err = "'" + raw + "' is not an integer";
        return false;
      }
      out = Value(i);
      return true;
    }
    case Kind::Real: {
      double d;
      if (!try_parse_real(raw, d)) {
        err = "'" + raw + "' is not a number";
        return false;
      }
      out = Value(d);
      return true;
    }
    case Kind::Str:
      out = Value(raw);
      return true;
  }
  err = "internal: unknown kind";
  return false;
}

}  // namespace

// Free functions

std::string to_string(Kind k) {
  switch (k) {
    case Kind::Bool: return "bool";
    case Kind::Int:  return "int";
    case Kind::Real: return "real";
    case Kind::Str:  return "string";
  }
  return "?";
}

std::string to_string(Source s) {
  switch (s) {
    case Source::Unset:   return "unset";
    case Source::Default: return "default";
    case Source::File:    return "file";
    case Source::Cmdline: return "cmdline";
    case Source::Runtime: return "runtime";
  }
  return "?";
}

// Full precision, so a dump round-trips
std::string value_to_string(const Value& v) {
  std::ostringstream ss;
  if (std::holds_alternative<bool>(v)) {
    ss << (std::get<bool>(v) ? "true" : "false");
  } else if (std::holds_alternative<long long>(v)) {
    ss << std::get<long long>(v);
  } else if (std::holds_alternative<double>(v)) {
    ss << std::setprecision(std::numeric_limits<double>::max_digits10)
       << std::get<double>(v);
  } else {
    ss << std::get<std::string>(v);
  }
  return ss.str();
}

namespace {
std::string compose(const std::string& app,
                    const std::vector<std::string>& problems) {
  std::ostringstream ss;
  ss << app << ": " << problems.size()
     << (problems.size() == 1 ? " problem" : " problems")
     << " with the parameters";
  for (const auto& p : problems) ss << "\n  " << p;
  return ss.str();
}
}  // namespace

ParamError::ParamError(std::string app, std::vector<std::string> problems)
    : std::runtime_error(compose(app, problems)),
      m_app(std::move(app)),
      m_problems(std::move(problems)) {}

void report(const ParamError& e, std::ostream& out) {
  out << e.what() << "\n\n"
      << "Run './" << e.app() << " --help' for the parameter list." << std::endl;
}

// SchemaImpl

namespace detail {

const Entry* SchemaImpl::find(const std::string& key) const {
  const auto it = index.find(key);
  return (it == index.end()) ? nullptr : &entries[it->second];
}

Entry& SchemaImpl::add(Entry e) {
  if (index.count(e.key))
    throw ParamError(app, {"parameter '" + e.key + "' registered twice"});
  index[e.key] = entries.size();
  entries.push_back(std::move(e));
  return entries.back();
}

}  // namespace detail

// Params

const Value& Params::raw(const std::string& key) const {
  const auto it = m_values.find(key);
  if (it != m_values.end()) return it->second;

  if (m_schema && !m_schema->find(key))
    throw ParamError(m_schema ? m_schema->app : "params",
                     {"parameter '" + key + "' is not declared by this app"});
  throw ParamError(m_schema ? m_schema->app : "params",
                   {"parameter '" + key + "' is not set for this configuration"});
}

void Params::store(const std::string& key, Value v, Source src) {
  m_values[key] = std::move(v);
  m_sources[key] = src;
}

void Params::set_checked(const std::string& key, Value v, Kind k) {
  const detail::Entry* e = m_schema ? m_schema->find(key) : nullptr;
  if (!e)
    throw ParamError(m_schema ? m_schema->app : "params",
                     {"cannot set '" + key + "': not declared by this app"});
  // Widening an integer literal into a real parameter is intentional and safe.
  if (e->kind != k && !(e->kind == Kind::Real && k == Kind::Int))
    throw ParamError(m_schema->app,
                     {"cannot set '" + key + "': declared as " + to_string(e->kind) +
                      ", given " + to_string(k)});
  if (e->kind == Kind::Real && k == Kind::Int)
    v = Value(static_cast<double>(std::get<long long>(v)));
  store(key, std::move(v), Source::Runtime);
}

bool Params::has(const std::string& key) const {
  return m_values.count(key) != 0;
}

bool Params::was_set(const std::string& key) const {
  const auto it = m_sources.find(key);
  if (it == m_sources.end()) return false;
  return it->second == Source::File || it->second == Source::Cmdline;
}

Source Params::source(const std::string& key) const {
  const auto it = m_sources.find(key);
  return (it == m_sources.end()) ? Source::Unset : it->second;
}

bool Params::as_bool(const std::string& key, const Value& v) const {
  if (!std::holds_alternative<bool>(v))
    throw ParamError(m_schema ? m_schema->app : "params",
                     {"parameter '" + key + "' is not a boolean"});
  return std::get<bool>(v);
}

long long Params::as_int(const std::string& key, const Value& v) const {
  if (!std::holds_alternative<long long>(v))
    throw ParamError(m_schema ? m_schema->app : "params",
                     {"parameter '" + key + "' is not an integer"});
  return std::get<long long>(v);
}

double Params::as_real(const std::string& key, const Value& v) const {
  if (std::holds_alternative<double>(v))    return std::get<double>(v);
  if (std::holds_alternative<long long>(v)) return static_cast<double>(std::get<long long>(v));
  throw ParamError(m_schema ? m_schema->app : "params",
                   {"parameter '" + key + "' is not a number"});
}

std::string Params::as_str(const std::string& key, const Value& v) const {
  if (!std::holds_alternative<std::string>(v))
    throw ParamError(m_schema ? m_schema->app : "params",
                     {"parameter '" + key + "' is not a string"});
  return std::get<std::string>(v);
}

void Params::fail_range(const std::string& key, long long x,
                        const std::string& why) const {
  throw ParamError(m_schema ? m_schema->app : "params",
                   {"parameter '" + key + "' = " + std::to_string(x) + " " + why});
}

void Params::print(std::ostream& out) const {
  if (!m_schema) return;
  std::size_t w = 0;
  for (const auto& e : m_schema->entries) w = std::max(w, e.key.size());
  for (const auto& e : m_schema->entries) {
    const auto it = m_values.find(e.key);
    if (it == m_values.end()) continue;
    out << "  " << std::left << std::setw(static_cast<int>(w)) << e.key
        << " = " << value_to_string(it->second)
        << "  [" << to_string(source(e.key)) << "]" << std::endl;
  }
}

void Params::print() const { print(std::cout); }

void Params::write_to(const std::string& filename) const {
  std::ofstream f(filename);
  if (!f)
    throw ParamError(m_schema ? m_schema->app : "params",
                     {"could not open '" + filename + "' for writing"});
  if (m_schema) {
    for (const auto& e : m_schema->entries) {
      const auto it = m_values.find(e.key);
      if (it == m_values.end()) continue;
      f << e.key << "=" << value_to_string(it->second) << "\n";
    }
  }
}

void Params::dump(const std::string& folder) const {
  write_to(folder + "/params.dat");
}

void Params::dump(const std::string& folder, const double t) const {
  write_to(folder + "/params_from_t" + std::to_string(t) + ".dat");
}

// Schema

Schema::Schema(std::string app_name, std::string positional_doc)
    : m_impl(std::make_shared<detail::SchemaImpl>()) {
  m_impl->app = std::move(app_name);
  m_impl->positional_doc = std::move(positional_doc);
}

const std::string& Schema::app() const { return m_impl->app; }

detail::Entry& Schema::entry(const std::string& key) {
  const auto it = m_impl->index.find(key);
  if (it == m_impl->index.end())
    throw ParamError(m_impl->app, {"parameter '" + key + "' is not registered"});
  return m_impl->entries[it->second];
}

Schema& Schema::token_choices(const std::string& key, const std::string& sep,
                              std::vector<std::string> allowed) {
  entry(key).choices_sep = sep;
  return choices(key, std::move(allowed));
}

Schema& Schema::choices(const std::string& key, std::vector<std::string> allowed) {
  entry(key).choices = std::move(allowed);
  return *this;
}

Schema& Schema::check(Pred ok, std::string message) {
  m_impl->checks.push_back({std::move(ok), std::move(message)});
  return *this;
}

Schema& Schema::warn(Pred trigger, std::string message) {
  m_impl->warns.push_back({std::move(trigger), std::move(message)});
  return *this;
}

Schema& Schema::finalize(std::function<void(Params&)> f) {
  m_impl->finalizers.push_back(std::move(f));
  return *this;
}

Schema& Schema::strict_file(bool on) {
  m_impl->strict_file = on;
  return *this;
}

Params Schema::parse(int argc, char* argv[]) const {
  std::vector<std::string> args;
  args.reserve(static_cast<std::size_t>(argc));
  for (int i = 0; i < argc; ++i) args.emplace_back(argv[i]);
  return parse(args);
}

Params Schema::parse(const std::vector<std::string>& args) const {
  std::vector<std::string> problems;
  Params p;
  p.m_schema = m_impl;

  // --- 1. tokenize -------------------------------------------------------
  std::map<std::string, std::string> cmd;
  bool positional_seen = false;
  for (std::size_t i = 1; i < args.size(); ++i) {
    const std::string a = trim(args[i]);
    if (a.empty()) continue;
    if (a == "--check" || a == "--dry-run") { p.m_check_only = true; continue; }
    if (a == "--help" || a == "-h") { p.m_help = true; continue; }
    const auto eq = a.find('=');
    if (eq == std::string::npos) {
      // First bare argument is the interpolator file
      if (!positional_seen && a.rfind("-", 0) != 0) {
        positional_seen = true;
        p.m_positional = a;
        continue;
      }
      problems.push_back("unrecognised argument '" + a + "' (expected key=value)");
      continue;
    }
    const std::string key = trim(a.substr(0, eq));
    const std::string val = trim(a.substr(eq + 1));
    if (cmd.count(key))
      problems.push_back("parameter '" + key + "' given more than once");
    cmd[key] = val;
  }

  // --help short-circuits everything: the caller prints the schema and exits,
  // so it must work even when required parameters are absent.
  if (p.m_help) return p;

  if (!positional_seen && !m_impl->positional_doc.empty())
    problems.push_back("no input file given, expected " + m_impl->positional_doc);

  // --- 2. validate command-line keys -------------------------------------
  for (const auto& kv : cmd) {
    const detail::Entry* e = m_impl->find(kv.first);
    if (!e) {
      const std::string s = suggest(kv.first, m_impl->entries);
      problems.push_back("unknown parameter '" + kv.first + "'" +
                         (s.empty() ? "" : "; did you mean '" + s + "'?"));
    } else if (e->is_runtime) {
      problems.push_back("parameter '" + kv.first +
                         "' is computed by the program and cannot be given here");
    }
  }

  // --- 3. defaults -------------------------------------------------------
  for (const auto& e : m_impl->entries)
    if (e.def) p.store(e.key, *e.def, Source::Default);

  // --- 4. checkpoint file, if restarting ---------------------------------
  const auto rf = cmd.find("restart_folder");
  if (rf != cmd.end() && !rf->second.empty() && m_impl->find("restart_folder")) {
    const std::string path = rf->second + "/Checkpoints/params.dat";
    std::ifstream in(path);
    if (!in) {
      problems.push_back("restart file '" + path + "' does not exist");
    } else {
      std::set<std::string> seen;
      for (std::string line; std::getline(in, line);) {
        const auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        const std::string key = trim(line.substr(0, eq));
        const std::string val = trim(line.substr(eq + 1));
        if (!seen.insert(key).second)
          problems.push_back("parameter '" + key + "' given more than once in " + path);
        const detail::Entry* e = m_impl->find(key);
        if (!e) {
          const std::string msg = "unknown parameter '" + key + "' in " + path;
          if (m_impl->strict_file) problems.push_back(msg);
          else std::cerr << "Warning: " << msg << ", ignoring." << std::endl;
          continue;
        }
        Value v;
        std::string err;
        if (!convert(*e, val, v, err)) problems.push_back(key + ": " + err +
                                                          " (in " + path + ")");
        else p.store(key, std::move(v), Source::File);
      }
    }
  }

  // --- 5. command line wins ----------------------------------------------
  for (const auto& kv : cmd) {
    const detail::Entry* e = m_impl->find(kv.first);
    if (!e || e->is_runtime) continue;  // already reported
    Value v;
    std::string err;
    if (!convert(*e, kv.second, v, err)) problems.push_back(kv.first + ": " + err);
    else p.store(kv.first, std::move(v), Source::Cmdline);
  }

  // --- 6. required, including conditionals --------------------------------
  for (const auto& e : m_impl->entries) {
    if (e.required) {
      if (p.has(e.key)) continue;
      problems.push_back("missing required parameter '" + e.key + "'" +
                         (e.doc.empty() ? "" : "\n      " + e.key + ": " + e.doc));
    } else if (e.required_when) {
      // a default does not satisfy a condition -- it must be given explicitly
      if (p.was_set(e.key)) continue;
      bool needed = false;
      try {
        needed = e.required_when(p);
      } catch (const ParamError& ex) {
        problems.push_back("could not evaluate the condition for '" + e.key +
                           "': " + ex.problems().front());
        continue;
      }
      if (needed)
        problems.push_back("missing required parameter '" + e.key + "'" +
                           (e.required_why.empty() ? "" : "\n      required because " +
                                                          e.required_why) +
                           (e.doc.empty() ? "" : "\n      " + e.key + ": " + e.doc));
    }
  }

  // Stop here if anything is missing: finalizers and checks below assume a
  // fully resolved set and would otherwise throw on the first absent key.
  if (!problems.empty()) throw ParamError(m_impl->app, std::move(problems));

  // --- 7. finalize (the clamps), then value constraints -------------------
  for (const auto& f : m_impl->finalizers) f(p);

  for (const auto& e : m_impl->entries) {
    if (e.choices.empty() || !p.has(e.key)) continue;
    const std::string v = p.get<std::string>(e.key);
    const std::string tok = e.choices_sep.empty()
                          ? v : v.substr(0, v.find(e.choices_sep));
    if (std::find(e.choices.begin(), e.choices.end(), tok) == e.choices.end()) {
      std::ostringstream ss;
      ss << "parameter '" << e.key << "' = '" << v << "' is not one of:";
      for (const auto& c : e.choices)
        ss << " " << c << e.choices_sep << (e.choices_sep.empty() ? "" : "*");
      problems.push_back(ss.str());
    }
  }

  for (const auto& c : m_impl->checks)
    if (!c.pred(p)) problems.push_back(c.message);

  if (!problems.empty()) throw ParamError(m_impl->app, std::move(problems));

  for (const auto& w : m_impl->warns)
    if (w.pred(p)) std::cerr << "Warning: " << w.message << std::endl;

  return p;
}

std::string Schema::help() const {
  std::ostringstream ss;
  ss << "Usage: " << m_impl->app << " " << m_impl->positional_doc
     << " [key=value ...]\n\n";
  std::size_t w = 0;
  for (const auto& e : m_impl->entries) w = std::max(w, e.key.size());

  auto emit = [&](const std::string& title, auto pick) {
    bool any = false;
    for (const auto& e : m_impl->entries) {
      if (!pick(e)) continue;
      if (!any) { ss << title << "\n"; any = true; }
      ss << "  " << std::left << std::setw(static_cast<int>(w)) << e.key
         << "  " << to_string(e.kind);
      if (!e.choices.empty()){
        ss << " {";
        for (std::size_t i = 0; i < e.choices.size(); ++i)
          ss << (i ? ", " : "") << e.choices[i] << e.choices_sep
             << (e.choices_sep.empty() ? "" : "*");
        ss << "}";
      }
      if (e.def) ss << " (default " << value_to_string(*e.def) << ")";
      if (!e.doc.empty()) ss << " -- " << e.doc;
      if (!e.required_why.empty()) ss << " [required when " << e.required_why << "]";
      ss << "\n";
    }
    if (any) ss << "\n";
  };

  emit("Required:", [](const detail::Entry& e) { return e.required; });
  emit("Required in some configurations:",
       [](const detail::Entry& e) { return static_cast<bool>(e.required_when); });
  emit("Optional:", [](const detail::Entry& e) {
    return !e.required && !e.required_when && !e.is_runtime;
  });
  emit("Computed by the program (not accepted here):",
       [](const detail::Entry& e) { return e.is_runtime; });
  return ss.str();
}

void Schema::validate_self() const {
  std::vector<std::string> problems;

  for (const auto& e : m_impl->entries) {
    if (e.def) {
      const Kind actual =
          std::holds_alternative<bool>(*e.def)      ? Kind::Bool :
          std::holds_alternative<long long>(*e.def) ? Kind::Int  :
          std::holds_alternative<double>(*e.def)    ? Kind::Real : Kind::Str;
      if (actual != e.kind)
        problems.push_back("'" + e.key + "' declared as " + to_string(e.kind) +
                           " but its default is " + to_string(actual));
    }
    if (e.required && e.def)
      problems.push_back("'" + e.key + "' is required but also has a default");
    if (e.required && e.required_when)
      problems.push_back("'" + e.key + "' is both required and conditional");
    if (!e.choices.empty() && e.kind != Kind::Str)
      problems.push_back("'" + e.key + "' has choices but is not a string");
  }

  // Predicates may only read keys that are always present, otherwise they can
  // throw during parse.  Probe them against a Params holding exactly those.
  Params probe;
  probe.m_schema = m_impl;
  for (const auto& e : m_impl->entries)
    if (e.def) probe.store(e.key, *e.def, Source::Default);
    else if (e.required) probe.store(e.key, [&] {
      switch (e.kind) {
        case Kind::Bool: return Value(false);
        case Kind::Int:  return Value(0LL);
        case Kind::Real: return Value(0.0);
        default:         return Value(std::string());
      }
    }(), Source::Default);

  auto probe_pred = [&](const Pred& pred, const std::string& what) {
    if (!pred) return;
    try {
      (void)pred(probe);
    } catch (const ParamError& ex) {
      problems.push_back(what + " reads a parameter that may not be set: " +
                         ex.problems().front());
    }
  };

  for (const auto& e : m_impl->entries)
    probe_pred(e.required_when, "the condition for '" + e.key + "'");
  for (const auto& c : m_impl->checks) probe_pred(c.pred, "check '" + c.message + "'");
  for (const auto& w : m_impl->warns) probe_pred(w.pred, "warning '" + w.message + "'");

  if (!problems.empty()) throw ParamError(m_impl->app, std::move(problems));
}

Params parse_or_exit(const Schema& s, int argc, char* argv[]) {
  try {
    Params p = s.parse(argc, argv);
    if (p.help_requested()) {
      std::cout << s.help();
      std::exit(0);
    }
    return p;
  } catch (const ParamError& e) {
    report(e, std::cerr);
    std::exit(2);
  }
}

}  // namespace partrac
