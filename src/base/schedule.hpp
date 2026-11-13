#ifndef __SCHEDULE_HPP
#define __SCHEDULE_HPP

// The schedules of dump_intv, stat_intv and checkpoint_intv: comma-separated
// segments in absolute time t,
//   X         every X to the end of the run
//   X:U       every X while below U, from the previous segment's U (the run's
//             start time for the first segment): t = start + n X; and U itself
//   logN:S:U  N events a decade from t = S to U: t_k = S 10^(k/N), k = 0, 1, ...;
//             U is an event when a t_k equals it; else the last t_k below U
//             ends the segment
//   logN:U    the same from the previous segment's U
// Each event time is taken at the nearest step, a time half way between two
// steps at the later one, and the times on one step are one event, the switch
// between two segments too. No event comes after the last segment's U. The
// run's start time t0 (a march's xn0) is an event too, where the first segment
// starts later: the first dump holds the start positions; on a restart t0 is
// the original run's. A uniform X is floored at one step. Event times are a
// function of t, so a restart continues the schedule.
// A single bare number X is every floor(X/dt) steps, at least one, from step 0
// (at_interval); 0 is off.

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include "intervals.hpp"

struct ScheduleSegment {
  bool log = false;
  double X = 0.;         // uniform: the interval
  long long N = 0;       // log: events a decade
  bool has_S = false;    // log: S given, else the previous segment's U
  double S = 0.;
  bool has_U = false;    // until; a bare X has none
  double U = 0.;
  std::string text;      // as written
};

struct Schedule {
  std::vector<ScheduleSegment> segs;
  // A single bare number: the step interval of at_interval
  bool plain() const { return segs.size() == 1 && !segs[0].log && !segs[0].has_U; }
  bool off() const { return plain() && !(segs[0].X > 0.); }
};

// The most events a decade
constexpr long long schedule_max_N = 1000;

// What is wrong with a schedule; each has one rule, schedule_rule
enum class ScheduleError { None, Syntax, Negative, NonPositive, LogN, LogS, LogStart,
                           BareNotLast, Order };

struct ScheduleProblem {
  ScheduleError code = ScheduleError::None;
  std::string where;     // the segment, or the schedule
  explicit operator bool() const { return code != ScheduleError::None; }
};

inline std::string schedule_rule(const ScheduleError e){
  switch (e){
    case ScheduleError::None: return "";
    case ScheduleError::Syntax:
      return "not a schedule: comma-separated segments X, X:U, logN:S:U or logN:U, "
             "with X, S and U numbers";
    case ScheduleError::Negative: return "intervals cannot be negative";
    case ScheduleError::NonPositive:
      return "an interval X in a schedule must be positive; 0 alone turns the output off";
    case ScheduleError::LogN:
      return "N in logN must be an integer from 1 to " + std::to_string(schedule_max_N);
    case ScheduleError::LogS:
      return "a logN segment must start above 0: S, or the previous segment's U";
    case ScheduleError::LogStart:
      return "logN:U starts at the previous segment's U, so it cannot be the first segment";
    case ScheduleError::BareNotLast:
      return "a bare interval X runs to the end of the run, so only the last segment can be one";
    case ScheduleError::Order:
      return "the untils U must increase strictly along the schedule, and lie above "
             "each logN:S:U's S";
  }
  return "";
}

// key: where: rule; empty if no problem
inline std::string schedule_message(const std::string& key, const ScheduleProblem& p){
  return p ? key + ": " + p.where + ": " + schedule_rule(p.code) : "";
}

namespace schedule_detail {

// A real, the whole string, finite: as the parameters parse one
inline bool parse_real(const std::string& s, double& out){
  if (s.empty()) return false;
  char* end = nullptr;
  const double d = std::strtod(s.c_str(), &end);
  if (end == s.c_str() || *end != '\0' || !std::isfinite(d)) return false;
  out = d;
  return true;
}

inline std::string strip(const std::string& s){
  const char* ws = " \t\n\r\f\v";
  const auto b = s.find_first_not_of(ws);
  if (b == std::string::npos) return "";
  return s.substr(b, s.find_last_not_of(ws) - b + 1);
}

inline std::vector<std::string> split(const std::string& s, const char sep){
  std::vector<std::string> out;
  std::string::size_type b = 0;
  while (true){
    const auto e = s.find(sep, b);
    out.push_back(strip(s.substr(b, e == std::string::npos ? std::string::npos : e - b)));
    if (e == std::string::npos) return out;
    b = e + 1;
  }
}

inline ScheduleProblem problem(const ScheduleError e, const std::string& where){
  return ScheduleProblem{e, "'" + where + "'"};
}

// Round-off in times of size a and b, and a sliver of a step
inline double tol(const double step, const double a, const double b){
  return 1e-9*step + 1e-13*(std::abs(a) + std::abs(b));
}

}  // namespace schedule_detail

// text into s; the problem if it is not a schedule
inline ScheduleProblem parse_schedule(const std::string& text, Schedule& s){
  using namespace schedule_detail;
  s.segs.clear();
  for (const std::string& seg : split(text, ',')){
    ScheduleSegment g;
    g.text = seg;
    const std::vector<std::string> parts = split(seg, ':');
    if (seg.empty())
      return problem(ScheduleError::Syntax, text);
    if (parts[0].rfind("log", 0) == 0){
      double n;
      if (parts.size() < 2 || parts.size() > 3 || !parse_real(parts[0].substr(3), n))
        return problem(ScheduleError::Syntax, seg);
      if (!(n >= 1.) || n != std::floor(n) || n > double(schedule_max_N))
        return problem(ScheduleError::LogN, seg);
      g.log = true;
      g.N = static_cast<long long>(n);
      g.has_S = parts.size() == 3;
      g.has_U = true;
      if ((g.has_S && !parse_real(parts[1], g.S)) || !parse_real(parts.back(), g.U))
        return problem(ScheduleError::Syntax, seg);
    }
    else {
      g.has_U = parts.size() == 2;
      if (parts.size() > 2 || !parse_real(parts[0], g.X) || (g.has_U && !parse_real(parts[1], g.U)))
        return problem(ScheduleError::Syntax, seg);
    }
    s.segs.push_back(g);
  }

  if (s.plain()){
    if (s.segs[0].X < 0.)
      return problem(ScheduleError::Negative, text);
    return {};
  }
  for (std::size_t i = 0; i < s.segs.size(); ++i){
    const ScheduleSegment& g = s.segs[i];
    if (!g.has_U && i + 1 < s.segs.size())
      return problem(ScheduleError::BareNotLast, g.text);
    if (!g.log && !(g.X > 0.))
      return problem(ScheduleError::NonPositive, g.text);
    if (g.log && !g.has_S && i == 0)
      return problem(ScheduleError::LogStart, g.text);
    if (g.log && !((g.has_S ? g.S : s.segs[i - 1].U) > 0.))
      return problem(ScheduleError::LogS, g.text);
    if (g.log && g.has_S && !(g.U > g.S))
      return problem(ScheduleError::Order, g.text);
    if (g.has_U && i > 0 && !(g.U > s.segs[i - 1].U))
      return problem(ScheduleError::Order, g.text);
  }
  return {};
}

// The steps' times: t_start at step it_start, a step dt
struct StepGrid {
  double t_start = 0.;
  long long it_start = 0;
  double dt = 1.;
  double time(const long long k) const { return t_start + (k - it_start)*dt; }
  // The nearest step; half way, the later, whichever way round-off went
  long long snap(const double t) const {
    return it_start + static_cast<long long>(std::floor((t - t_start)/dt + 0.5 + 1e-6));
  }
};

// The steps a schedule has events at, on a grid; events before the grid's
// start are passed (a restart's, which may have another step)
class EventClock {
public:
  static constexpr long long never = LLONG_MAX;

  EventClock() = default;   // no events

  // s parsed; start is the run's start time t0
  EventClock(const Schedule& s, const StepGrid& grid, const double start)
    : grid_(grid), start_(start) {
    using namespace schedule_detail;
    if (s.off() || s.segs.empty())
      return;
    if (s.plain()){
      plain_ = true;
      intv_ = s.segs[0].X;
      return;
    }
    const double step = grid.dt;
    for (std::size_t i = 0; i < s.segs.size(); ++i){
      const ScheduleSegment& g = s.segs[i];
      Piece p;
      p.log = g.log;
      p.U = g.U;
      p.ends_at_U = g.has_U;
      if (!g.log){
        p.a = i == 0 ? start : s.segs[i - 1].U;
        p.X = std::max(g.X, step);
        if (g.has_U){
          // a + n X below U, then U
          const double t = tol(step, g.U, p.a);
          long long n = std::max(-1LL, static_cast<long long>(std::ceil((g.U - p.a)/p.X)) - 1);
          while (n >= 0 && p.raw(n) >= g.U - t) --n;
          while (p.raw(n + 1) < g.U - t) ++n;
          p.J = n + 1;
        }
      }
      else {
        p.a = g.has_S ? g.S : s.segs[i - 1].U;
        p.N = g.N;
        const double t = tol(step, g.U, p.a);
        long long J = std::max(0LL, static_cast<long long>(std::floor(p.N*std::log10(g.U/p.a))));
        while (p.raw(J + 1) <= g.U + t) ++J;
        while (J > 0 && p.raw(J) > g.U + t) --J;
        p.J = J;
        p.ends_at_U = std::abs(p.raw(J) - g.U) <= t;
      }
      pieces_.push_back(p);
    }
    // The run's start
    Piece t0;
    t0.a = start;
    t0.X = 1.;
    t0.J = 0;
    pieces_.push_back(t0);
  }

  bool on() const { return plain_ || !pieces_.empty(); }
  bool plain() const { return plain_; }

  // Whether step it has an event; it rises along a run
  bool at(const long long it){
    if (plain_)
      return at_interval(static_cast<Uint>(it), intv_, grid_.dt);
    if (pieces_.empty())
      return false;
    if (it < from_ || it > next_){
      from_ = it;
      next_ = next(it);
    }
    return it == next_;
  }

  // The first step from k on with an event, never if none; a schedule's
  long long next(const long long k) const {
    double t;
    return first(k, passed(), t);
  }

  // The events from the run's start on before step k: those at the steps
  // before it, and on a restart's grid the times before its start
  long long count_before(const long long k) const {
    if (plain_){
      const long long m = static_cast<long long>(std::min<Uint>(steps_per(intv_, grid_.dt), LLONG_MAX/2));
      return k <= 0 ? 0 : (k - 1)/m + 1;
    }
    const double t_min = start_ - schedule_detail::tol(grid_.dt, start_, 0.);
    long long c = 0;
    double t;
    for (long long s = first(grid_.snap(start_), t_min, t); s != never; s = first(s + 1, t_min, t)){
      if (!(s < k || (s == k && t < passed())))
        break;
      ++c;
    }
    return c;
  }

private:
  struct Piece {
    bool log = false;
    double a = 0.;          // uniform: start; log: S
    double X = 0.;          // uniform
    long long N = 0;        // log
    double U = 0.;
    long long J = -1;       // the last index; -1: none
    bool ends_at_U = false; // the last is U

    double raw(const long long j) const {
      return log ? a*std::pow(10., static_cast<double>(j)/static_cast<double>(N)) : a + j*X;
    }
    double time(const long long j) const { return (ends_at_U && j == J) ? U : raw(j); }

    // The first step from k on with an event at t_min or later, and its time
    long long next(const StepGrid& g, const long long k, const double t_min, double& t) const {
      const double tk = std::max(g.time(k), t_min);
      double guess = 0.;
      if (log)
        guess = tk > a ? std::floor(N*std::log10(tk/a)) : 0.;
      else
        guess = tk > a ? std::floor((tk - a)/X) : 0.;
      if (!(guess < 1e18))
        return EventClock::never;
      long long j = static_cast<long long>(guess);
      if (J >= 0) j = std::min(j, J);
      while (j > 0 && g.snap(time(j - 1)) >= k && time(j - 1) >= t_min) --j;
      while ((J < 0 || j <= J) && (g.snap(time(j)) < k || time(j) < t_min)) ++j;
      if (J >= 0 && j > J)
        return EventClock::never;
      t = time(j);
      return g.snap(t);
    }
  };

  // Times before the grid's start are passed
  double passed() const { return grid_.t_start - schedule_detail::tol(grid_.dt, grid_.t_start, 0.); }

  // The first step from k on with an event at t_min or later, and the
  // earliest time on it
  long long first(const long long k, const double t_min, double& t) const {
    long long best = never;
    t = 0.;
    for (const Piece& p : pieces_){
      double tp;
      const long long s = p.next(grid_, k, t_min, tp);
      if (s < best || (s == best && tp < t)){
        best = s;
        t = tp;
      }
    }
    return best;
  }

  StepGrid grid_;
  double start_ = 0.;
  bool plain_ = false;
  double intv_ = 0.;
  std::vector<Piece> pieces_;
  long long from_ = LLONG_MAX, next_ = LLONG_MIN;   // next_ is the first event from from_ on
};

#endif
