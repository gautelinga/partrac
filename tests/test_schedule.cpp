// The schedules of dump_intv, stat_intv and checkpoint_intv (schedule.hpp):
// parsing and its messages, the event steps on a grid whose step does not
// divide the intervals, the merging of times on one step, the end of a
// schedule, restarts, and the parameter files.
#include <catch2/catch.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "Params.hpp"
#include "intervals.hpp"
#include "schedule.hpp"
#include "AppParams.hpp"
#include "RunLoop.hpp"

namespace {

// Steps from 0 to T, with the given grid
std::vector<long long> steps(const std::string& text, const double dt, const double T,
                             const StepGrid* grid = nullptr, const long long from = 0,
                             const double start = 0.){
  Schedule s;
  REQUIRE_FALSE(parse_schedule(text, s));
  const StepGrid g = grid ? *grid : StepGrid{0., 0, dt};
  EventClock c(s, g, start);
  std::vector<long long> out;
  for (long long k = c.next(from); k != EventClock::never && g.time(k) <= T + dt/2; k = c.next(k + 1))
    out.push_back(k);
  return out;
}

// Times snapped to the nearest step of dt, half way the later, merged
std::vector<long long> snapped(const std::vector<double>& times, const double dt){
  std::vector<long long> out;
  for (const double t : times)
    out.push_back(static_cast<long long>(std::floor(t/dt + 0.5 + 1e-6)));
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

std::vector<long long> concat(std::vector<long long> a, const std::vector<long long>& b){
  for (const long long k : b)
    if (a.empty() || k > a.back()) a.push_back(k);
  return a;
}

ScheduleError error_of(const std::string& text){
  Schedule s;
  return parse_schedule(text, s).code;
}

// The interval keys and what they need
partrac::Schema mini(){
  partrac::Schema s("mini");
  add_run_params(s, "current time");
  s.require<double>("dt", "timestep");
  s.opt<double>("t0", 0.0, "start time");
  s.require<Uint>("Nrw", "particles");
  s.opt<Uint>("Nrw_max", 1, "most particles");
  add_intervals(s, "dt");
  return s;
}

partrac::Params parse_mini(std::vector<std::string> more){
  std::vector<std::string> args{"mini", "in.dat", "dt=0.003", "Nrw=1"};
  args.insert(args.end(), more.begin(), more.end());
  return mini().parse(args);
}

// The problems of a parse that must fail
std::vector<std::string> problems(const std::vector<std::string>& more){
  try {
    parse_mini(more);
  } catch (const partrac::ParamError& e){
    return e.problems();
  }
  FAIL("parsed");
  return {};
}

struct TempDir {
  std::filesystem::path path;
  TempDir(){
    path = std::filesystem::temp_directory_path()
           / ("partrac_schedule_" + std::to_string(std::random_device{}()));
    std::filesystem::create_directories(path / "Checkpoints");
  }
  ~TempDir(){ std::error_code ec; std::filesystem::remove_all(path, ec); }
};

}  // namespace

// --- parsing -------------------------------------------------------------------

TEST_CASE("Schedules parse into their segments", "[schedule]"){
  Schedule s;
  REQUIRE_FALSE(parse_schedule("10", s));
  CHECK(s.plain());
  CHECK_FALSE(s.off());
  CHECK(s.segs[0].X == 10.);

  REQUIRE_FALSE(parse_schedule("0", s));
  CHECK(s.off());

  REQUIRE_FALSE(parse_schedule("0.1:1, log5:100 ,10", s));
  REQUIRE(s.segs.size() == 3);
  CHECK_FALSE(s.plain());
  CHECK((!s.segs[0].log && s.segs[0].X == 0.1 && s.segs[0].has_U && s.segs[0].U == 1.));
  CHECK((s.segs[1].log && s.segs[1].N == 5 && !s.segs[1].has_S && s.segs[1].U == 100.));
  CHECK((!s.segs[2].log && s.segs[2].X == 10. && !s.segs[2].has_U));

  REQUIRE_FALSE(parse_schedule("log10:0.1:100,10", s));
  CHECK((s.segs[0].log && s.segs[0].N == 10 && s.segs[0].has_S && s.segs[0].S == 0.1 && s.segs[0].U == 100.));

  // An until on the last: the schedule ends there
  REQUIRE_FALSE(parse_schedule("1e-2:1", s));
  CHECK((s.segs.size() == 1 && !s.plain() && s.segs[0].X == 0.01));

  for (const std::string ok : {"0.1:1,1:10,10", "log10:0.1:100,10", "0.1:1,log5:100,10", "log5:1:50",
                               "0.5:2,log3:20,log7:200:1000,100", "0.3:1", "0.001:1", "1:0.5",
                               "log1000:1:10"})
    CHECK(error_of(ok) == ScheduleError::None);
}

TEST_CASE("A schedule that breaks a rule is refused with that rule", "[schedule]"){
  using E = ScheduleError;
  for (const std::string bad : {"", "abc", "0.1:1,,10", "1:2:3", "log5", "log5:1:2:3", "logx:1:10",
                                "0.1:abc", "0.1:1,", "nan", "inf:10", "log:1:10"})
    CHECK(error_of(bad) == E::Syntax);
  CHECK(error_of("-1") == E::Negative);
  for (const std::string bad : {"0:1,10", "0.1:1,0", "-0.1:1", "0.1:1,-10"})
    CHECK(error_of(bad) == E::NonPositive);
  for (const std::string bad : {"log0:1:10", "log2.5:1:10", "log-1:1:10", "log1001:1:10"})
    CHECK(error_of(bad) == E::LogN);
  // S at or below 0, given or the previous U
  for (const std::string bad : {"log5:0:10", "log5:-1:10", "0.5:-1,log5:10", "0.5:0,log5:10"})
    CHECK(error_of(bad) == E::LogS);
  CHECK(error_of("log5:10,1") == E::LogStart);
  CHECK(error_of("10,0.1:1") == E::BareNotLast);
  for (const std::string bad : {"0.1:1,1:0.5", "log5:10:5", "1:10,log5:20:10", "0.1:1,0.1:1"})
    CHECK(error_of(bad) == E::Order);
}

TEST_CASE("The parameters name the key, the segment and the rule a schedule breaks", "[schedule]"){
  const auto said = [](const std::vector<std::string>& more, const std::string& message){
    const auto p = problems(more);
    INFO(more[0]);
    REQUIRE(p.size() == 1);
    CHECK(p[0] == message);
  };
  using E = ScheduleError;
  said({"dump_intv=0.1:1,0.1:abc"}, "dump_intv: '0.1:abc': " + schedule_rule(E::Syntax));
  said({"dump_intv=0.1:1,,10"}, "dump_intv: '0.1:1,,10': " + schedule_rule(E::Syntax));
  said({"stat_intv=-1"}, "stat_intv: '-1': intervals cannot be negative");
  said({"checkpoint_intv=0:1,10"}, "checkpoint_intv: '0:1': " + schedule_rule(E::NonPositive));
  said({"dump_intv=0.1:1,log0:10"}, "dump_intv: 'log0:10': " + schedule_rule(E::LogN));
  said({"dump_intv=log5:0:10"}, "dump_intv: 'log5:0:10': " + schedule_rule(E::LogS));
  said({"stat_intv=log5:10"}, "stat_intv: 'log5:10': " + schedule_rule(E::LogStart));
  said({"dump_intv=10,0.1:1"}, "dump_intv: '10': " + schedule_rule(E::BareNotLast));
  said({"checkpoint_intv=0.1:1,1:0.5"}, "checkpoint_intv: '1:0.5': " + schedule_rule(E::Order));
  // Each key its own problem
  CHECK(problems({"dump_intv=abc", "stat_intv=log5:10"}).size() == 2);
  // Not refused: an interval that does not divide its segment, one below a
  // step, a start off the intervals
  CHECK_NOTHROW(parse_mini({"dump_intv=0.3:1", "stat_intv=0.001:1", "checkpoint_intv=0.1:1", "t0=0.05"}));

  // The messages say what is allowed
  CHECK(schedule_rule(E::Syntax).find("logN:S:U") != std::string::npos);
  CHECK(schedule_rule(E::NonPositive).find("0 alone") != std::string::npos);
  CHECK(schedule_rule(E::LogN).find("1000") != std::string::npos);
}

// --- events ----------------------------------------------------------------------

TEST_CASE("A bare number is every floor(X/dt) steps, at least one, from step 0", "[schedule]"){
  for (const double dt : {0.1, 0.003, 0.05, 1e-7})
    for (const double X : {0.25, 0.3, 0.1, 10., 1e-8, 1e9}){
      Schedule s;
      REQUIRE_FALSE(parse_schedule(partrac::value_to_string(partrac::Value(X)), s));
      EventClock c(s, StepGrid{0., 0, dt}, 0.);
      REQUIRE(c.plain());
      for (long long it = 0; it < 2000; ++it)
        REQUIRE(c.at(it) == at_interval(Uint(it), X, dt));
      // The count every chunk is counted in
      const long long m = static_cast<long long>(steps_per(X, dt));
      CHECK(c.count_before(0) == 0);
      if (m > 1000000) continue;
      CHECK(c.count_before(1) == 1);
      CHECK(c.count_before(3*m) == 3);
      CHECK(c.count_before(3*m + 1) == 4);
    }
  // 0 is off
  Schedule s;
  REQUIRE_FALSE(parse_schedule("0", s));
  EventClock c(s, StepGrid{0., 0, 0.1}, 0.);
  CHECK_FALSE(c.on());
  CHECK_FALSE(c.at(0));
}

TEST_CASE("Uniform segments snap start + n X to the nearest step", "[schedule]"){
  const double dt = 0.003;   // 0.1/dt = 33.3...
  const auto s = steps("0.1:1,1:10,10", dt, 40.);
  // Every 0.1 to 1: 33 or 34 steps; every 1 to 10; every 10; the switches once
  const std::vector<long long> first{0, 33, 67, 100, 133, 167, 200, 233, 267, 300, 333};
  const std::vector<long long> second{333, 667, 1000, 1333, 1667, 2000, 2333, 2667, 3000, 3333};
  const std::vector<long long> third{3333, 6667, 10000, 13333};
  CHECK(s == concat(concat(first, second), third));
  CHECK(s.size() == 11 + 9 + 3);
}

TEST_CASE("X:U is start + n X below U, and U", "[schedule]"){
  const double dt = 0.003;
  // 0.3 does not divide 1: 0, 0.3, 0.6, 0.9, then 1; nothing after
  CHECK(steps("0.3:1", dt, 2.) == snapped({0., 0.3, 0.6, 0.9, 1.}, dt));
  // Below a step: every step to 1, then every 1
  const auto fine = steps("0.001:1,1", dt, 3.);
  CHECK(fine.size() == 334 + 2);
  CHECK(std::vector<long long>(fine.begin(), fine.begin() + 3) == std::vector<long long>{0, 1, 2});
  // U before a whole X: the start and U
  CHECK(steps("1:0.5", dt, 2.) == snapped({0., 0.5}, dt));
  // From t0 = 0.05: 0.05, 0.15, ..., 0.95, then 1
  std::vector<double> t;
  for (int n = 0; n < 10; ++n)
    t.push_back(0.05 + 0.1*n);
  t.push_back(1.);
  CHECK(steps("0.1:1", dt, 2., nullptr, 0, 0.05) == snapped(t, dt));
  // A U on round-off of a multiple is one event with the last start + n X
  CHECK(steps("0.1:0.3,0.1:0.7", dt, 1.) == snapped({0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7}, dt));
}

TEST_CASE("A time half way between two steps is at the later", "[schedule]"){
  // 0.15 a step of 0.1: 1.5, 3, 4.5, ... steps, the halves up
  CHECK(steps("0.15:1.05,0.15", 0.1, 1.5)
        == std::vector<long long>{0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15});
  CHECK(steps("0.3:0.9,0.3", 0.2, 1.6) == std::vector<long long>{0, 2, 3, 5, 6, 8});
  const StepGrid g{0., 0, 0.1};
  CHECK(g.snap(0.05) == 1);
  CHECK(g.snap(0.25) == 3);
  CHECK(g.snap(0.35) == 4);
  CHECK(g.snap(0.0499) == 0);
}

TEST_CASE("Log segments are geometric, snapped, and merged on one step", "[schedule]"){
  std::vector<double> logs;
  for (int k = 0; k <= 30; ++k)
    logs.push_back(0.1*std::pow(10., k/10.));
  {
    // No two on one step at dt = 0.003; 100 is t_30 and the tail's start;
    // t0 = 0 first, then nothing to S
    const double dt = 0.003;
    const auto s = steps("log10:0.1:100,10", dt, 130.);
    const std::vector<long long> head{0, 33, 42, 53, 67, 84, 105, 133, 167, 210, 265, 333};
    CHECK(std::vector<long long>(s.begin(), s.begin() + 12) == head);
    CHECK(s == concat(concat({0}, snapped(logs, dt)), {33333, 36667, 40000, 43333}));
    CHECK(s.size() == 1 + 31 + 3);
  }
  {
    // At dt = 0.05, 0.1259 and 0.1585 are both step 3: one event
    const double dt = 0.05;
    const auto s = steps("log10:0.1:100,10", dt, 130.);
    CHECK(std::vector<long long>(s.begin(), s.begin() + 7) == std::vector<long long>{0, 2, 3, 4, 5, 6, 8});
    CHECK(s.size() == 1 + 30 + 3);
    CHECK(s == concat(concat({0}, snapped(logs, dt)), {2200, 2400, 2600}));
  }
  {
    // Every 0.1 to 1, 5 a decade from 1 to 100, every 10
    const double dt = 0.003;
    const auto s = steps("0.1:1,log5:100,10", dt, 130.);
    const std::vector<long long> first{0, 33, 67, 100, 133, 167, 200, 233, 267, 300, 333};
    const std::vector<long long> second{333, 528, 837, 1327, 2103, 3333, 5283, 8373, 13270, 21032, 33333};
    CHECK(s == concat(concat(first, second), {33333, 36667, 40000, 43333}));
  }
  {
    // Many a decade on few steps: 0.001 to 1 at 10 a decade, dt = 0.003;
    // t0 and 0.001 on step 0, one event
    const double dt = 0.003;
    const auto s = steps("log10:0.001:1", dt, 1.);
    std::vector<double> t;
    for (int k = 0; k <= 30; ++k)
      t.push_back(0.001*std::pow(10., k/10.));
    CHECK(s == snapped(t, dt));
    CHECK(s.size() < 31);
    CHECK(std::vector<long long>(s.begin(), s.begin() + 4) == std::vector<long long>{0, 1, 2, 3});
  }
}

TEST_CASE("A schedule that ends at an until has no event after it", "[schedule]"){
  const double dt = 0.003;
  // Every 0.1 to 1, every 1 to 10, then nothing to T = 100
  const auto s = steps("0.1:1,1:10", dt, 100.);
  CHECK(s.back() == 3333);
  CHECK(s.size() == 20);
  // 5 a decade from 1 to 50: the last is 10^(8/5) = 39.8, not 50; t0 first
  const auto l = steps("log5:1:50", dt, 100.);
  CHECK(l.back() == std::llround(std::pow(10., 8./5.)/dt));
  CHECK(l.front() == 0);
  CHECK(l.size() == 1 + 9);
  // A bare tail goes on from 50
  const auto t = steps("log5:1:50,10", dt, 100.);
  CHECK(std::vector<long long>(t.end() - 6, t.end())
        == std::vector<long long>{16667, 20000, 23333, 26667, 30000, 33333});
  CHECK(t[l.size() - 1] == l.back());
}

TEST_CASE("A restart continues the schedule from any step", "[schedule]"){
  const double dt = 0.003;
  for (const std::string text : {"0.1:1,1:10,10", "log10:0.1:100,10", "0.1:1,log5:100,10"}){
    INFO(text);
    const auto all = steps(text, dt, 130.);
    Schedule s;
    REQUIRE_FALSE(parse_schedule(text, s));
    // Mid uniform, mid log, on an event: the same steps after it, and the
    // count before it
    for (const long long it : {150LL, 1500LL, 5283LL, 6000LL, 20000LL, 35000LL}){
      EventClock c(s, StepGrid{0., 0, dt}, 0.);
      std::vector<long long> after;
      for (long long k = it; k <= 43333; ++k)
        if (c.at(k)) after.push_back(k);
      const auto from = std::lower_bound(all.begin(), all.end(), it);
      CHECK(after == std::vector<long long>(from, all.end()));
      CHECK(c.count_before(it) == from - all.begin());
    }
  }
  // A restart with another step: the same times, snapped to the new steps
  const std::string text = "0.1:1,log5:100,10";
  Schedule s;
  REQUIRE_FALSE(parse_schedule(text, s));
  const StepGrid resumed{20.001, 6667, 0.002};
  EventClock c(s, resumed, 0.);
  std::vector<double> t;
  for (long long k = c.next(6667); k != EventClock::never && resumed.time(k) <= 130.0015; k = c.next(k + 1))
    t.push_back(resumed.time(k));
  const std::vector<double> want{std::pow(10., 7./5.), std::pow(10., 8./5.), std::pow(10., 9./5.),
                                 100., 110., 120., 130.};
  REQUIRE(t.size() == want.size());
  for (std::size_t i = 0; i < t.size(); ++i)
    CHECK(std::abs(t[i] - want[i]) <= 0.001 + 1e-9);
}

TEST_CASE("Every schedule has an event at the run's start, a bare number from step 0", "[schedule]"){
  const double dt = 0.003;
  // t0 = 0, then 20 a decade from 0.1, or from 1
  for (const std::string text : {"log20:0.1:80,10", "log20:1:80,10"}){
    const auto s = steps(text, dt, 100.);
    CHECK(s[0] == 0);
    CHECK(s[1] == std::llround((text[6] == '0' ? 0.1 : 1.)/dt));
  }
  // From t0 = 0.5: the start, then 1, ...; nothing before t0
  const auto s = steps("log20:1:80,10", dt, 100., nullptr, LLONG_MIN/4, 0.5);
  CHECK(std::vector<long long>(s.begin(), s.begin() + 2) == std::vector<long long>{167, 333});
  // The start on an event already: one
  CHECK(steps("log20:0.001:80,10", dt, 1.)[0] == 0);
  CHECK(steps("log20:0.001:80,10", dt, 1.)[1] == 1);
  CHECK(steps("0.1:1,10", dt, 1.)[1] == 33);
  // A restart that moved the grid has no event at its own start, nor at t0
  Schedule sc;
  REQUIRE_FALSE(parse_schedule("log20:1:80,10", sc));
  const StepGrid resumed{0.4, 200, 0.002};
  EventClock c(sc, resumed, 0.);
  CHECK(resumed.time(c.next(200)) == Approx(1.0).margin(0.001));
  CHECK(c.count_before(200) == 1);
  // A bare number is untouched: from t0 on its steps
  REQUIRE_FALSE(parse_schedule("0.25", sc));
  EventClock b(sc, StepGrid{0., 0, 0.1}, 0.);
  CHECK((b.at(0) && !b.at(1) && b.at(2)));
}

TEST_CASE("Dumps are counted from the run's start", "[schedule]"){
  // t0 = 0.5, 10 a decade from 0.1: the times below t0 are not counted (the
  // chunks of dump_chunk_size)
  const double dt = 0.003;
  Schedule s;
  REQUIRE_FALSE(parse_schedule("log10:0.1:2,0.5", s));
  const StepGrid g{0.5, 0, dt};
  EventClock c(s, g, 0.5);
  CHECK(c.count_before(0) == 0);
  CHECK(c.next(0) == 0);
  CHECK(c.count_before(1) == 1);
  // 0.5 (0.501 on its step), 0.631, 0.794
  CHECK(c.count_before(g.snap(1.)) == 3);
  CHECK(c.count_before(g.snap(3.) + 1) == 10);
}

TEST_CASE("A restart with another step passes the times before it", "[schedule]"){
  // dt 0.003 to 0.2 at t = 0.453: 0.4 is nearest the restart's step, but
  // was dumped before it
  Schedule s;
  REQUIRE_FALSE(parse_schedule("0.2:1,10", s));
  const StepGrid resumed{0.453, 151, 0.2};
  EventClock c(s, resumed, 0.);
  CHECK(resumed.snap(0.4) == 151);
  CHECK_FALSE(c.at(151));
  CHECK(c.next(151) == 152);
  // 0, 0.2 and 0.4 came before
  CHECK(c.count_before(151) == 3);
  // A time from the restart on is not passed
  REQUIRE_FALSE(parse_schedule("0.453:1,10", s));
  EventClock d(s, resumed, 0.);
  CHECK(d.at(151));
}

// --- the keys --------------------------------------------------------------------

TEST_CASE("dump_intv, stat_intv and checkpoint_intv each take a schedule", "[schedule]"){
  partrac::Params p = parse_mini({"dump_intv=0.1:1,log5:100,10", "stat_intv=log10:0.1:100,10",
                                  "checkpoint_intv=0.1:1,1:10,10"});
  Run run{p, nullptr, RunFolders{}, {}, 0., 130., false, false, 0.};
  const StepGrid g = step_grid(run, 0.003);
  const auto times_of = [&](const std::string& key){
    std::vector<long long> out;
    EventClock c = event_clock(run, key, g);
    for (long long k = 0; k <= 43333; ++k)
      if (c.at(k)) out.push_back(k);
    return out;
  };
  CHECK(times_of("dump_intv") == steps("0.1:1,log5:100,10", 0.003, 130.));
  CHECK(times_of("stat_intv") == steps("log10:0.1:100,10", 0.003, 130.));
  CHECK(times_of("checkpoint_intv") == steps("0.1:1,1:10,10", 0.003, 130.));

  // On a restart at a time off the steps from t0, the grid starts there
  p.set<std::string>("restart_folder", "somewhere");
  p.set<Uint>("it", 500);
  p.set<double>("t", 1.6);
  const StepGrid r = step_grid(run, 0.003);
  CHECK((r.t_start == 1.6 && r.it_start == 500));
  EventClock c = event_clock(run, "dump_intv", r);
  // 10^(1/5) = 1.585 came before 1.6
  CHECK(r.time(c.next(500)) == Approx(std::pow(10., 2./5.)).margin(0.0015));
}

TEST_CASE("A schedule goes through the parameter file unchanged, a bare number as a number", "[schedule]"){
  TempDir d;
  const std::string sched = "0.1:1,log5:100,10";
  partrac::Params p = parse_mini({"dump_intv=" + sched, "stat_intv=0.001", "checkpoint_intv=log10:0.1:100,10"});
  CHECK(p.get<std::string>("dump_intv") == sched);
  // A bare number floored at the step, written as a number
  CHECK(p.get<std::string>("stat_intv") == partrac::value_to_string(partrac::Value(0.003)));
  p.dump((d.path / "Checkpoints").string());
  std::ifstream f(d.path / "Checkpoints" / "params.dat");
  const std::string text((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
  CHECK(text.find("\ndump_intv=" + sched + "\n") != std::string::npos);
  CHECK(text.find("\ncheckpoint_intv=log10:0.1:100,10\n") != std::string::npos);
  CHECK(text.find("\nstat_intv=0.0030000000000000001\n") != std::string::npos);

  partrac::Params q = mini().parse({"mini", "in.dat", "restart_folder=" + d.path.string()});
  CHECK(q.get<std::string>("dump_intv") == sched);
  CHECK(q.get<std::string>("checkpoint_intv") == "log10:0.1:100,10");
  CHECK(std::stod(q.get<std::string>("stat_intv")) == 0.003);

  // The defaults and a bare 0 as the numbers were
  partrac::Params r = parse_mini({"stat_intv=0"});
  CHECK(r.get<std::string>("dump_intv") == "100");
  CHECK(r.get<std::string>("checkpoint_intv") == "1000");
  CHECK(r.get<std::string>("stat_intv") == "0");
  // The clamp runs once, after the file and the command line
  TempDir e;
  std::ofstream(e.path / "Checkpoints" / "params.dat") << "dump_intv=0.1\nNrw=10\ndt=0.003\n";
  partrac::Params s = mini().parse({"mini", "in.dat", "restart_folder=" + e.path.string(), "dt=1.0"});
  CHECK(s.get<std::string>("dump_intv") == "1");
}
