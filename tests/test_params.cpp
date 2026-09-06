#include <catch2/catch.hpp>

#include <chrono>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "Params.hpp"

using partrac::ParamError;
using partrac::Params;
using partrac::Schema;
using partrac::Source;

namespace {

using Uint = std::size_t;

// args[0] is the program name, args[1] the positional interpolator file, as in
// the real invocation `partrac <file> key=value ...`.
std::vector<std::string> cli(std::vector<std::string> pairs) {
  std::vector<std::string> a{"app", "mesh.dat"};
  for (auto& p : pairs) a.push_back(std::move(p));
  return a;
}

Schema basic() {
  Schema s("testapp");
  s.opt<double>("Dm", 1.0, "molecular diffusivity");
  s.opt<double>("dt", 1.0, "timestep");
  s.require<Uint>("Nrw", "number of particles");
  s.opt<bool>("verbose", false, "chatty output");
  s.opt<std::string>("init_mode", "uniform_x", "initial distribution");
  return s;
}

// Collects every problem message into one blob for substring assertions.
std::string all_problems(const ParamError& e) {
  std::string s;
  for (const auto& p : e.problems()) s += p + "\n";
  return s;
}

struct TempDir {
  std::filesystem::path path;
  TempDir() {
    static int counter = 0;
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    path = std::filesystem::temp_directory_path() /
           ("partrac_test_" + std::to_string(stamp) + "_" +
            std::to_string(counter++));
    std::filesystem::create_directories(path / "Checkpoints");
  }
  ~TempDir() {
    std::error_code ec;
    std::filesystem::remove_all(path, ec);
  }
  void write_checkpoint(const std::string& contents) const {
    std::ofstream f(path / "Checkpoints" / "params.dat");
    f << contents;
  }
};

}  // namespace

// ---------------------------------------------------------------------------
// Defaults and required
// ---------------------------------------------------------------------------

TEST_CASE("defaults are applied when a parameter is absent", "[params]") {
  Params p = basic().parse(cli({"Nrw=100"}));
  REQUIRE(p.get<double>("Dm") == 1.0);
  REQUIRE(p.source("Dm") == Source::Default);
  REQUIRE(p.was_set("Dm") == false);
  REQUIRE(p.get<Uint>("Nrw") == 100);
  REQUIRE(p.source("Nrw") == Source::Cmdline);
  REQUIRE(p.was_set("Nrw") == true);
}

TEST_CASE("a missing required parameter is an error naming the key", "[params]") {
  REQUIRE_THROWS_AS(basic().parse(cli({})), ParamError);
  try {
    basic().parse(cli({}));
  } catch (const ParamError& e) {
    REQUIRE(all_problems(e).find("Nrw") != std::string::npos);
    REQUIRE(all_problems(e).find("number of particles") != std::string::npos);
  }
}

TEST_CASE("reading a parameter that is not declared throws", "[params]") {
  Params p = basic().parse(cli({"Nrw=100"}));
  REQUIRE_THROWS_AS(p.get<double>("nonexistent"), ParamError);
}

// ---------------------------------------------------------------------------
// Unknown keys -- the bug this whole change exists to fix
// ---------------------------------------------------------------------------

TEST_CASE("an unknown parameter is rejected, not silently ignored", "[params]") {
  try {
    basic().parse(cli({"Nrw=100", "bogus=3"}));
    FAIL("expected ParamError");
  } catch (const ParamError& e) {
    REQUIRE(all_problems(e).find("unknown parameter 'bogus'") != std::string::npos);
  }
}

TEST_CASE("a near-miss key suggests the intended one", "[params]") {
  Schema s = basic();
  s.opt<double>("refine_intv", 100.0, "refinement interval");
  try {
    s.parse(cli({"Nrw=100", "refine_int=10"}));
    FAIL("expected ParamError");
  } catch (const ParamError& e) {
    REQUIRE(all_problems(e).find("did you mean 'refine_intv'") != std::string::npos);
  }
}

// A one-character edit is always offered; the limit only rules out the noisy
// matches, so 'nx' gets nothing here but would still suggest 'Lx' if the
// schema had it.
TEST_CASE("a transposed key is suggested, a distant one is not", "[params]") {
  try {
    basic().parse(cli({"Nwr=100"}));
    FAIL("expected ParamError");
  } catch (const ParamError& e) {
    REQUIRE(all_problems(e).find("did you mean 'Nrw'") != std::string::npos);
  }
  try {
    basic().parse(cli({"Nrw=1", "nx=0"}));
    FAIL("expected ParamError");
  } catch (const ParamError& e) {
    REQUIRE(all_problems(e).find("did you mean") == std::string::npos);
  }
}

TEST_CASE("every problem is reported at once, not just the first", "[params]") {
  try {
    basic().parse(cli({"bogus1=1", "bogus2=2", "Dm=notanumber"}));
    FAIL("expected ParamError");
  } catch (const ParamError& e) {
    const std::string all = all_problems(e);
    REQUIRE(all.find("bogus1") != std::string::npos);
    REQUIRE(all.find("bogus2") != std::string::npos);
    REQUIRE(all.find("Dm") != std::string::npos);
    REQUIRE(all.find("Nrw") != std::string::npos);  // still reports the missing one
    REQUIRE(e.problems().size() >= 4);
  }
}

TEST_CASE("a bare argument is rejected", "[params]") {
  // parse_cmd used to skip anything without '=', so --help and typos like
  // "Nrw 100" were silently ignored.
  REQUIRE_THROWS_AS(basic().parse(cli({"Nrw=100", "--nonsense"})), ParamError);
  REQUIRE_THROWS_AS(basic().parse(cli({"Nrw=100", "stray"})), ParamError);
}

TEST_CASE("a duplicated key is an error", "[params]") {
  REQUIRE_THROWS_AS(basic().parse(cli({"Nrw=100", "Nrw=200"})), ParamError);
}

// ---------------------------------------------------------------------------
// Type conversion
// ---------------------------------------------------------------------------

TEST_CASE("trailing garbage in a number is rejected", "[params]") {
  // stodouble used `ss >> d` without checking eof(), so "10meters" parsed as 10.
  REQUIRE_THROWS_AS(basic().parse(cli({"Nrw=100", "Dm=10meters"})), ParamError);
  REQUIRE_THROWS_AS(basic().parse(cli({"Nrw=100", "Dm=abc"})), ParamError);
  REQUIRE_THROWS_AS(basic().parse(cli({"Nrw=100", "Dm="})), ParamError);
}

TEST_CASE("scientific notation still works for integers", "[params]") {
  // commands.txt relies on Nrw=1e6 and Nrw_max=1e7.
  Params p = basic().parse(cli({"Nrw=1e6"}));
  REQUIRE(p.get<Uint>("Nrw") == 1000000u);
}

TEST_CASE("a non-integral value for an integer parameter is rejected", "[params]") {
  REQUIRE_THROWS_AS(basic().parse(cli({"Nrw=3.5"})), ParamError);
  REQUIRE_THROWS_AS(basic().parse(cli({"Nrw=1e20"})), ParamError);
}

TEST_CASE("a negative value for an unsigned parameter is rejected", "[params]") {
  // Nrw_max = -1 on a size_t wrapped to 2^64-1 and was then handed to
  // ParticleSet::resize(), which is a guaranteed bad_alloc.
  Schema s("testapp");
  s.require<Uint>("Nrw_max", "maximum particles");
  Params p = s.parse(cli({"Nrw_max=-1"}));
  REQUIRE_THROWS_AS(p.get<Uint>("Nrw_max"), ParamError);
}

TEST_CASE("large integers survive without truncation", "[params]") {
  // n_accepted/n_declined are long int but were parsed via stoint -> int.
  Schema s("testapp");
  s.opt<long long>("n_accepted", 0, "accepted steps");
  Params p = s.parse(cli({"n_accepted=3000000000"}));
  REQUIRE(p.get<long long>("n_accepted") == 3000000000LL);
  REQUIRE_THROWS_AS(p.get<int>("n_accepted"), ParamError);  // does not fit
}

TEST_CASE("booleans accept the usual spellings and reject nonsense", "[params]") {
  for (const std::string v : {"true", "True", "TRUE", "1", "yes", "on"})
    REQUIRE(basic().parse(cli({"Nrw=1", "verbose=" + v})).get<bool>("verbose"));
  for (const std::string v : {"false", "False", "FALSE", "0", "no", "off"})
    REQUIRE(!basic().parse(cli({"Nrw=1", "verbose=" + v})).get<bool>("verbose"));
  // stobool() used to return false for anything it did not recognise.
  REQUIRE_THROWS_AS(basic().parse(cli({"Nrw=1", "verbose=ture"})), ParamError);
  REQUIRE_THROWS_AS(basic().parse(cli({"Nrw=1", "verbose=maybe"})), ParamError);
}

TEST_CASE("asking for the wrong type throws", "[params]") {
  Params p = basic().parse(cli({"Nrw=100"}));
  REQUIRE_THROWS_AS(p.get<double>("init_mode"), ParamError);
  REQUIRE_THROWS_AS(p.get<std::string>("Dm"), ParamError);
}

// ---------------------------------------------------------------------------
// Precedence, restart files, and the clamps
// ---------------------------------------------------------------------------

namespace {
Schema restartable() {
  Schema s = basic();
  s.opt<std::string>("restart_folder", "", "checkpoint to resume from");
  s.opt<double>("dump_intv", 100.0, "dump interval");
  s.finalize([](Params& p) {
    p.set<double>("dump_intv",
                  std::max(p.get<double>("dump_intv"), p.get<double>("dt")));
  });
  return s;
}
}  // namespace

TEST_CASE("the command line beats the checkpoint file", "[params]") {
  TempDir d;
  d.write_checkpoint("Dm=1\ndt=0.5\nNrw=50\n");
  Params p = restartable().parse(
      cli({"restart_folder=" + d.path.string(), "Dm=2"}));
  REQUIRE(p.get<double>("Dm") == 2.0);
  REQUIRE(p.source("Dm") == Source::Cmdline);
  REQUIRE(p.get<double>("dt") == 0.5);        // only the file set this
  REQUIRE(p.source("dt") == Source::File);
  REQUIRE(p.get<Uint>("Nrw") == 50);          // required, satisfied by the file
}

TEST_CASE("the clamps run once, after the file and the command line", "[params]") {
  // parse_file never applied these, which is exactly why the old code had to
  // re-run parse_cmd as a third pass.
  TempDir d;
  d.write_checkpoint("dump_intv=0.1\nNrw=10\n");
  Params p = restartable().parse(
      cli({"restart_folder=" + d.path.string(), "dt=1.0"}));
  REQUIRE(p.get<double>("dump_intv") == 1.0);
}

TEST_CASE("an unknown key in a checkpoint file is an error by default", "[params]") {
  TempDir d;
  d.write_checkpoint("Nrw=10\nwrite_mode=hdf5\n");
  REQUIRE_THROWS_AS(
      restartable().parse(cli({"restart_folder=" + d.path.string()})),
      ParamError);
}

TEST_CASE("strict_file(false) downgrades unknown file keys to a warning", "[params]") {
  TempDir d;
  d.write_checkpoint("Nrw=10\nwrite_mode=hdf5\n");
  Schema s = restartable();
  s.strict_file(false);
  Params p = s.parse(cli({"restart_folder=" + d.path.string()}));
  REQUIRE(p.get<Uint>("Nrw") == 10);
}

TEST_CASE("a missing restart file is reported", "[params]") {
  REQUIRE_THROWS_AS(restartable().parse(cli({"restart_folder=/no/such/place"})),
                    ParamError);
}

// ---------------------------------------------------------------------------
// Conditional requirements
// ---------------------------------------------------------------------------

namespace {
bool needs_La(const Params& p) {
  const std::string m = p.get<std::string>("init_mode");
  return m.rfind("ellipsoid", 0) == 0 || m.rfind("sheet", 0) == 0 ||
         m.rfind("strip", 0) == 0;
}

Schema conditional() {
  Schema s = basic();
  s.require_if<double>("La", needs_La, "init_mode is a strip, sheet or ellipsoid",
                       "principal extent of the initial condition");
  return s;
}
}  // namespace

TEST_CASE("a conditionally required parameter is enforced", "[params]") {
  // La defaulting to 0.0 divides by zero in EllipsoidInitializer and silently
  // collapses every particle onto one point for strip_*.
  try {
    conditional().parse(cli({"Nrw=100", "init_mode=ellipsoid_xy"}));
    FAIL("expected ParamError");
  } catch (const ParamError& e) {
    const std::string all = all_problems(e);
    REQUIRE(all.find("La") != std::string::npos);
    REQUIRE(all.find("init_mode is a strip, sheet or ellipsoid") != std::string::npos);
  }
  REQUIRE_NOTHROW(conditional().parse(
      cli({"Nrw=100", "init_mode=ellipsoid_xy", "La=2.0"})));
}

TEST_CASE("a conditional parameter is simply unset when not required", "[params]") {
  Params p = conditional().parse(cli({"Nrw=100", "init_mode=uniform_x"}));
  REQUIRE(p.has("La") == false);
  REQUIRE_THROWS_AS(p.get<double>("La"), ParamError);
  REQUIRE(p.get_or<double>("La", -1.0) == -1.0);
}

TEST_CASE("a conditional parameter with a default must still be given", "[params]") {
  // partrac reads inject_intv unconditionally to precompute the interval, so it
  // needs a value even when inject is off -- but a default must not satisfy the
  // condition, or inject=true with no interval is a modulo by zero.
  Schema s = basic();
  s.opt<bool>("inject", false, "inject new particles");
  s.require_if<double>("inject_intv", 0.0,
                       [](const Params& p) { return p.get<bool>("inject"); },
                       "inject is set", "injection interval");

  Params p = s.parse(cli({"Nrw=1"}));
  REQUIRE(p.get<double>("inject_intv") == 0.0);   // readable when not required
  REQUIRE_THROWS_AS(s.parse(cli({"Nrw=1", "inject=true"})), ParamError);
  REQUIRE_NOTHROW(s.parse(cli({"Nrw=1", "inject=true", "inject_intv=0.5"})));
}

TEST_CASE("value checks and choices are enforced", "[params]") {
  Schema s = basic();
  s.opt<std::string>("exit_plane", "none", "exit plane");
  s.opt<double>("Ln", 0.0, "exit plane position");
  s.opt<double>("filter_intv", 0.0, "filter interval");
  s.choices("exit_plane", {"none", "x", "y", "z"});
  s.check(
      [](const Params& p) {
        const std::string e = p.get<std::string>("exit_plane");
        if (e != "x" && e != "y" && e != "z") return true;
        return p.get<double>("Ln") > 0.0 && p.get<double>("filter_intv") > 0.0;
      },
      "exit_plane requires Ln > 0 and filter_intv > 0");

  // Currently an integer modulo-by-zero (SIGFPE) at partrac.cpp:357.
  REQUIRE_THROWS_AS(s.parse(cli({"Nrw=1", "exit_plane=x"})), ParamError);
  REQUIRE_THROWS_AS(s.parse(cli({"Nrw=1", "exit_plane=w"})), ParamError);
  REQUIRE_NOTHROW(
      s.parse(cli({"Nrw=1", "exit_plane=x", "Ln=1.0", "filter_intv=0.1"})));
}

// ---------------------------------------------------------------------------
// Runtime (program-computed) parameters
// ---------------------------------------------------------------------------

TEST_CASE("runtime parameters are rejected on the command line", "[params]") {
  Schema s = basic();
  s.runtime<double>("Lx", 0.0, "domain size, from the interpolator");
  try {
    s.parse(cli({"Nrw=1", "Lx=200"}));
    FAIL("expected ParamError");
  } catch (const ParamError& e) {
    REQUIRE(all_problems(e).find("computed by the program") != std::string::npos);
  }
}

TEST_CASE("runtime parameters round-trip through a dump", "[params]") {
  // Lx/Ly/Lz used to be written by write_params_to_file but had no set_param
  // entry, so the code could not read back its own output.
  TempDir d;
  Schema s = basic();
  s.opt<std::string>("restart_folder", "", "checkpoint to resume from");
  s.runtime<double>("Lx", 0.0, "domain size");

  Params p = s.parse(cli({"Nrw=7"}));
  p.set<double>("Lx", 2133.33);
  REQUIRE(p.source("Lx") == Source::Runtime);
  p.dump((d.path / "Checkpoints").string());

  Params q = s.parse(cli({"restart_folder=" + d.path.string()}));
  REQUIRE(q.get<double>("Lx") == Approx(2133.33));
  REQUIRE(q.get<Uint>("Nrw") == 7);
}

TEST_CASE("setting an undeclared or mistyped parameter throws", "[params]") {
  Params p = basic().parse(cli({"Nrw=1"}));
  REQUIRE_THROWS_AS(p.set<double>("bogus", 1.0), ParamError);
  REQUIRE_THROWS_AS(p.set<std::string>("Dm", "x"), ParamError);
  REQUIRE_NOTHROW(p.set<double>("Dm", 2.0));
  REQUIRE(p.get<double>("Dm") == 2.0);
}

// ---------------------------------------------------------------------------
// Round-trip fidelity
// ---------------------------------------------------------------------------

TEST_CASE("doubles survive a dump/parse round trip exactly", "[params]") {
  // write_param used the default 6 significant digits, so a checkpoint written
  // at t = 1234567.89 came back as 1.23457e+06.
  TempDir d;
  Schema s("testapp");
  s.opt<std::string>("restart_folder", "", "checkpoint to resume from");
  s.opt<double>("t", 0.0, "current time");
  s.opt<double>("Dm", 0.0, "diffusivity");

  const double t = 1234567.89;
  const double Dm = 1.0 / 3.0;
  Params p = s.parse(cli({}));
  p.set<double>("t", t);
  p.set<double>("Dm", Dm);
  p.dump((d.path / "Checkpoints").string());

  Params q = s.parse(cli({"restart_folder=" + d.path.string()}));
  REQUIRE(q.get<double>("t") == t);    // bit-identical, not approximately
  REQUIRE(q.get<double>("Dm") == Dm);
}

// ---------------------------------------------------------------------------
// Help and schema self-validation
// ---------------------------------------------------------------------------

TEST_CASE("--help short-circuits validation and lists every key", "[params]") {
  Params p = conditional().parse(cli({"--help"}));
  REQUIRE(p.help_requested());

  const std::string h = conditional().help();
  for (const std::string k : {"Dm", "dt", "Nrw", "verbose", "init_mode", "La"})
    REQUIRE(h.find(k) != std::string::npos);
  REQUIRE(h.find("Required:") != std::string::npos);
}

TEST_CASE("--check is recognised as a dry run", "[params]") {
  Params p = basic().parse(cli({"Nrw=1", "--check"}));
  REQUIRE(p.check_only());
  REQUIRE(!basic().parse(cli({"Nrw=1"})).check_only());
}

TEST_CASE("a well-formed schema validates", "[params]") {
  REQUIRE_NOTHROW(basic().validate_self());
  REQUIRE_NOTHROW(conditional().validate_self());
  REQUIRE_NOTHROW(restartable().validate_self());
}

TEST_CASE("registering the same key twice is caught", "[params]") {
  Schema s("testapp");
  s.opt<double>("Dm", 1.0, "diffusivity");
  REQUIRE_THROWS_AS(s.opt<double>("Dm", 2.0, "again"), ParamError);
}

TEST_CASE("a predicate reading a conditional key is caught by validate_self",
          "[params]") {
  // This is the trap with predicates: 'La' may be unset, so a condition that
  // reads it would throw mid-parse.  Catch it at test time instead.
  Schema s = basic();
  s.require_if<double>("La", needs_La, "init_mode needs it", "extent");
  s.require_if<double>("Lb",
                       [](const Params& p) { return p.get<double>("La") > 1.0; },
                       "La is large", "second extent");
  REQUIRE_THROWS_AS(s.validate_self(), ParamError);
}
