#include <catch2/catch.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>

#include "Params.hpp"
#include "AnalyticInterpol.hpp"
#include "loader_params.hpp"

// Every parameter file of the examples parses against the schema of the loader
// that reads it, so the key lists in loader_params.hpp and in the expressions'
// add_params cover what the examples use; a key they miss would stop those
// cases at startup. An OpenFOAM case's partrac_params.dat is mode openfoam's.
// A plain dolfin file does not say which of tet, triangle and fenics reads it,
// so it is checked against tet's keys; what a mode of its own declares is not
// caught here, neither a tet-only key in a file a triangle run reads nor
// renumber_cells, which only fenics takes. Directories in PARTRAC_PARAM_DIRS
// (colon separated) are scanned too, for case collections kept outside the
// repository.

namespace {

namespace fs = std::filesystem;

// The schema of the loader that reads the file
partrac::Schema schema_for(const fs::path& file) {
  const std::string name = file.filename().string();
  const std::string path = file.string();
  if (name == "felbm_params.dat")
    return felbm_schema();
  if (name == "partrac_params.dat")
    return openfoam_schema();
  if (name == "expr_params.dat") {
    const ExprKind* kind = find_expr_kind(partrac::peek_file(path, "expression"));
    if (!kind)
      throw partrac::ParamError(path, {"unknown expression"});
    return expression_schema(*kind);
  }
  // the mode is not in the file; the two freq modes declare the same keys
  if (!partrac::peek_file(path, "freqstamps").empty())
    return simplex_freq_schema("trianglefreq");
  const std::string u = partrac::peek_file(path, "u");
  if (u.size() > 5 && u.substr(u.size() - 5) == ".xdmf")
    return xdmf_schema("xdmftet");
  // the mode is not in the file; tet is the one the examples are written for
  return dolfin_h5_schema("tet");
}

std::vector<fs::path> param_files() {
  std::vector<std::string> roots{std::string(PARTRAC_SOURCE_DIR) + "/data_example"};
  if (const char* extra = std::getenv("PARTRAC_PARAM_DIRS")) {
    std::string s(extra);
    for (std::size_t a = 0, b; a <= s.size(); a = b + 1) {
      b = s.find(':', a);
      if (b == std::string::npos) b = s.size();
      if (b > a) roots.push_back(s.substr(a, b - a));
    }
  }
  std::vector<fs::path> files;
  for (const auto& root : roots) {
    for (const auto& e : fs::recursive_directory_iterator(root, fs::directory_options::skip_permission_denied)) {
      const std::string n = e.path().filename().string();
      if (e.is_regular_file() && (n == "dolfin_params.dat" || n == "expr_params.dat" || n == "felbm_params.dat"
                                  || n == "partrac_params.dat"))
        files.push_back(e.path());
    }
  }
  return files;
}

}  // namespace

TEST_CASE("every example parameter file parses against its loader's keys", "[params][interpol]") {
  const auto files = param_files();
  REQUIRE(!files.empty());
  for (const auto& f : files) {
    INFO(f.string());
    REQUIRE_NOTHROW(schema_for(f).parse_file(f.string()));
  }
}

TEST_CASE("mode openfoam refuses the keys the case itself says", "[params][interpol]") {
  // The periodic axes are the case's cyclic patches; a periodic_* given
  // could only agree with them or contradict them
  // A file of this process's own, so concurrent runs do not share it
  const auto parse = [](const std::string& text){
    const std::string path = (std::filesystem::temp_directory_path()
                              / ("partrac_openfoam_params_" + std::to_string(::getpid()) + ".dat")).string();
    { std::ofstream f(path); f << text; }
    const auto p = openfoam_schema().parse_file(path);
    std::filesystem::remove(path);
    return p;
  };
  REQUIRE_NOTHROW(parse("velocity_field=U\nsplit=6\n"));
  REQUIRE(parse("").get<std::string>("split") == "12");
  REQUIRE(parse("").get<bool>("wall_upgrade"));
  REQUIRE_THROWS_AS(parse("periodic_x=true\n"), partrac::ParamError);
  REQUIRE_THROWS_AS(parse("periodic_z=false\n"), partrac::ParamError);
  REQUIRE_THROWS_AS(parse("split=5\n"), partrac::ParamError);
  REQUIRE_THROWS_AS(parse("mesh=mesh.h5\n"), partrac::ParamError);
  // The phase field is named, never switched on by include_phi
  REQUIRE(parse("phase_field=alpha.water\n").get<std::string>("phase_field") == "alpha.water");
  REQUIRE_FALSE(parse("").has("phase_field"));
  REQUIRE_THROWS_AS(parse("include_phi=true\n"), partrac::ParamError);
}

TEST_CASE("wall_p2 is refused by name, pointing at wall_upgrade", "[params][interpol]") {
  // The key is gone without an alias: a file that still holds it must say
  // what replaced it, in every schema that had it
  const std::string path = (std::filesystem::temp_directory_path()
                            / ("partrac_wall_p2_" + std::to_string(::getpid()) + ".dat")).string();
  const auto message = [&](const partrac::Schema& s, const std::string& text){
    { std::ofstream f(path); f << text; }
    std::string all;
    try { s.parse_file(path); }
    catch (const partrac::ParamError& e){ for (const auto& p : e.problems()) all += p + "\n"; }
    std::filesystem::remove(path);
    return all;
  };
  const std::string h5 = "velocity_space=P2\npressure_space=P1\ntimestamps=t.dat\nmesh=m.h5\n";
  const std::string xdmf = "u=u.xdmf\np=p.xdmf\n";
  for (const auto& [schema, base] : std::vector<std::pair<partrac::Schema, std::string>>{
         {dolfin_h5_schema("tet"), h5}, {dolfin_h5_schema("triangle"), h5},
         {xdmf_schema("xdmftet"), xdmf}, {openfoam_schema(), ""}}){
    for (const std::string v : {"edge", "none"}){
      const std::string m = message(schema, base + "wall_p2=" + v + "\n");
      REQUIRE(m.find("wall_p2") != std::string::npos);
      REQUIRE(m.find("wall_upgrade") != std::string::npos);
    }
    REQUIRE(message(schema, base + "wall_upgrade=true\n").empty());
    schema.validate_self();
  }
  // The defaults: true where the velocity is P1 by construction
  { std::ofstream f(path); f << h5; }
  REQUIRE_FALSE(dolfin_h5_schema("tet").parse_file(path).get<bool>("wall_upgrade"));
  REQUIRE_FALSE(dolfin_h5_schema("triangle").parse_file(path).get<bool>("wall_upgrade"));
  { std::ofstream f(path); f << xdmf; }
  REQUIRE(xdmf_schema("xdmftriangle").parse_file(path).get<bool>("wall_upgrade"));
  REQUIRE(xdmf_schema("xdmftet").parse_file(path).get<bool>("wall_upgrade"));
  std::filesystem::remove(path);
}

TEST_CASE("A flow's time interval may be one instant but not run backwards", "[params][interpol]") {
  // t_min and t_max bound the times a run can start at; t_max before t_min
  // leaves none, and is refused when the file is read
  const std::string path = (std::filesystem::temp_directory_path()
                            / ("partrac_t_interval_" + std::to_string(::getpid()) + ".dat")).string();
  const auto parses = [&](const partrac::Schema& s, const std::string& text){
    { std::ofstream f(path); f << text; }
    bool ok = true;
    try { s.parse_file(path); }
    catch (const partrac::ParamError& e){
      ok = false;
      bool named = false;
      for (const auto& p : e.problems()) named = named || p.find("t_max is before t_min") != std::string::npos;
      REQUIRE(named);
    }
    std::filesystem::remove(path);
    return ok;
  };
  const ExprKind* sine = find_expr_kind("sine_flow");
  REQUIRE(sine);
  std::string expr;
  {
    std::ifstream in(std::string(PARTRAC_SOURCE_DIR) + "/data_example/sine_flow/expr_params.dat");
    for (std::string line; std::getline(in, line); )
      if (line.rfind("t_min=", 0) != 0 && line.rfind("t_max=", 0) != 0) expr += line + "\n";
  }
  const std::string freq = "freqstamps=f.dat\nmesh=m.h5\nvelocity_space=P1\npressure_space=P1\n";
  for (const auto& [schema, base] : std::vector<std::pair<partrac::Schema, std::string>>{
         {expression_schema(*sine), expr}, {simplex_freq_schema("trianglefreq"), freq}}){
    REQUIRE(parses(schema, base + "t_min=0\nt_max=2\n"));
    REQUIRE(parses(schema, base + "t_min=1.5\nt_max=1.5\n"));
    REQUIRE_FALSE(parses(schema, base + "t_min=2\nt_max=1\n"));
    schema.validate_self();
  }
}
