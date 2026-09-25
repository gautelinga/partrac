#include <catch2/catch.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
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
  REQUIRE(parse("").get<std::string>("wall_p2") == "edge");
  REQUIRE_THROWS_AS(parse("periodic_x=true\n"), partrac::ParamError);
  REQUIRE_THROWS_AS(parse("periodic_z=false\n"), partrac::ParamError);
  REQUIRE_THROWS_AS(parse("split=5\n"), partrac::ParamError);
  REQUIRE_THROWS_AS(parse("mesh=mesh.h5\n"), partrac::ParamError);
  // The phase field is named, never switched on by include_phi
  REQUIRE(parse("phase_field=alpha.water\n").get<std::string>("phase_field") == "alpha.water");
  REQUIRE_FALSE(parse("").has("phase_field"));
  REQUIRE_THROWS_AS(parse("include_phi=true\n"), partrac::ParamError);
}
