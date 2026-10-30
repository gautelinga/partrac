// Host side of the OpenFOAM reader: the plugin opened on the first case read,
// from PARTRAC_OPENFOAM_PLUGIN, beside the app or the build tree

#include "openfoam_load.hpp"

#include <cstdlib>
#include <string>

#include "Error.hpp"

#if defined(PARTRAC_OPENFOAM_PLUGIN)
#include <dlfcn.h>
#include <filesystem>
#include <vector>
#endif

namespace openfoam_load {

namespace {

using ReadCase = int (*)(const char*, const char*, CaseData*, std::string*);
using ReadField = int (*)(const CaseData*, const char*, const char*, FieldData*, std::string*);
using HasEmpty = int (*)(const char*, int*, std::string*);
using Layout = std::uint64_t (*)();

struct Entries {
  ReadCase read_case = nullptr;
  ReadField read_field = nullptr;
  HasEmpty has_empty = nullptr;
};

#if defined(PARTRAC_OPENFOAM_PLUGIN)

namespace fs = std::filesystem;

template<typename F>
F symbol(void* h, const char* name, const std::string& path){
  void* s = dlsym(h, name);
  if (!s) partrac::fail(path, ": no ", name, " in the OpenFOAM reader; it is not this build's "
                        "(PARTRAC_OPENFOAM_PLUGIN names the reader to load)");
  return reinterpret_cast<F>(s);
}

// The build's OpenFOAM environment; refuse another
void set_project(const char* name, const std::string& build, const bool is_dir){
  const char* env = std::getenv(name);
  if (env && *env){
    std::error_code ec;
    const bool same = env == build || (is_dir && fs::equivalent(env, build, ec));
    if (!same)
      partrac::fail("the environment's OpenFOAM is not the one this build reads cases with: ", name, " is ", env,
                    ", the build's ", build, "; unset it, or source the build's OpenFOAM");
  }
  setenv(name, build.c_str(), 1);
}

// The plugin: the variable's, else the first that exists of those beside the
// app and the build's
std::string plugin_path(){
  const char* env = std::getenv("PARTRAC_OPENFOAM_PLUGIN");
  if (env && *env) return env;
  std::vector<fs::path> tried;
  std::error_code ec;
  const fs::path exe = fs::read_symlink("/proc/self/exe", ec);
  if (!ec){
    const fs::path name = fs::path(PARTRAC_OPENFOAM_PLUGIN).filename(), bin = exe.parent_path();
    for (const fs::path& d : {bin/PARTRAC_OPENFOAM_LIBDIR_REL, bin/".."/"lib", bin/"lib", bin})
      tried.push_back((d/name).lexically_normal());
  }
  tried.push_back(PARTRAC_OPENFOAM_PLUGIN);
  std::string list;
  for (const fs::path& p : tried){
    if (fs::exists(p, ec)) return p.string();
    list += (list.empty() ? "" : ", ") + p.string();
  }
  partrac::fail("cannot find the OpenFOAM reader (looked for ", list, "); set PARTRAC_OPENFOAM_PLUGIN "
                "to its path");
}

const Entries& entries(){
  static Entries e;
  if (e.read_case) return e;
  // OpenFOAM's static initialisation reads its etc/controlDict
  set_project("WM_PROJECT_DIR", PARTRAC_WM_PROJECT_DIR, true);
  set_project("WM_PROJECT_VERSION", PARTRAC_WM_PROJECT_VERSION, false);
  const std::string path = plugin_path();
  void* h = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
  if (!h) partrac::fail("cannot load the OpenFOAM reader: ", dlerror(), " (PARTRAC_OPENFOAM_PLUGIN "
                        "names the reader to load)");
  const std::uint64_t stamp = symbol<Layout>(h, "partrac_openfoam_layout", path)();
  if (stamp != layout_stamp)
    partrac::fail(path, ": an OpenFOAM reader of another build (layout ", stamp, ", this build's ",
                  layout_stamp, "); rebuild it, or set PARTRAC_OPENFOAM_PLUGIN to this build's");
  e.read_field = symbol<ReadField>(h, "partrac_openfoam_read_field", path);
  e.has_empty = symbol<HasEmpty>(h, "partrac_openfoam_has_empty", path);
  e.read_case = symbol<ReadCase>(h, "partrac_openfoam_read_case", path);
  return e;
}

#else

const Entries& entries(){
  partrac::fail("mode openfoam needs a build with PARTRAC_ENABLE_OPENFOAM=ON.");
}

#endif

}  // namespace

bool available(){
#if defined(PARTRAC_OPENFOAM_PLUGIN)
  return true;
#else
  return false;
#endif
}

CaseData read_case(const std::string& dir, const std::string& velocity_field){
  const Entries& e = entries();
  CaseData c;
  std::string err;
  if (e.read_case(dir.c_str(), velocity_field.c_str(), &c, &err) != 0) partrac::fail(err);
  return c;
}

FieldData read_field(const CaseData& c, const std::string& time, const std::string& name){
  const Entries& e = entries();
  FieldData f;
  std::string err;
  if (e.read_field(&c, time.c_str(), name.c_str(), &f, &err) != 0) partrac::fail(err);
  return f;
}

bool has_empty_patches(const std::string& dir){
  const Entries& e = entries();
  int out = 0;
  std::string err;
  if (e.has_empty(dir.c_str(), &out, &err) != 0) partrac::fail(err);
  return out != 0;
}

}  // namespace openfoam_load
