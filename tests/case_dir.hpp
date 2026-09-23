// Temporary files and directories of the unit tests. A name carries the
// process id and a count, so concurrent runs and cases never share one.
#pragma once

#include <atomic>
#include <filesystem>
#include <string>
#include <system_error>

#include <unistd.h>

// A fresh path under the temporary directory, ending in tag
inline std::filesystem::path temp_path(const std::string& tag){
  static std::atomic<unsigned> count{0};
  return std::filesystem::temp_directory_path()
    / ("partrac_" + std::to_string(::getpid()) + "_" + std::to_string(count++) + "_" + tag);
}

// A case directory of its own per test, removed with its contents
struct CaseDir {
  std::filesystem::path path;
  explicit CaseDir(const std::string& tag) : path(temp_path(tag)) {
    std::filesystem::remove_all(path);
    std::filesystem::create_directories(path);
  }
  ~CaseDir(){ std::error_code ec; std::filesystem::remove_all(path, ec); }
  CaseDir(const CaseDir&) = delete;
  CaseDir& operator=(const CaseDir&) = delete;
  std::string file(const std::string& name) const { return (path / name).string(); }
  std::string params() const { return file("dolfin_params.dat"); }
};
