#include <cstdlib>
#include <filesystem>
#include <iostream>
#include "Error.hpp"
#include "files.hpp"

std::string create_folder(const std::string& folder){
  std::error_code ec;
  std::filesystem::create_directory(folder, ec);   // no error when it exists
  if (ec || !std::filesystem::is_directory(folder)){
    partrac::fail("cannot make the folder ", folder, ec ? ": " + ec.message() : ": a file is in the way");
  }
  return folder;
}

void verify_file_exists(const std::string& infilename){
  if (!std::filesystem::exists(infilename)){
    partrac::fail("no such file: ", infilename);
  }
}
