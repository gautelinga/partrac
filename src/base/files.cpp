#include <cstdlib>
#include <filesystem>
#include <iostream>
#include "Error.hpp"
#include "files.hpp"

std::string create_folder(const std::string& folder){
  if (!std::filesystem::is_directory(folder)){
    std::filesystem::create_directory(folder);
  }
  return folder;
}

void verify_file_exists(const std::string& infilename){
  if (!std::filesystem::exists(infilename)){
    partrac::fail("no such file: ", infilename);
  }
}
