#include <iostream>
#include "param_print.hpp"

// template this...
void print_param(const std::string& key, const double val){
  std::cout << key << " = " << val << std::endl;
}

void print_param(const std::string& key, const int val){
  std::cout << key << " = " << val << std::endl;
}

void print_param(const std::string& key, const std::string& val){
  std::cout << key << " = " << val << std::endl;
}

void write_param(std::ofstream &ofile, const std::string& key, double val){
  ofile << key << "=" << val << std::endl;
}

void write_param(std::ofstream &ofile, const std::string& key, int val){
  ofile << key << "=" << val << std::endl;
}

void write_param(std::ofstream &ofile, const std::string& key, long int val){
  ofile << key << "=" << val << std::endl;
}

void write_param(std::ofstream &ofile, const std::string& key, const std::string& val){
  ofile << key << "=" << val << std::endl;
}
