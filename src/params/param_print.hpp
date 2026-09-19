#ifndef __PARAM_PRINT_HPP
#define __PARAM_PRINT_HPP

#include <fstream>
#include <string>

void print_param(const std::string& key, const double val);
void print_param(const std::string& key, const int val);
void print_param(const std::string& key, const std::string& val);
void write_param(std::ofstream &ofile, const std::string& key, double val);
void write_param(std::ofstream &ofile, const std::string& key, int val);
void write_param(std::ofstream &ofile, const std::string& key, long int val);
void write_param(std::ofstream &ofile, const std::string& key, const std::string& val);

#endif
