#ifndef __IO_HPP
#define __IO_HPP

#include <iostream>
#include <vector>
#include <map>
#include <filesystem>

#include "H5Cpp.h"


//using namespace std;
using namespace H5;

std::map<std::string, std::string> load_settings(std::string casefile);
std::vector<std::string> get_files(const std::string& s);
std::vector<std::vector<std::string>> load_grid(std::string infile);
std::vector<std::vector<std::string>> load_fields(std::string infile);
void print_grid(bool** grid, const int nx, const int ny);
void copy_arr(double*** f,
	      double*** f_prev,
	      bool** grid,
	      const int nx,
	      const int ny,
	      const int nc);
void copy_arr(double** a_from,
	      double** a_to,
	      bool** grid,
	      const int nx,
	      const int ny);
void dump2file(std::string filename,
	       double** rho,
	       double** m_x,
	       double** m_y,
	       bool** grid,
	       const int nx,
	       const int ny);
std::string create_folder(const std::string folder);
void verify_file_exists(const std::string infilename);
void print_param(const std::string key, const double val);
void print_param(const std::string key, const int val);
void print_param(const std::string key, const std::string val);
void write_param(std::ofstream &ofile, std::string key, double val);
void write_param(std::ofstream &ofile, std::string key, int val);
void write_param(std::ofstream &ofile, std::string key, long int val);
void write_param(std::ofstream &ofile, std::string key, std::string val);
void load_field(H5File &h5file,
                double*** u,
                const std::string field,
                const int nx,
                const int ny,
                const int nz);
void load_int_field(H5File &h5file,
                    int*** u,
                    const std::string field,
                    const int nx,
                    const int ny,
                    const int nz);
void load_h5(const std::string h5filename,
             double*** ux,
             double*** uy,
             double*** uz,
             double*** rho,
             double*** p,
             const int nx,
             const int ny,
             const int nz,
             const bool verbose);

#endif
