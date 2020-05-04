#include <iostream>
#include <vector>
#include <map>
#include "H5Cpp.h"

#ifndef __IO_HPP
#define __IO_HPP

using namespace std;
using namespace H5;

map<string, string> load_settings(string casefile);
vector<string> get_files(const string& s);
vector<vector<string>> load_grid(string infile);
vector<vector<string>> load_fields(string infile);
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
void dump2file(string filename,
	       double** rho,
	       double** m_x,
	       double** m_y,
	       bool** grid,
	       const int nx,
	       const int ny);
string create_folder(const string folder);
void verify_file_exists(const string infilename);
void print_param(const string key, const double val);
void print_param(const string key, const int val);
void print_param(const string key, const string val);
void write_param(ofstream &ofile, string key, double val);
void write_param(ofstream &ofile, string key, int val);
void write_param(ofstream &ofile, string key, long int val);
void write_param(ofstream &ofile, string key, string val);
void load_field(H5File &h5file,
		double*** u,
		const string field,
		const int nx,
		const int ny,
		const int nz);
void load_int_field(H5File &h5file,
		    int*** u,
		    const string field,
		    const int nx,
		    const int ny,
		    const int nz);
void load_h5(const string h5filename,
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
