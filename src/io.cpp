#include <iostream>
#include <vector>
#include <map>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <filesystem>
#include "H5Cpp.h"
#include "io.hpp"

using namespace std;

map<string, string> load_settings(string casefile){
  ifstream input(casefile);
  string line;
  string el;
  map<string, string> settings;
  
  while (getline(input, line)){
    size_t pos = line.find("=");

    // cout << line.substr(0, pos) << ": " << line.substr(pos+1) << endl;
    settings[line.substr(0, pos)] = line.substr(pos+1);
  }
  return settings;
}

vector<string> get_files(const string& s)
{
    vector<string> r;
    for(auto& p : filesystem::directory_iterator(s))
        if (!p.is_directory())
            r.push_back(p.path().string());
    return r;
}

vector<vector<string>> load_grid(string infile){
  ifstream input(infile);

  vector<char> x;
  string line;

  vector<vector<string>> sites;
  int nx = 0;
  while (getline(input, line)){
    // cout << line << endl;
    vector<string> sites_loc;
    boost::split(sites_loc, line, boost::is_any_of(" "));
    int sites_loc_size = sites_loc.size();
    nx = max(nx, sites_loc_size);
    sites.push_back(sites_loc);
  }
  int ny = sites.size();

  cout << "Size of grid: " << nx << " x " << ny << endl;
  for (int iy=0; iy < ny; ++iy){
    int nx_loc = sites[iy].size();
    for (int j=0; j < nx-nx_loc; ++j){
      sites[iy].push_back("s");
    }
  }
  return sites;
}

vector<vector<string>> load_fields(string infile){
  ifstream input(infile);

  vector<char> x;
  string line;

  vector<vector<string>> data;
  int nx = 0;
  while (getline(input, line)){
    // cout << line << endl;
    vector<string> data_loc;
    boost::split(data_loc, line, boost::is_any_of(" "));
    int data_loc_size = data_loc.size();
    nx = max(nx, data_loc_size);
    data.push_back(data_loc);
  }
  return data;
}

void print_grid(bool** grid, const int nx, const int ny){
  for (int iy=0; iy < ny; ++iy){
    for (int ix=0; ix < nx; ++ix){
      cout << grid[iy][ix];
    }
    cout << endl;
  }
}

void copy_arr(double*** f,
	      double*** f_prev,
	      bool** grid,
	      const int nx,
	      const int ny,
	      const int nc){
  for (int iy=0; iy < ny; ++iy){
    for (int ix=0; ix < nx; ++ix){
      if (grid[iy][ix]){
	for (int ic=0; ic < nc; ++ic){
	  f_prev[iy][ix][ic] = f[iy][ix][ic];
	}
      }
    }
  }
}

void copy_arr(double** a_from,
	      double** a_to,
	      bool** grid,
	      const int nx,
	      const int ny){
  for (int iy=0; iy < ny; ++iy){
    for (int ix=0; ix < nx; ++ix){
      if (grid[iy][ix]){
	a_to[iy][ix] = a_from[iy][ix];
      }
    }
  }
}

void dump2file(string filename,
	       double** rho,
	       double** m_x,
	       double** m_y,
	       bool** grid,
	       const int nx,
	       const int ny){
  ofstream outfile(filename);
  outfile.precision(17);
  for (int iy=0; iy < ny; ++iy){
    for (int ix=0; ix < nx; ++ix){
      if (grid[iy][ix]){
	outfile << ix << " " << iy << " " << rho[iy][ix] << " " << m_x[iy][ix] << " " << m_y[iy][ix] << endl;
      }
    }
  }
  outfile.close();
}

string create_folder(const string folder){
  if (!filesystem::is_directory(folder)){
    filesystem::create_directory(folder);
  }
  return folder;
}

void verify_file_exists(const string infilename){
  if (!filesystem::exists(infilename)){
    std::cout << "No such file: " << infilename << std::endl;
    exit(0);
  }
}

void print_param(const string key, const double val){
  std::cout << key << " = " << val << std::endl;
}

void print_param(const string key, const int val){
  std::cout << key << " = " << val << std::endl;
}

void print_param(const string key, const string val){
  std::cout << key << " = " << val << std::endl;
}

void write_param(ofstream &ofile, string key, double val){
  ofile << key << "=" << val << std::endl;
}

void write_param(ofstream &ofile, string key, int val){
  ofile << key << "=" << val << std::endl;
}

void write_param(ofstream &ofile, string key, string val){
  ofile << key << "=" << val << std::endl;
}

void load_field(H5File &h5file,
		double*** u,
		const string field,
		const int nx, const int ny, const int nz){
  DataSet dset = h5file.openDataSet(field);
  DataSpace dspace = dset.getSpace();
  vector<double> Uv(nx*ny*nz);
  dset.read(Uv.data(), PredType::NATIVE_DOUBLE,
	    dspace, dspace);
  for (int ix=0; ix<nx; ++ix){
    for (int iy=0; iy<ny; ++iy){
      for (int iz=0; iz<nz; ++iz){
  	u[ix][iy][iz] = Uv[nx*ny*iz+nx*iy+ix];
      }
    }
  }
}

void load_field(H5File &h5file,
		int*** u,
		const string field,
		const int nx, const int ny, const int nz){
  DataSet dset = h5file.openDataSet(field);
  DataSpace dspace = dset.getSpace();
  vector<int> Uv(nx*ny*nz);
  dset.read(Uv.data(), PredType::NATIVE_INT,
	    dspace, dspace);
  for (int ix=0; ix<nx; ++ix){
    for (int iy=0; iy<ny; ++iy){
      for (int iz=0; iz<nz; ++iz){
  	u[ix][iy][iz] = Uv[nx*ny*iz+nx*iy+ix];
      }
    }
  }
}

void load_h5(const string h5filename,
	     double*** ux,
	     double*** uy,
	     double*** uz,
	     double*** rho,
	     double*** p,
	     const int nx,
	     const int ny,
	     const int nz,
	     const bool verbose){
  // Assert that h5 file exists
  verify_file_exists(h5filename);
  if (verbose)
    cout << "Opening " << h5filename << endl;
  H5File h5file(h5filename, H5F_ACC_RDONLY);
  load_field(h5file, ux, "u_x", nx, ny, nz);
  load_field(h5file, uy, "u_y", nx, ny, nz);
  load_field(h5file, uz, "u_z", nx, ny, nz);
  load_field(h5file, rho, "density", nx, ny, nz);
  load_field(h5file, p, "pressure", nx, ny, nz);
  h5file.close();
}
