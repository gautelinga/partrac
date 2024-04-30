#ifndef __XDMF_HELPERS_HPP
#define __XDMF_HELPERS_HPP

#include "typedefs.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "dolfin/io/HDF5Interface.h"

void read_dataset_scalar(std::string& h5filename, std::string& field, std::vector<double>& data);

void reorder_indices(std::vector<double>& data_, const std::vector<double>& xdata, const std::vector<Uint>& j2i, const int dim);

void read_dataset_vector(std::string& h5filename_u, std::string& field, std::vector<double>& data, const int dim);

std::vector<std::pair<double, std::vector<std::string>>> parse_xdmf(const std::string& xdmffilename, 
                                                                    std::string& h5filename,
                                                                    std::string& topology_path, 
                                                                    std::string& geometry_path);

std::vector<std::pair<double, std::vector<std::string>>> parse_xdmf(const std::string& xdmffilename);

#endif