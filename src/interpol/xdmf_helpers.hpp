#ifndef __XDMF_HELPERS_HPP
#define __XDMF_HELPERS_HPP

// The XDMF side of the XDMF loaders: the xml a field's .xdmf holds (the mesh's
// datasets and one grid per stamp) and the values of one stamp. Nothing here
// needs dolfin.

#include <string>
#include <utility>
#include <vector>

// The first ncols of each row of a rank-2 dataset: dolfin writes a field as
// rows of columns and pads a 2D vector with a third column no one reads
void read_dataset_columns(const std::string& h5filename, const std::string& field,
                          std::vector<double>& data, const int ncols);

// Each grid's time and (file, dataset), with the mesh's datasets of the first
std::vector<std::pair<double, std::vector<std::string>>> parse_xdmf(const std::string& xdmffilename,
                                                                    std::string& h5filename,
                                                                    std::string& topology_path,
                                                                    std::string& geometry_path);

std::vector<std::pair<double, std::vector<std::string>>> parse_xdmf(const std::string& xdmffilename);

#endif
