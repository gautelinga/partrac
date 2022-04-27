#ifndef __EXP_TYPEDEFS_HPP
#define __EXP_TYPEDEFS_HPP

#include <vector>
#include <list>
#include <array>
#include <map>
#include <eigen3/Eigen/Dense>
#include <memory>

#include "H5Cpp.h"

typedef std::size_t Uint;

typedef double Real;
typedef Eigen::Vector3d Vector;
typedef Eigen::Matrix3d Matrix;
#define MPI_Real MPI_DOUBLE

//typedef H5::H5File H5File;
//typedef H5::DataSpace DataSpace;
//typedef H5::DataSet DataSet;

#endif