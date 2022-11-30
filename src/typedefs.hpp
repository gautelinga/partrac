#ifndef __TYPEDEFS_HPP
#define __TYPEDEFS_HPP

#include <vector>
#include <list>
#include <array>
#include <map>
#include <eigen3/Eigen/Dense>
#include <memory>
#include "H5Cpp.h"

typedef std::size_t Uint;

class EdgeType {
public:
    EdgeType (const std::array<Uint, 2>& a, const double b) : first(a), second(b) {};
    EdgeType (const std::array<Uint, 2>& a, const double b, const double tau, const double rho_prev) : first(a), second(b), tau(tau), rho_prev(rho_prev) {};
    std::array<Uint, 2> first;
    double second;
    double tau = 0.0;
    double rho_prev = 1.0;
};

//typedef std::vector<std::pair<std::array<Uint, 2>, double>> EdgesType;
typedef std::vector<EdgeType> EdgesType; 
typedef std::vector<std::pair<std::array<Uint, 3>, double>> FacesType;
typedef std::list<Uint> FacesListType;
typedef std::list<Uint> EdgesListType;
typedef std::list<Uint> NodesListType;
typedef std::vector<FacesListType> Edge2FacesType;
typedef std::vector<EdgesListType> Node2EdgesType;
typedef std::vector<std::map<Uint, double>> InteriorAnglesType;
typedef Eigen::Vector3d Vector3d;
typedef Eigen::Matrix3d Matrix3d;

typedef std::shared_ptr<H5::H5File> H5FilePtr;

#endif
