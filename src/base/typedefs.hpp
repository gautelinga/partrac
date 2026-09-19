#ifndef __TYPEDEFS_HPP
#define __TYPEDEFS_HPP

#include <vector>
#include <list>
#include <array>
#include <map>
#include <eigen3/Eigen/Dense>
#include <memory>

typedef std::size_t Uint;

// Inline everything in per-particle loops
#define PARTRAC_HOT_LOOP __attribute__((flatten))

class EdgeType {
public:
    EdgeType (const std::array<Uint, 2>& a, const double b) : first(a), second(b) {};
    EdgeType (const std::array<Uint, 2>& a, const double b, const double tau, const double rho_prev) : first(a), second(b), tau(tau), rho_prev(rho_prev) {};
    std::array<Uint, 2> first;
    double second;
    double tau = 0.0;
    double rho_prev = 1.0;
};

class FaceType {
public:
    FaceType (const std::array<Uint, 3>& a, const double b) : first(a), second(b) {};
    FaceType (const std::array<Uint, 3>& a, const double b, const double tau, const double rho_prev) : first(a), second(b), tau(tau), rho_prev(rho_prev) {};
    std::array<Uint, 3> first;
    double second;
    double tau = 0.0;
    double rho_prev = 1.0;
};

//typedef std::vector<std::pair<std::array<Uint, 2>, double>> EdgesType;
typedef std::vector<EdgeType> EdgesType; 
//typedef std::vector<std::pair<std::array<Uint, 3>, double>> FacesType;
typedef std::vector<FaceType> FacesType;
typedef std::list<Uint> FacesListType;
typedef std::vector<Uint> EdgesListType;   // template edge -> live edge
typedef std::vector<Uint> NodesListType;   // template node -> live node
// Adjacency tables: the faces of an edge, the edges of a node
typedef std::vector<Uint> AdjRowType;
typedef std::vector<AdjRowType> Edge2FacesType;
typedef std::vector<AdjRowType> Node2EdgesType;
typedef std::vector<std::map<Uint, double>> InteriorAnglesType;
typedef Eigen::Vector3d Vector3d;
typedef Eigen::Matrix3d Matrix3d;

#endif
