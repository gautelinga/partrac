#ifndef __MPIWRAP_HPP
#define __MPIWRAP_HPP

#include <mpi.h>
#include "typedefs.hpp"

class MPIwrap {
public:
  MPIwrap(int argc, char* argv[]);
  ~MPIwrap();
  void finalize();
  int rank() const { return m_rank; };
  int size() const { return m_size; };
  //MPI_Comm& comm() { return MPI_COMM_WORLD; };
  //MPI_Info& info() { return MPI_INFO_NULL; };
  std::vector<int> gather(int);
  int scatter(std::vector<int>&);
  //std::vector<T> gatherVec(std::vector<T>& vec, MPI_Datatype type);
  //std::vector<int> gatherVec(std::vector<int>& vec);
  EdgesType wrapEdges(EdgesType& edges, int id_offset);
  FacesType wrapFaces(FacesType& faces, int id_offset);
  void barrier() { MPI_Barrier(MPI_COMM_WORLD); };
  template <typename T> T sum(T &a, MPI_Datatype type){
    T a_out;
    MPI_Reduce(&a, &a_out, 1, type, MPI_SUM, 0, MPI_COMM_WORLD);
    return a_out;
  }
  template <typename T> T allsum(T &a, MPI_Datatype type){
    T a_out;
    MPI_Allreduce(&a, &a_out, 1, type, MPI_SUM, MPI_COMM_WORLD);
    return a_out;
  }
private:
  //MPI_Comm m_comm = MPI_COMM_WORLD;
  //MPI_Info m_info  = MPI_INFO_NULL;
  int m_rank;
  int m_size;
};

MPIwrap::MPIwrap(int argc, char* argv[]){
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &m_size);
}

MPIwrap::~MPIwrap(){
  int flag = 0;
  MPI_Finalized(&flag);
  if ( !flag ){
    finalize();
  }
}

void MPIwrap::finalize(){
  MPI_Finalize();
}

std::vector<int> MPIwrap::gather(int value){
  std::vector<int> values(size());
  MPI_Gather(&value, 1, MPI_INT, values.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  return values;
}

/*
std::vector<int> MPIwrap::gatherVec(std::vector<int>& vec){
  int nelements = (int)vec.size();
  std::vector<int> num_elem_ = gather(nelements);

  // Displacements in the receive buffer for MPI_GATHERV
  int* disps = new int[size()];
  for (int i = 0; i < size(); i++)
    disps[i] = (i > 0) ? (disps[i-1] + num_elem_[i-1]) : 0;
  
  int n_total = disps[size()-1] + num_elem_[size()-1];
  std::vector<int> all_vec;
  if (m_rank == 0)
    all_vec.resize(n_total);
  MPI_Gatherv(vec.data(), nelements, MPI_INT,
              all_vec.data(), num_elem_.data(), disps, MPI_INT, 0, m_comm);
  return all_vec;
}*/

template <typename T>
std::vector<T> gather_vector(MPIwrap& mpi, std::vector<T>& vec, MPI_Datatype type){
  int nelements = (int)vec.size();
  std::vector<int> num_elem_ = mpi.gather(nelements);

  // Displacements in the receive buffer for MPI_GATHERV
  int* disps = new int[mpi.size()];
  for (int i = 0; i < mpi.size(); i++)
    disps[i] = (i > 0) ? (disps[i-1] + num_elem_[i-1]) : 0;
  
  int n_total = disps[mpi.size()-1] + num_elem_[mpi.size()-1];
  std::vector<T> all_vec;
  if (mpi.rank() == 0)
    all_vec.resize(n_total);
  MPI_Gatherv(vec.data(), nelements, type,
              all_vec.data(), num_elem_.data(), disps, type, 0, MPI_COMM_WORLD);
  return all_vec;
}

EdgesType MPIwrap::wrapEdges(EdgesType& edges, int id_offset){
  std::vector<Uint> ids_;
  std::vector<double> vals_;
  for (auto &edge: edges){
    ids_.push_back(edge.first[0]+id_offset);
    ids_.push_back(edge.first[1]+id_offset);
    vals_.push_back(edge.second);
  }
  auto all_ids_ = gather_vector<Uint>(*this, ids_, MPI_UNSIGNED_LONG);
  auto all_vals_ = gather_vector<double>(*this, vals_, MPI_DOUBLE);

  auto length = all_ids_.size()/2;
  assert(length == all_vals_.size());
  //all_val_ = gatherVec(val_);
  EdgesType all_edges;
  for (Uint i=0; i<length; ++i){
    all_edges.push_back({{all_ids_[2*i], all_ids_[2*i+1]}, all_vals_[i]});
  }
  return all_edges;
}

FacesType MPIwrap::wrapFaces(FacesType& faces, int id_offset){
  std::vector<Uint> ids_;
  std::vector<double> vals_;
  for (auto &face: faces){
    ids_.push_back(face.first[0]+id_offset);
    ids_.push_back(face.first[1]+id_offset);
    ids_.push_back(face.first[2]+id_offset);
    vals_.push_back(face.second);
  }
  auto all_ids_ = gather_vector<Uint>(*this, ids_, MPI_UNSIGNED_LONG);
  auto all_vals_ = gather_vector<double>(*this, vals_, MPI_DOUBLE);

  auto length = all_ids_.size()/3;
  assert(length == all_vals_.size());

  FacesType all_faces;
  for (Uint i=0; i<length; ++i){
    all_faces.push_back({{all_ids_[3*i], all_ids_[3*i+1], all_ids_[3*i+2]}, all_vals_[i]});
  }
  return all_faces;
}

int MPIwrap::scatter(std::vector<int>& values){
  assert(values.size() == size());
  int value;
  MPI_Scatter(values.data(), 1, MPI_INT, &value, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return value;
}

#endif