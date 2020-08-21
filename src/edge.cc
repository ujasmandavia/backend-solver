#include "../vertex.hh"
#include "../edge.hh"
#include <iostream>

namespace SLAMSolver{
Edge::Edge(const int error_dimension, const int num_vertices) : error_dimension_(error_dimension), num_vertices_(num_vertices) {
  static unsigned long uuid_edge = 0;
  id_ = uuid_edge++;

  //Resize the vertex pointer vector to number of vertices
  vertex_interface_ptr_.resize(num_vertices_);

  //Resize the error Vector
  error_.resize(error_dimension_,1)

  //Resize the Jacobian Matrix
  jacobian_.resize(num_vertices_);

  //Set the information Matrix
  Eigen::MatrixXd information(error_dimension_,error_dimension_);
  information.setIdentity();
  information_ = information;
}

std::shared_ptr<Vertex> Edge::get_vertex_interface(const int i){
  //Check the vertex index
  if (i >= num_vertices_){
    std::cout << "[ERROR] Access Vertex Number out of Bound" << std::endl;
    exit(0);
  }
  return vertices_interface_ptr_[i];
}

double Edge::chi2(){
  return error_.transpose() * information_ * error_;
}
} //namespace SLAMSolver
