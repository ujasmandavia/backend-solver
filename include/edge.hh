#pragma once

#include <iostream>
#include <memory>
#include "eigen_types.h"

class Vertex;

/*
* @brief Base class for all edges
*/
class Edge{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  /*
  * @brief Constructor of Edge
  * @param error_dimension: Dimension of error vector
  * @param num_vertices: Number of vertices connected by this edge
  */
  Edge(const int error_dimension, const int num_vertices);

  //Virtual destructor of edge class
  virtual ~Edge();

  // @return id of the edge
  unsigned long id()const { return id_; }

  // @return Dimensions of error vector
  int error_dimension() const { return error_dimension_; }

  //@return Number of vertices connected to this edge
  int num_vertices() const { return num_vertices_; }

  /*
  * @brief get Vertex class pointer (not BaseVertex template)\
  * @param i the index of vertex
  * @return shared pointer to vertex
  */
  std::shared_ptr<Vertex> get_vertex_interface(const int i);

  /*
  * @brief pure virtual function of computing error, every Edge class needs to implement this
  */
  virtual void compute_error() = 0;

  /*
  * @brief pure virtual function of computing jacobian, every Edge class needs to implement this
  */
  virtual void compute_jacobian() = 0;

  //return cost
  // error_.transpose() * information_ * error_;
  double chi2();

  // @brief set information matrix in the edge
  void set_information(const Eigen::MatrixXd &information){ information_ = information; }

  // @return information matric
  Eigen::MatrixXd information() const { return information_; }

  //Error vector and Jacobian Matrix
  Eigen::VectorXd error_;                      //error
  std::vector<Eigen::MatrixXd> jacobian_;      // jacobian

protected:

  const int error_dimension_;   //Dimenion of error
  const int num_vertices_;  //number of vertices connected to this edge
  unsigned long id_;    //id of each edge
  std::vector<std::shared_ptr<Vertex>> vertices_interface_ptr_; //interface pointer of type vertex
  Eigen::MatrixXd information_;  //information matrix
}; //class Edge
