#pragma once

#include "eigen_types.h"

namespace SLAMSolver{
/*
 *  Base Interface class for all vertex
 */
 class Vertex{
 public:
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
   /*
    * @brief Constructor of vertex
    * @param minimal_dimension: Dimension of minimal parametrization
    */
  Vertex(const int minimal_dimension);

  // @brief virtual destructor
  virtual ~Vertex();

  // @return Dimension of minimal parametrization
  int minimal_dimension() const { return minimal_dimension_; }

  // @return ID of the vertex
  unsigned long id() const { return id_; }

  /*
  * @brief pure virtual function of updating the parameter( over parametrization)
  * @param delta update value needed to add to parameter
  */
  virtual void add(const Eigen::VectorXd &delta) = 0;

  /*
  * @brief Set the vertex to be fixed or changable
  * @param fixed true if we want to set this vertex fixed
  */

  void set_fixed(bool fixed = true){
    fixed_ = fixed;
  }

  // @return True if the Vertex is fixed
  bool is_fixed() const { return fixed_; }

protected:
  const int minimal_dimension_;  //Dimension of minimal parametrization
  unsigned long id_;  //vertex id

  bool fixed_ = false;  //true if the vertex is fixed 
 };  //class Vertex
} //namespace SLAMSolver
