#pragma once

#include "problem.hh"
#include "eigen_types.h"

namespace SLAMSolver{
/*
 * @brief Base class for all kind of solvers
*/
class BaseSolver{
public:
  using IDType = Problem::IDType;
  /*
   * @brief Constructor for BaseSolver
   * @param problem_ptr The problem need to be solved
  */
  BaseSolver(std::shared_ptr<Problem> problem_ptr);

  virtual ~BaseSolver() = default;

  /*
   * @brief Update all vertices' parameters using the current delta_x_
  */
  void update_vertices_with_delta_x();

  /*
   * @brief Rollback all the vertices' parameters for delta_x_ change
  */
  void rollback_vertices();

  /*
   * @brief Computer current problems total cost
   * @return total_cost
  */
  double compute_total_cost();

  /*
   * @brief pure virtual function for computing each vertex's paramter vector's starting index
   * Different solver need to implement their own starategy for ordering their vertices
  */
  virtual void computer_vertices_index() = 0;

  /*
   * @brief Virtual function for building the normal equation given the vertices' position
   * compute hessian_ and g_
   * Base class implemented a common version of building normal equation
  */
  virtual void build_solve_structure();

  /*
   * pure virtual function for solving the normal equ
   * hessian_ * delta_x_ = g_
  */
  virtual void solve_delta_x() = 0;

public:
  Eigen::MatrixXd hessian_;
  Eigen::VectorXd delta_x_;
  Eigen::VectorXd g_;

protected:
  int problem_param_size_ = 0;
  //Hash map from vertex's ID to vertex's start index in problem's
  std::map<IDType,int> vertex_id_map_start_index_;
  std::shared_ptr<Problem> problem_ptr_;
}; //class BaseSolver
} //namespace SLAMSolver
