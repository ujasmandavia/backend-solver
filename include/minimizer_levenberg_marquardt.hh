#pragma once
#include "../base_solver.h"

struct LevenbergMarquardtConfig{
  double min_norm_delta_x = 1e-9;
  double tau = 1e-6;
  int max_continuous_fail_update = 10;
  double min_cost_function_gradient = 1e-9;
}; //struct

namespace SLAMSolver{

/*
 * @brief Non-Linear least square minimizer with Levenberg-Marquardt method
*/
class MinimizerLevenbergMarquardt{
public:
  /*
   * @brief Constructor
   * @param solver_ptr A solver pointer containing the pointer
   * @param config Configuration of the minimizer
  */
  MinimizerLevenbergMarquardt(std::shared_ptr<BaseSolver> solver_ptr, LevenbergMarquardtConfig& config);

  /*
   * @brief minimize Minimize the non linear least square cost
   * @param max_iterations maximum number of iterations
   * @return
  */
  bool minimize(const int max_iterations);

private:
  /*
   * Compute initial lambda value used in LM
   * @return lambda value
  */
  double compute_initial_lambda();

  /*
   * @brief Add lambda to the hessian matrix
   * @lambda lambda value hessian += lambda * Identity
  */
  void add_lambda_to_hessian(const double lambda);

  /*
   * @brief remove added lambda from hessian matrix in solver ptr
   * @param lambda Lambda value
  */
  void remove_lambda_from_hessian(const double lambda);

  /*
   * @brief compute gain factor used in LM
   * @param cost_before_update least square cost before updating vertices using delta_x in solver_ptr
   * @param lambda Lambda value
   * @return compute gain factor
  */
  double compute_gain_factor(const double cost_before_update, const double lambda);

private:
  ///base class pointer to solver
  std::shared_ptr<BaseSolver> solver_ptr_;
  ///COnfiguration for LM algorithm
  LevenbergMarquardtConfig config_;
}; //class MinimizerLevenbergMarquardt

} //namespace SLAMSolver
