#include "minimizer_levenberg_marquardt.hh"
#include <iostream>
#include <fstream>

namespace SLAMSolver{

MinimizerLevenbergMarquardt::MinimizerLevenbergMarquardt(std::shared_ptr<BaseSolver> solver_ptr, const LevenbergMarquardtConfig& config) :
                                                         solver_ptr_(solver_ptr), config_(config) {}

double MinimizerLevenbergMarquardt::compute_initial_lambda(){
  //Check Hessian matrix
  if(solver_ptr_->hessian_.cols() != solver_ptr_->hessian_.rows()){
    std::cerr << "[ERROR] Hessian Matrix is not a square!" << std::endl;
    exit(0);
  }

  double max_hessian_diag = 0;
  for(int i = 0; i < solver_ptr_->hessian_.cols(); i++){
    max_hessian_diag = std::max(max_hessian_diag, solver_ptr_->hessian_(i,i));
  }
  return config_.tau * max_hessian_diag;
}

void MinimizerLevenbergMarquardt::add_lambda_to_hessian(const double lambda){
  //Check hessian matrix
  if(solver_ptr_->hessian_.cols() != solver_ptr_->hessian_.rows()){
    std::cerr << "[ERROR] Hessian Matrix is not a square!" << std::endl;
    exit(0);
  }

  for(int i = 0; i < solver_ptr_->hessian_.cols(); i++) {
    solver_ptr_->hessian_(i,i) += lambda;
  }
}

void MinimizerLevenbergMarquardt::remove_lambda_from_hessian(const double lambda) {
  //Check Hessian Matrix
  if(solver_ptr_->hessian_.cols() != solver_ptr_->hessian_.rows()){
    std::cout<<"[ERROR] Hessian Matrix is not Square"<<std::endl;
    exit(0);
  }
  for(int i = 0; i < solver_ptr_->hessian_.cols(); i++) {
    solver_ptr_->hessian_(i,i) -= lambda;
  }
}

void MinimizerLevenbergMarquardt::compute_gain_factor(const double cost_before_update, const double lambda){
  //compute the cost after update
  double cost_after_update = solver_ptr_->compute_total_cost();
  //Compute gauss newton approximation cost
  double gn_cost = solver_ptr_->delta_x_.transpose() * (lambda * solver_ptr_->delta_x_ + solver_ptr_->g_) + 0.001;
  return (cost_before_update - cost_after_update) / gn_cost;
}

bool MinimizerLevenbergMarquardt::minimize(const int max_iterations){
  int iter = 0;   //optimization iterations
  double v = 2;

  //Gauss-Newton, build Hessian Matrix and b vector in normal equation
  solver_ptr_->computer_vertices_index();
  solver_ptr_->build_solve_structure();

  //Compute initial value of lambda in LM
  double lambda = compute_initial_lambda();
  std::cout << "Initial Lambda: " << lambda << std::endl;
  //Threshold for stopping criteria
  bool found = solver_ptr_->g_.lpNorm<Eigen::Infinity>() < config_.min_cost_function_gradient;    //TODO
  while(!found && iter < max_iterations){
    iter++;
    bool update_successful = false;
    int fail_count = 0;
    while(!update_successful){
      //1. Add lambda to the hessian and construct the normal equation
      std::cout << "The current lambda is: " << lambda << std::endl;
      add_lambda_to_hessian(lambda);

      //2. Solver linear system to obtain delta_x
      //(Hessian + lambda*I) * delta_x  = g
      solver_ptr_->solve_delta_x();
      std::cout << "Delta x update norm" << std::endl;
      std::cout << solver_ptr_->delta_x_.norm() << std::endl;
      if(solver_ptr_->delta_x_.norm() <= config_.min_norm_delta_x || fail_count > 10){
        //Stop criteria: delta_x is too small
        found = true;
        return false;
      }else{
        // 1. Cache total cost before update
        double cost_before_update = solver_ptr_->compute_total_cost();
        std::cout<<"Cost: "<<cost_before_update<<std::endl;
        //2. Update all the vertices using delta_x;
        solver_ptr_->update_vertices_with_delta_x();
        //3. compute gain factor
        double gain = compute_gain_factor(cost_before_update,lambda);
        if(gain > 0){
          std::cout << "A good step" << std::endl;
          //A good update, cost goes down
          update_successful = true;
          fail_count = 0;
          // 1. Check if satisfied the first stopping criteria
          if (solver_ptr_->g_.lpNorm<Eigen::Infinity>() < config_.min_cost_function_gradient) {
            found = true;
            continue;
          }
          // 2. Compute new heassian and g with Guass Newton method
          solver_ptr_->build_solve_structure();
          //3. update lambda
          double factor = 1 - std::pow(2 * gain - 1, 3);
          std::cout<<gain<<std::endl;
          std::cout<<factor<<std::endl;
          lambda = lambda * std::max(0.333, factor);
          v = 2;
        }else{
          std::cout << "Not a good step, roll back" << std::endl;
          //Not a good step
          update_successful = false;
          fail_count++;
          // 1.roll back to original parameter
          solver_ptr_->rollback_verteices();
          // 2.Remove previous lambda from Hessian Matrix
          remove_lambda_from_hessian(lambda);
          // 3. Update lambda and v;
          lambda *= v;
          v *= 2;
        }
      }
    }
  }
  return true;
}

} //namespace SLAMSolver
