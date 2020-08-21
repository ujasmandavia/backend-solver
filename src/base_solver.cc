#include "../base_solver.h"

namespace SLAMSolver{
BaseSolver::BaseSolver(std::shared_ptr<Problem> problem_ptr) : problem_ptr_(problem_ptr) {}

void BaseSolver::update_vertices_with_delta_x(){
  for(auto it = problem_ptr_->vertex_id_map_vertex_ptr_.begin(); it != problem_ptr_->vertex_id_map_vertex_ptr_.end(); it++){
    auto v_ptr = it->second;
    int start_index = vertex_id_map_start_index_[it->first];
    int dim = v_ptr->minimal_dimension();
    //Fetch this vertex's corresponding update vector from delta_x
    Eigen::VectorXd delta = delta_x_.segment(start_index,dim);
    //UPpdate the vertex through plus interface
    v_ptr->add(delta);
  }
}

void BaseSolver::rollback_vertices(){
  for(auto it = problem_ptr_->vertex_id_map_vertex_ptr_.begin(); it != problem_ptr_->vertex_id_map_vertex_ptr_.end(); it++){
    auto v_ptr = it->second;
    int start_index = vertex_id_map_start_index_[it->first];
    int dim = v_ptr->minimal_dimension();
    Eigen::VectorXd delta = delta_x_.segment(start_index,dim);
    v_ptr->add(-delta);
  }
}

double  BaseSolver::compute_total_cost(){
  return problem_ptr_->computer_cost();
}

void BaseSolver::build_solve_structure(){
  //Make Hessian matrix and g vector in Optimization Normal Equnation
  //Normal equation J.transpose() * J * delta_x = -J.transpose() * error
  //Initial Hessian matrix and b vector are all zeros
  //Accumulate value on Hessian and b vector
  hessian_ = Eigen::MatrixXd::Zero(problem_param_size_,problem_param_size_);
  g_ = Eigen::VectorXd::Zero(problem_param_size_);

  //For every edge in the problem
  for(auto& edge:problem_ptr_->edge_id_map_edge_ptr_){
    //compute errors and error function's jacobian
    edge.second->compute_error();
    edge.second->compute_jacobian();
    //Row block in hessian
    for(size_t i=0; i<edge.second->num_vertices(); i++){
      std::shared_ptr<Vertex> v_ptr_i = edge.second->get_vertex_interface(i);
      if(v_ptr_i->is_fixed())
        continue;

      Eigen::MatrixXd& jacobian_i = edge.second->jacobian_[i];
      int start_index_i = vertex_id_map_start_index_[v_ptr_i->id()];
      int dim_j = v_ptr_i->minimal_dimension();
      Eigen::MatrixXd JtW = jacobian_i.transpose() * edge.second->information();

      //Column block in Hessian
      //start from i to exploit the symmetry in the Hessian Matrix
      for(size_t j=i; j<edge.second->num_vertices(); j++){
        std::shared_ptr<Vertec> v_ptr_j = edge.second->get_vertex_interface();
        if(v_ptr_j->is_fixed())
          continue;

        Eigen::MatrixXd& jacobian_j = edge.second->jacobian_[j];
        int start_index_j = vertex_id_map_start_index_[v_ptr_j->id()];
        int dim_j = v_ptr_j->minimal_dimension();
        Eigen::MatrixXd hessian = JtW * jacobian_j;
        hessian_.block(start_index_i,start_index_j,dim_i,dim_j).noalias() += hessian;
        //Symmetry block cross the diagonal line

        if(j != i)
          hessian_.block(start_index_j,start_index_i,dim_j,dim_i).noalias() += hessian.transpose();
      }
      g_.segment(start_index_i, dim_i).noalias() -= JtW * edge.second->error_;
    }
  }
  delta_x_ = Eigen::VectorXd::Zero(problem_param_size_);  // initial delta_x = 0_n;
}
} //namespace SLAMSolver
