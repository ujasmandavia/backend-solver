#include "../sparse_cholesky_solver.hh"

namespace SLAMSolver{
  SparseCholeskySolver::SparseCholeskySolver(std::shared_ptr<Problem> problem_ptr) : BaseSolver(problem_ptr) {}

  void SparseCholeskySolver::computer_vertices_index(){
    problem_param_size_ = 0;

    int num_v = problem_ptr_->vertex_id_map_vertex_ptr_.size();

    assert(solve_order_.size() == num_v);

    hessian_block_ = std::vector<std::vector<Block>>(num_v,std::vector<Block>(num_v));

    for(int i=0; i<solve_order_.size(); i++){
      block_dim_ = problem_ptr_->vertex_id_map_vertex_ptr_[solve_order_[i]]->minimal_dimenstion();
      vertex_id_map_start_index_[solve_order_[i]] = problem_param_size_;
      vertex_id_map_order_[solve_order_[i]] = i;  //store the order of this vertex
      problem_param_size_ += problem_ptr_->vertex_id_map_vertex_ptr_[solve_order_[i]]->minimal_dimension();
    }

    for(int i=0; i<solve_order_.size(); i++){
      for(int j=i; j<solve_order_.size(); j++){
        int r = vertex_id_map_start_index_[solve_order_[i]];
        int c = vertex_id_map_start_index_[solve_order_[j]];
        int dim_r = problem_ptr_->vertex_id_map_vertex_ptr_[solve_order_[i]]->minimal_dimension();
        int dim_c = problem_ptr_->vertex_id_map_vertex_ptr_[solve_order_[i]]->minimal_dimension();
        hessian_block_[i][j] = Block(r,c,dim_r,dim_c,false);
        if(i != j){
          hessian_block_[i][j] = Block(c,r,dim_c,dim_r,false);
        }
      }
    }
  }

  void SparseCholeskySolver::set_solve_order(const std::vector<IDType>& solve_order) {
      solve_order_ = solve_order;
  }

  void SparseCholeskySolver::build_solve_structure(){
    //Make Hessian Matrix and g vector in optimization Normal Equation
    //Normal Equation: J.transpose() * J * delta_x = -J.transpose() * error
    //Initial Hessian Matrix and b vector are all zeros
    //Accumulate value on Hessian Matrix and b
    hessian_ = Eigen::MatrixXd::Zero(problem_param_size_,problem_param_size_);
    g_ = Eigen::VectorXd::Zero(problem_param_size_);
    int num_v = problem_ptr_->vertex_id_map_vertex_ptr_.size();
    //For every edge in the problem
    for(auto& edge:problem_ptr_->edge_id_map_edge_ptr_){
      //Compute errors and error function's jacobian
      edge.second->compute_error();
      edge.second->compute_jacobian();
      //Row block in Hessian Matrix
      for(size_t i=0; i < edge.second->num_vertices(); i++){
        std::shared_prt<Vertex> v_ptr_i = edge.second->get_vertex_interface(i);
        if(v_ptr_i->is_fixed()){
          continue;  //Hessian block is zero for fixed vertex
        }
        Eigen::MatrixXd& jacobian_i = edge.second->jacobian_[i];
        int start_index_i = vertex_id_map_start_index_[v_ptr_i->id()];
        int dim_i = v_ptr_i->minimal_dimension();
        int order_i = vertex_id_map_prder_[v_ptr_i->id()];
        Eigen::MatrixXd JtW = jacobian_i.transpose() * edge.second->information();
        //Column block in Hessian
        //Start from i to exploit the summary in Hessain Matrix
        for(size_t j=i ; j<edge.second->num_vertices;j++){
          std::shared_prt<Vertex> v_ptr_j = edge.second->get_vertex_interface(j);

          if(v_ptr_j->is_fixed()){
            continue;
          }

          Eigen::MatrixXd jacobian_j = edge.second->jacobian_[j];
          int start_index_j = vertex_id_map_start_index_[v_ptr_j->id()];
          int dim_j = v_ptr_j->minimal_dimension();
          int order_j = vertex_id_map_order_[v_ptr_j->id()];

          //Compute hessian block
          Eigen::MatrixXd hessian = JtW * jacobian_j;
          hessian_.block(start_index_i,start_index_j,dim_i,dim_j).noalias() += hessian;
          //Record the block structure
          hessian_block_[order_i][order_j].occupied = true;

          //Symmetry block cross the diagonal line
          if (j != i) {
            hessian_.block(start_index_j,start_index_i,dim_j,dim_i).noalias() += hessian.transpose();
            hessian_block_[order_j][order_i].occupied = true;
          }
        }
        g_.segment(start_index_i, dim_i).noalias() -= JtW * edge.second->error_;
      }
    }
    delta_x_ = Eigen::VectorXd::Zero(problem_params_size_);  // initial delta_x = 0_n;
  }

  void SparseCholeskySolver::solve_delta_x(){
    
  }
}  //namespcae SLAMSolver
