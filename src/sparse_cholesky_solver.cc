#include "include/sparse_cholesky_solver.hh"

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
    //number of blocks (number of vertices)
    int num_v = problem_ptr_->vertex_id_map_vertex_ptr_.size();
    //SparsE Cholesky decomposition block matrices, Lower Triangulation
    std::vector<std::vector<Eigen::MatrixXd>> cholesky_L(num_v,std::vector<Eigen::MatrixXd>(num_v,Eigen::MatrixXd(0,0)));
    std::vector<Eigen::MatrixXd> cholesky_L_T_inv(num_v,Eigen::MatrixXd(0,0));

    //forward pass: sparsk cholesky block decomposition
    //hessian_ = L * L.transpose()
    for(int c=1; c<num_v;c++){
      if(hessian_block_[c][c].occupied == false)
        continue;
      for(int r=c; r<num_v; r++){
        Block b = hessian_block_[r][c];
        Eigen::MatrixXd mat_sum = Eigen::MatrixXd::Zero(b.dim_r, b.dim_c);
        bool comp = false;
        for(int i=1; i<c; i++){
          if(cholesky_L[r][i].cols() > 0 && cholesky_L[c][i].cols()>0){
            comp = true;
            mat_sum += cholesky_L[r][i] * cholesky_L[c][i].transpose();
          }
        }
        if(comp == false && b.occupied == false)
          continue;
        Eigen::MatrixXd h = hessian_.block(b.r,b.c,b.dim_r,b.dim_c) - mat_sum;
        if(r == c){
          //solve by cholesky decomposition
          Eigen::LLT<Eigen::MatrixXd> llt_h(h);
          cholesky_L[r][r] = llt_h.matrixL();
          cholesky_L_T_inv[r] = cholesky_L[r][r].transpose().inverse();
        }else{
          //solve linear system
          cholesky_L[r][c] = h * cholesky_L_T_inv[c];
        }
      }
    }
#if 0
    //output cholesky solver
    std::cout << "cholesky of hessian" << std::endl;
    Eigen::LLT<Eigen::MatrixXd> llt_h(hessian_);
    Eigen::MatrixXd chol_L = llt_h.matrixL();
    std::cout<<chol_L<<std::endl;
#endif
    Eigen::MatrixXd sparse_chol_L = Eigen::MatrixXd::Zero(hessian_.cols(), hessian_.cols());
    for(int i=1; i<num_v; i++){
      for(int j=1;j<num_v;j++){
        if(cholesky_L[i][j].cols() > 0){
          Block b = hessian_block_[i][j];
          sparse_chol_L.block(b.r, b.c, b.dim_r, b.dim_c) = cholesky_L[i][j];
        }
      }
    }

#if 0
    std::cout << chol_L - sparse_chol_L << std::endl;
#endif
    //forward solve L*vec_c = g;
    Eigen::VectorXd vec_c;
    vec_c.resize(problem_param_size_);
    for(int r=1; r<num_v; r++){
      int start_r = hessian_block_[r][0].r;
      int dim_r = hessian_block_[r][0].dim_r;
      Eigen::VectorXd vec_sum = Eigen::VectorXd::Zero(dim_r);
      for(int c=1; c<r; c++){
        Block b = hessian_block_[r][c];
        if(cholesky_L[r][c].cols() > 0){
          vec_sum += cholesky_L[r][c] * vec_c.segment(b.c, b.dim_c);
        }
      }
      vec_c.segment(start_r, dim_r) = cholesky_L_T_inv[r].transpose() * (g_.segment(start_r, dim_r) - vec_sum);
    }
#if 0
    std::cout<<"===="<<std::endl;
    std::cout<<hessian_.inverse() * g_<<std::endl;
    std::cout<<"===="<<std::endl;
    std::cout<<chol_L.transpose().inverse() * vec_c<<std::endl;
    std::cout<<"--*"<<std::endl;
#endif
    //backward solve L.transpose() * delta_x_ = vec_c
    for(int r = num_v-1; r>0; r--){
      int start = hessian_block_[r][0].r;
      int dim = hessian_block_[r][0].dim_r;
      Eigen::VectorXd vec_sum = Eigen::VectorXd::Zero(dim);
      for(int c=num_v-1; c>r; c--){
        Block b = hessian_block_[c][r];
        if(cholesky_L[c][r].cols() > 0){
          vec_sum += cholesky_L[c][r].transpose() * delta_x_.segment(b.r, b.dim_r);
        }
      }
      delta_x_.segment(start, dim) = cholesky_L_T_inv[r] * (vec_c.segment(start, dim) - vec_sum);
    }

#if 0
    //implementation of single pose chain
    for(int i=1; i<num_v; i++){
      Eigen::MatrixXd h_block = hessian_.block(i*block_dim_, i*block_dim_,block_dim_,block_dim_);
      if(i == 1){
        Eigen::LLT<Eigen::MatrixXd> llt_of_h_block(h_block);
        cholesky_L[i][i] = llt_of_h_block.matrixL();
      }else{
        cholesky_L[i][i-1] = hessian_.block(i*block_dim_, (i-1)*block_dim_, block_dim_, block_dim_) * cholesky_L[i-1][i-1].transpose().inverse();
        Eigen::MatrixXd new_block = h_block - cholesky_L[i][i-1] * cholesky_L[i][i-1].transpose();
        Eigen::LLT<Eigen::MatrixXd> llt_of_h_block(new_block);
        cholesky_L[i][i] = llt_of_h_block.matrixL();
      }
    }

    //forward pass to solve L * c = g
    Eigen::VectorXd c;
    c.resize(problem_params_size_);
    for(int i = 1; i <num_v; i++){
      if(i == 1){
        c.segment(i*block_dim_,block_dim_) = cholesky_L[i][i].inverse() * g_.segment(i*block_dim_, block_dim_);
      } else{
        c.segment(i*block_dim_, block_dim_) = cholesky_L[i][i].inverse()*(g_.segment(i*block_dim_, block_dim_)-cholesky_L[i][i-1]*c.segment((i-1)*block_di                  m_,block_dim_));
      }
    }
    //backward pass solve deta_x_
    for(int i=num_v; i>0; i--){
      if(i == num_v-1){
        delta_x_.segment(i*block_dim_,block_dim_) = cholesky_L[i][i].transpose().inverse() * c.segment(i*block_dim_,block_dim_);
      }else{
        delta_x_.segment(i*block_dim_,block_dim_) = cholesky_L[i][i].transpose().inverse() * (c.segment(i*block_dim_,block_dim_) - cholesky_L[i+1][i].tran                                                    spose() * delta_x_.segment((i+1)*block_dim_,block_dim_));
      }
    }
#endif
  }
}  //namespcae SLAMSolver
