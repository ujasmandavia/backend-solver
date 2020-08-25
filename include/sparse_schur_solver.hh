#pragma once

#include "../base_colver.hh"
#include <set>

namespace SLAMSolver{
///Parameter block in Hessian Matrix
struct Parameters{
  Parameters(int start_index, int minimal_dim) : start_index_(start_index), minimal_dim_(minimal_dim) {}
  int start_index_;
  int minimal_dim_;
};

/*
 * @brief Solver using sparse schur decomposition to solve normal equation
 * Assume there is no edge between any pair of vertices to be marginalized
 */
class SparseSchurSolver : public BaseSolver{
  public:
    ///@ brief Constructor passes problem_ptr to base solver
    SparseSchurSolver(std::shared_ptr<Problem> problem_ptr);

    ///@brief Reorder the vertices, make all parameters of marginalized vertices to bottom-right
    void compute_vertices_index() override;

    ///@brief Solve normal equation using sparse schur solver
    void solve_delta_x() override;

    /*
     * @brief Set vertices that need to be marginalized
     *  @param : marginalized_vertices A set containing id for all vertices that want to marginalize
     */
    void set_marginalized_vertices(std::set<IDType>& marginalized_vertices);

  public:
    ///@brief Compute sparse schur complements that need to be performed on hessian matrix and g vector
    void compute_schur_complements();
    int marginalized_param_size_; // size of parameters that want to marginalize
    std::set<IDType> marginalized_vertices_;
    /// represent the schur complement matrix block operation
    /// First in pair represent the parameters that need to be marginalized out
    //  Second in pair represent a set of parameters that need to complement when marginalizing out the pair.first
    std::vector<std::pair<Parameters, std::vector<Parameters>>> schur_complements_;
};
}  //namespace SLAMSolver
