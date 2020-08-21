#pragma once
#include "../base_vertex.hh"
#include "sophus/se3.hpp"

namespace SLAMSolver{

/*
 * @brief Graph optimization for 3D Pose as a vertex
*/
class Pose3DVertex : public BaseVertex<6,Sophus::SE3>{
public:
  /// @brief Default constructor
  /// @param combined True if we can combine the translation and rotation in optimization
  Pose3DVertex(const bool combined);

  /// @brief default destructor
  ~Pose3DVertex() = default;

  /// @brief Override the update function for 3D pose
  /// @note Use left multiplication of perturbation and separately update rotation and translation
  void add(Eigen::VectorXd& delta) override;

private:
  // True if we want to use full SE3
  // False if we want to optimize seperately on rotation and translation
  bool combined_ = true;
}

} //namespace`SLAMSolver
