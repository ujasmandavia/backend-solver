#pragma once
#include "../base_edge.hh"
#include "../pose3d_vertex.hh"
#include "../point3d_vertex.hh"

namespace SLAMSolver{

//Represent the relative pose information between the index -1 and 0
//pose0_SE3_pose1
class RelativePose3DEdge : public BaseEdge<6,Pose3DVertex,Pose3DVertex>{
public:
  RelativePose3DEdge(const Sophus::SE3& pose0_SE3_pose1);
  Eigen::Matrix<double,6,6> compute_left_jacobian_inv_of_SE3(const Eigen::Matrix<double,6,1>& error_se3);
  void compute_errors() override;
  void compute_jacobias() override;

private:
  Sophus::SE3 pose0_SE3_pose1_;
}; //class RelativePose3DEdge

} //namespace SLAMSolver
