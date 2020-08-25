#include "../relative_pose3d_edge.hh"

namespace SLAMSolver{

RelativePose3DEdge::RelativePose3DEdge(const Sophus::SE3 pose0_SE3_pose1) : pose0_SE3_pose1_(pose0_SE3_pose1) {}

void RelativePose3DEdge::compute_errors(){
  Sophus::SE3 frame_SE3_pose0 = get_vertex<0>->parameters();
  Sophus::SE3 frame_SE3_pose1 = get_vertex<1>->parameters();
  error_ = (pose0_SE3_pose1_.inverse() * frame_SE3_pose0.inverse() * frame_SE3_pose1).log();
}

Eigen:Matrix<double,6,6> RelativePose3DEdge::compute_left_jacobian_inv_of_SE3(const Eigen::Matrix<double,6,1>& error_se3){
  Eigen::Matrix<double,6,6> J_inv_left;
  J_inv_left.block(0,0,3,3) = Sophus::SO3::hat(error_se3.template tail<3>());
  J_inv_left.block(0,3,3,3) = Sophus::SO3::hat(error_se3.template head<3>());
  J_inv_left.block(3,0,3,3) = Eigen::Matrix3d::Zero();
  J_inv_left.block(3,3,3,3) = Sophus::SO3::hat(error_se3.template tail<3>());
  J_inv_left = 0.5 * J_inv_left + Eigen::Matrix<double,6,6>::Identity();
  return J_inv_left;
}

void RelativePose3DEdge::compute_jacobians(){
  Sophus::SE3 frame_SE3_pose0 = get_vertex<0>()->parameters();
  Sophus::SE3 frame_SE3_pose1 = get_vertex<1>()->parameters();
  Eigen::Matrix<double,6,1> errors_se3 = error_;
  Eigen::Matrix<double,6,6> J_inv_left = compute_left_jacobian_inv_of_SE3(errors_se3);
  jacobians_[0] = - J_inv_left * frame_SE3_pose1.inverse().Adj();
  jacobians_[1] = J_inv_left * frame_SE3_pose1.inverse().Adj();
}

} //namespace SLAMSolver
