#include "../reprojection_edge.hh"

namespace SLAMSolver{
ReprojectionEdge::ReprojectionEdge(const Eigen::Vector2d &uv)
: uv_(uv)
{
}

void ReprojectionEdge::compute_errors() {
  //Get point and camera pose
  Eigen::Vector3d point_world = get_vertex<0>()->parameters();
  Sophus::SE3 world_SE3_camera = get_vertex<1>()->parameters();
  //Transform point from world frame to camera frame
  Eigen::Vector3d point_camera = world_SE3_camera.inverse() * point_world;
  //project
  point_camera /= point_camera(2);
  errors_ = point_camera.segment(0,2) - uv_;
}

void ReprojectionEdge::compute_jacobians() {
  //Get point and camera pose
  Eigen::Vector3d point_world = get_vertex<0>()->parameters();
  Sophus::SE3 world_SE3_camera = get_vertex<1>()->parameters();
  //Transform point from world frame to camera frame
  Eigen::Vector3d point_camera = world_SE3_camera.inverse() * point_world;
  //Reprojection jacobians
  Eigen::Matrix<double,2,3> jacobian_point_world; // jacobian of error w.r.t point in world frame
  Eigen::Matrix<double,2,6> jacobian_pose; // jacobian of error w.r.t pose of camera
  //jacobian matrix of error w.r.t point in camera frame
  Eigen::Matrix<double,2,3> jacobian_error_point_camera;
  jacobian_error_point_camera<<1/point_camera(2), 0, -point_camera(0) / (point_camera(2) * point_camera(2)),
                               0, 1/point_camera(2), -point_camera(1) / (point_camera(2) * point_camera(2));
  jacobian_point_world = jacobian_error_point_camera * world_SE3_camera.so3().matrix().transpose();
  //Jacobian of point in camera frame w.r.t pose
  Eigen::Matrix<double,3,6> jacobian_point_camera_pose;
  //Lie algebra of SO3
  jacobian_point_camera_pose.leftCols<3>() = world_SE3_camera.so3().matrix().transpose() * Sophus::SO3d::hat(point_world - world_SE3_camera.translation());
  //Translation part
  jacobian_point_camera_pose.rightCols<3>() = -world_SE3_camera.so3().matrix().transpose();
  jacobian_pose = jacobian_error_point_camera * jacobian_point_camera_pose;
  jacobians_[0] = jacobian_point_world;
  jacobians_[1] = jacobian_pose;
}
} //namespace SLAMSolver
