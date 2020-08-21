#pragma once


#include "solver/base_edge.hh"
#include "solver/pose3d_vertex.hh"
#include "solver/point3d_vertex.hh"

namespace SLAMSolver {
/*!
 * @brief Reprojection Error Edge connecting one 3D keypoint and one camera pose
 */
class ReprojectionEdge : public BaseEdge<2, Point3DVertex, Pose3DVertex>{
public:
  /// @brief Constructor, initialize observation
  ReprojectionEdge(const Eigen::Vector2d& uv);
  /// @brief Default destructor
  ~ReprojectionEdge() = default;
  /// @brief override interface, compute reprojection error
  void compute_errors() override;
  /// @brief override interface, compute jacobian matrices for point and pose
  /// @note Use left perturbation in Pose3DVertex, so keep use left perturbation in computing jacobian
  void compute_jacobians() override;
private:
  /// This is the point on z_camera=1 plane, not the pixel coordinate in image
  /// Use this representation to save computation of transforming to pixel coordinate
  Eigen::Vector2d uv_;
};
} //namespace SLAMSolver
