#pragma once
#include "../base_vertex.hh"

namespace SLAMSolver{

/*
 * @brief graph optimization vertex representing 3D point position
*/
class Point3DVertex : public BaseVertex<3,Eigen::Vector3d>{
public:
  /// @brief Defauly constructor
  Point3DVertex() = default;
  /// @brief Default destructor
  ~Point3DVertex() = default;
  /// @brief trivial add
  void add(const Eigen::VectorXd& delta) override{
    parameters_(0) += delta(0);
    parameters_(1) += delta(1);
    parameters_(2) += delta(2);
  }
};

} //namespace SLAMSolver
