#include "../pose3d_vertex.hh"

namespace SLAMSolver{

Pose3DVertex::Pose3DVertex(const bool combined) : combined_(combined) {}

void Pose3DVertex::add(Eigen::VectorXd& delta){
  if(combined_){
    Eigen::Matrix<6,1> delta_ = delta.template head<6>();
    Sophus::SE3 delta_SE3 = Sophus::SE3::exp(delta_);
    parameters_ = delta_SE3 * parameters_;
  }else{
    //delta(0) delta(1) delta(2) is so3, lie algebra for rotation
    //delta(3) delta(4) delta(5) is pure translation if combined_ == false
    //construct delta SO3
    Sophus::SO3 delta_SO3(delta(0), delta(1), delta(2));
    //delta translation
    Eigen::Vector3d delta_translation(delta(3), delta(4), delta(5));
    //Left multiplication model, all edge connected to this Pose3DVertex also need to
    //use left multiplication model
    Sophus::SO3 new_SO3 = delta_SO3 * parameters_.so3();
    Eigen::Vector3d new_translation = delta_translation + parameters_.translation();
    parameters_ = Sophus::SE3(new_SO3, new_translation);
  }
}

} //namespace SLAMSolver
