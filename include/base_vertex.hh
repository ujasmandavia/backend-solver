#pragma once

#include "vertex.hh"

namespace SLAMSolver{

/*
 * @brief Template base class for vertex
 * @tparam D dimensions of minimal paraemtrization
 * @tparam ParamType Type for parameter's internal representation
*/
template<int D,typename ParamType>
class BaseVertex : public Vertex{
public:
  /// @brief Constructor for basevertex
  BaseVertex();
  /// @brief virtual destructor for BaseVertex;
  virtual ~BaseVertex();

  /// @return ParamType
  ParamType parameters() const;

  /// @brief set the parameters in the vertex
  void set_parameters(const ParamType& parameters);

protected:
  /// Storage of parameters in type of ParamType
  /// This is not exposed though vertex base class
  ParamType parameters_;
}; //class BaseVertex

template<int D, typename ParamType>
BaseVertex<D,ParamType>::BaseVertex() : Vertex(D) {}

template<int D, typename ParamType>
ParamType BaseVertex<D,ParamType>::parameters() const{
  return parameters_;
}

template<int D, typename ParamType>
void BaseVertex<D,ParamType>::set_parameters(const ParamType& parameters){
  parameters_ = parameters;
}
}  //namespace SLAMSolver
