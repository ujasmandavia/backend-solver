#pragma once

#include <iostream>
#include <tuple>
#include "edge.hh"

namespace SLAMSolver{
/*
* @brief The variadic template base class for all error edges
* @tparams D: Dimensions of error
* @tparam VertexTypes: Types of vertices edges are connected to
*/
template<int D, typename ... VertexTypes>
class BaseEdge : public Edge{
public:
  /*
  * @brief constructor for BaseEdge class for all error edges
  * with error dimensions and number or vertices
  */
  BaseEdge();

  // @brief virtual distructor
  virtual ~BaseEdge() {};

  /// @brief Number of Vertices mentioned in the template list
  static int NumVertices = sizeof ...(VertexTypes);

  /*
  * @brief set vertex tuple in a pointer
  * @tparam Index: Index of the vertex yu want to see
  * @tparam VertexType: Type of vertex, this can be deducted, no need to specify
  * @param vertex_ptr: The vertex pointer want to set
  */
  template<size_t Index, typename VertexType>
  void set_vertex(std::shared_ptr<VertexType> vertex_ptr);

  /*
  * @brief Get vertex pointer in a tuple
  * @tparam Index: The index of the param we want to get
  * @return vertex pointer with the deducted type
  */
  template<size_t Index>
  auto get_vertex();

protected:
  // A tuple containing all the vertices pointer that the edge is connecting to
  std::tuple<std::shared_ptr<VertexType>...> tuple_vertices_ptr_;
}; //class BaseEdge

template<int D, typename ... VertexTypes>
BaseEdge<D,VertexTypes ...>::BaseEdge() : Edge(D,NumVertices) {}

template<int D, typename ... VertexTypes>
template<size_t Index, typename ... VertexType>
void BaseEdge<D,VertexTypes ...>::set_vertex(std::shared_ptr<VertexType> vertex_ptr){
  //Set template vertex pointer in tuple
  std::get<Index>(tuple_vertices_ptr_) = vertex_ptr;

  //Set Pointer for vertex interface pointer
  vertices_interface_ptr_[Index] = vertex_ptr;
}

template<int D,typename ... VertexTypes>
template<size_t Index>
auto BaseEdge<D,VertexTypes ... >::get_vertex(){
  //Get the template vertex pointer from tuple
  return std::get<Index>(tuple_vertices_ptr_);
}

} //namespace SLAMSolver
