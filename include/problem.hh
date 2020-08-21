#pragma once

#include <iostream>
#include <unordered_map>
#include <map>
#include <memory>

#include "eigen_types.hh"
#include "edge.hh"
#include "vertex.hh"

namespace SLAMSolver{
/*
* @brief Problem class contaning all the vertices connected by edges
*/
class Problem{
public:
  using IDType = unsigned long;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  // @brief constructor
  Problem();

  // @brief destructor
  ~Problem();

  /*
  * @brief add_vertex : add vertex to the problem
  * @param vertex: vertex pointer
  * @return True if the vertex is added
  */
  bool add_vertex(std::shared_ptr<Vertex> vertex);

  /*
  * @brief remove_vertex : remove vertex from the problem
  * @param vertex: vertex pointer
  * @return True if the vertex is added
  */
  bool remove_vertex(std::shared_ptr<Vertex> vertex);

  /*
  * @brief add_edge : add edge to the problem
  * @param edge: edge pointer
  * @return True if the edge is added
  */
  bool add_edge(std::shared_ptr<Edge> edge);

  /*
  * @brief remove_edge : remove edge from the problem
  * @param edge: edge pointer
  * @return True if the edge is removed
  */
  bool remove_edge(std::shared_ptr<Edge> edge);

public:
  /*
  * @brief cost of the whole problem by summing cost from all edges
  * return accumulated total cost of the problem
  */
  double computer_cost();

  //Hash map for Vertex ID's to vertex pointer
  std::map<IDType,std::shared_ptr<Vertex>> vertex_id_map_vertex_ptr_;

  // Hash map from Edge's ID to Edge pointer
  std::map<IDType, std::shared_ptr<Edge>> edge_id_map_edge_ptr_;

  // Hash map from Vertex's ID to Edges that it connected to
  std::unordered_multimap<IDType, std::shared_ptr<Edge>>  vertex_id_map_edges_ptr_;
}
} //namespace SLAMSolver
