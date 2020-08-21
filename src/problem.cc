#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "problem.hh"

namespace SLAMSolver{

Problem::Problem() {}

Problem::~Problem() {}

bool Problem::add_vertex(std::shared_ptr<Vertex> vertex){
  if(vertex_id_map_vertex_ptr_.find(vertex->id()) != vertex_id_map_vertex_ptr_.end()){
    return false;
  }else{
    vertex_id_map_vertex_ptr_[vertex->id()] = vertex;
    return true;
  }
}

bool Problem::remove_vertex(std::shared_ptr<Vertex> vertex){
  //remove vertex and all the edges connected to it
  vertex_id_map_vertex_ptr_.erase(vertex->id());
  auto edges = vertex_id_map_edges_ptr_.equal_range(vertex->id());
  for(auto it = edges.first; it != edges.second; it++){
    remove_edge(it->second);
  }
}

bool Problem::add_edge(std::shared_ptr<Edge> edge){
  //Add edge interface pointer to hash map with its IDs as key
  if(edge_id_map_edge_ptr_.find(edge->id()) != edge_id_map_edge_ptr_.end()){
    return false;
  }else{
    edge_id_map_edge_ptr_[edge->id()] = edge;
    for(int i=0; i<edge->num_vertices(); i++){
      auto vertices_interface_ptr = edge->get_vertex_interface(i);
      add_vertex(vertex_interface_ptr);
      //Associate the vertex id with edge
      vertex_id_map_edges_ptr_.insert(std::pair<IDType,std::shared_ptr<Edge>>(vertex_interface_ptr->id(),edge));
    }
  }
}

bool Problem::remove_edge(std::shared_ptr<Edge> edge) {
  if(edge_id_map_edge_ptr_.find(edge->id()) == edge_id_map_edge_ptr_.end()){
    std::cerr<<"Edge want to delete not in Problem"<<std::endl;
  } else{
    for(int i = 0; i < edge->num_vertices(); i++){
      auto v = edge->get_vertex_interface(i);
      if(vertex_id_map_edges_ptr_.count(v->id()) == 1){
        //this edge is the only edge that the vertex connected to
        //delete the vertex
        remove_vertex(v);
      }
    }
    edge_id_map_edge_ptr_.erase(edge->id());
  }
}

double Problem::computer_cost(){
  //Compute total cost
  double total_cost = 0;
  for (auto &edge: edge_id_map_edge_ptr_) {
    //Compute errors and error function's jacobians
    edge.second->compute_errors();
    total_cost += edge.second->chi2();
  }

  return total_cost;
}

} //namespace SLAMSolver
