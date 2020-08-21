#include "../vertex.hh"

namespace SLAMSolver{
Vertex::Vertex(const int minimal_dimension) : minimal_dimension_(minimal_dimension) {
  static unsigned long uuid_vertex = 0;
  id_ = uuid_vertex++;
}

Vertex::~Vertex() {}

} //namespace SLAMSolver
