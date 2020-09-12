project(backend_solver)
cmake_minimum_required(VERSION 3.2)

set(DEFAULT_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-STD=C++14 -g -Wall")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wno-reorder" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -DNDEBUG -Wno-reorder -O2" CACHE STRING "" FORCE)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(Eigen REQUIRED)
include_directories(${EIGEN_INCLUDE_DIR})


# sophus
include_directories(${PROJECT_SOURCE_DIR}/thirdparty/Sophus)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_library(SLAMSolver SHARED
            include/edge.hh
            src/edge.cc
            include/base_edge.hh
            include/vertex.hh
            include/base_vertex.hh
            include/pose3d_vertex.hh
            src/pose3d_vertex.cc
            include/problem.hh
            src/problem.cc
            src/relative_pose3d_edge.cc
            include/relative_pose3d_edge.hh
            include/base_solver.hh
            src/plain_solver.cc
            include/plain_solver.hh
            src/base_solver.cc
            src/minimize_levenberg_marquardt.cc
            include/minimize_levenberg_marquardt.hh
            src/sparse_schur_solver.cc
            include/sparse_schur_solver.hh
            src/sparse_cholesky_solver.cc
            include/sparse_choleksy_solver.hh
            src/vertex.cc
            src/reprojection_edge.cc
            include/reprojection_edge.hh
            )
                                                                                                                                                  sparse_add_subdirectory(app)
