cmake_minimum_required(VERSION 3.1)
project(libigl-pupa)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 ON)

find_package(LIBIGL REQUIRED QUIET)


include_directories(extern/voronoi/src)
include_directories(extern)

add_executable(remeshing main.cpp weighted_cvd.cpp extern/clipper.cpp )
target_link_libraries(remeshing igl::core igl::opengl igl::opengl_glfw igl::triangle )
