cmake_minimum_required(VERSION 2.8)
add_executable(MAC Main.cpp Field.cpp)

project(MAC)

find_package(OpenCV REQUIRED)
find_package(GLUT)
find_package(OpenGL)
include_directories(.)
add_definitions(-std=c++11)
set(CMAKE_CXX_FLAGS "-Wno-deprecated-declarations -O2")
if(OpenCV_FOUND)
  include_directories(${OpenCV_INCLUDE_DIRS})

  target_link_libraries(MAC ${OpenCV_LIBS})
endif(OpenCV_FOUND)
target_link_libraries(MAC ${GLUT_LIBRARY} ${OPENGL_LIBRARY})
