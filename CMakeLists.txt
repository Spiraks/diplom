cmake_minimum_required(VERSION 3.0.0)
project(FDTD VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake) 
find_package(SFML REQUIRED COMPONENTS graphics window system) 

include_directories(${SFML_INCLUDE_DIR}) 

add_executable(3D src/3D.cpp src/FDTD.cpp) 
target_link_libraries(3D ${SFML_LIBRARIES} ${SFML_DEPENDENCIES}) 
target_link_libraries(3D sfml-graphics)
target_compile_options(3D PUBLIC -O2)


