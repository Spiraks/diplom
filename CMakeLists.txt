cmake_minimum_required(VERSION 3.0.0)
project(FDTD VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake) 
# find_package(SFML REQUIRED COMPONENTS graphics window system) 
# include_directories(${SFML_INCLUDE_DIR})
# find_package( OpenCV REQUIRED )
find_package(OpenMP REQUIRED)
include_directories(${OpenMP_INCLUDE_DIRS})
# Установка флагов компилятора
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

# Добавление флагов для offloading на GPU
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -foffload=nvptx-none")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-offload")

add_executable(3D src/3D.cpp src/FDTD.cpp src/Mur.cpp)
if(OpenMP_CXX_FOUND)
    target_link_libraries(3D PUBLIC OpenMP::OpenMP_CXX)
endif()
add_compile_options(-fcf-protection=none)
target_compile_options(3D PUBLIC -O2)
