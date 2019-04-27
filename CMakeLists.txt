cmake_minimum_required(VERSION 3.0)

project(ElecSim VERSION 1.0.0)

# Configure the dependencies
find_package(EDepSim)

add_executable(fineTrans fineTrans.C)
target_link_libraries(fineTrans PUBLIC EDepSim::edepsim_io)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++11")

set(EIGEN_INC_DIR /home/guang/work/Pandora/Eigen3)
include_directories(${EIGEN_INC_DIR})

#set(SOURCE_FILES kalman.cpp kalman-test.cpp)
#add_executable(kalman-test ${SOURCE_FILES})

