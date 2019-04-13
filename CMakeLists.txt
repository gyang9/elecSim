cmake_minimum_required(VERSION 3.0)

project(ElecSim VERSION 1.0.0)

# Configure the dependencies
find_package(EDepSim)

add_executable(fineTrans fineTrans.C)
target_link_libraries(fineTrans PUBLIC EDepSim::edepsim_io)
