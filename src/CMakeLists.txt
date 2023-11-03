cmake_minimum_required(VERSION 3.20)
project(ugap_simulation)

set(CMAKE_CXX_STANDARD 14)

add_executable(ugap_simulation
               main.cpp
               Voxel.cpp
               Voxel.h
               common.h
        )