cmake_minimum_required(VERSION 3.17)
project(network)

add_definitions(-DLOCAL)
set(CMAKE_CXX_STANDARD 17)

add_executable(network2 network_improved.cpp)
