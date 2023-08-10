cmake_minimum_required(VERSION 3.20)
project(TBB-test)

find_package(OpenMP REQUIRED)
add_executable(RLA_parallel RLA_parallel.cpp)

target_link_libraries(RLA_parallel PRIVATE OpenMP::OpenMP_CXX)
