cmake_minimum_required(VERSION 3.8)
project(reliabilitylibs)

set(SRCS graph.cpp)

add_library(reliabilitylib STATIC ${SRCS})

target_compile_features(reliabilitylib PUBLIC cxx_std_11)
#target_include_directories(reliabilitylib INTERFACE .)