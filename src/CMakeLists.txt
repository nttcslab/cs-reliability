cmake_minimum_required(VERSION 3.8)
project(reliability)

add_subdirectory(./mylib)
add_executable(main main.cpp)
add_executable(main_single main_single.cpp)
add_executable(tdzdd tdzdd.cpp)
add_executable(tdzdd_single tdzdd_single.cpp)

target_link_libraries(main reliabilitylib)
target_link_libraries(main_single reliabilitylib)
target_link_libraries(tdzdd reliabilitylib)
target_link_libraries(tdzdd_single reliabilitylib)
