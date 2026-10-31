cmake_minimum_required(VERSION 3.18)

find_package(benchmark REQUIRED)

add_executable(bench_noBinaries bench_noBinaries.cpp)

target_link_libraries(bench_noBinaries
    base9
    yaml-cpp
    ${GSL_LIBRARIES}
    benchmark::benchmark
)

target_include_directories(bench_noBinaries PRIVATE
    ${CMAKE_SOURCE_DIR}/base9
    ${CMAKE_BINARY_DIR}/include
    ${GSL_INCLUDE_DIRS}
)
