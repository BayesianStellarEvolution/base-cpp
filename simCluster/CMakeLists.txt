cmake_minimum_required(VERSION 2.8...3.26.0)

file(GLOB SIMCLUSTER_FILES *.cpp *.hpp)

include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_executable(simCluster ${SIMCLUSTER_FILES})

target_link_libraries(simCluster base9 m yaml-cpp)

install(TARGETS simCluster
        DESTINATION bin)
