cmake_minimum_required(VERSION 2.8...3.26.0)

file(GLOB SCATTERCLUSTER_FILES *.cpp *.hpp)

include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_executable(scatterCluster ${SCATTERCLUSTER_FILES})

target_link_libraries(scatterCluster base9 m yaml-cpp)

install(TARGETS scatterCluster
        DESTINATION bin)
