cmake_minimum_required(VERSION 2.8...3.26.0)

file(GLOB SAMPLEWDMASS_FILES *.cpp *.hpp)

include_directories( ${YAML_CPP_INCLUDE_DIRS} "/usr/local/include" )

add_executable(sampleMass ${SAMPLEWDMASS_FILES})

target_link_libraries(sampleMass base9 m yaml-cpp ${GSL_LOC} ${GSLCBLAS_LOC})

install(TARGETS sampleMass
        DESTINATION bin)
