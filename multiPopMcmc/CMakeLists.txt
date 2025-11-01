cmake_minimum_required(VERSION 2.8...3.26.0)

file(GLOB MULTIPOPMCMC_FILES *.cpp *.hpp)

include_directories( ${YAML_CPP_INCLUDE_DIRS} "/usr/local/include" )

add_executable(multiPopMcmc ${MULTIPOPMCMC_FILES})

target_link_libraries(multiPopMcmc base9 m GSL::gsl GSL::gslcblas yaml-cpp)

install(TARGETS multiPopMcmc
        DESTINATION bin)
