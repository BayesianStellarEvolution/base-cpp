cmake_minimum_required(VERSION 3.18...3.30)

file(GLOB BASE9_FILES *.cpp *.hpp MsRgbModels/*.cpp MsRgbModels/*.hpp WdCoolingModels/*.cpp WdCoolingModels/*.hpp WdAtmosphereModels/*.cpp WdAtmosphereModels/*.hpp IO/*.cpp IO/*.hpp)

include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_library(base9 STATIC ${BASE9_FILES})

set_target_properties(base9 PROPERTIES POSITION_INDEPENDENT_CODE true)
target_link_libraries(base9 dl sleef)
target_include_directories(base9 PRIVATE ${CMAKE_BINARY_DIR}/include)

install(TARGETS base9
        DESTINATION lib)
