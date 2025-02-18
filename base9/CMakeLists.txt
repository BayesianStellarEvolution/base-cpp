cmake_minimum_required(VERSION 2.8...3.26.0)

if ("${CMAKE_HOST_SYSTEM_PROCESSOR}" MATCHES "arm64" OR
    "${CMAKE_HOST_SYSTEM_PROCESSOR}" MATCHES "aarch64")
  file(GLOB BASE9_FILES *.cpp *.hpp MsRgbModels/*.cpp MsRgbModels/*.hpp WdCoolingModels/*.cpp WdCoolingModels/*.hpp WdAtmosphereModels/*.cpp WdAtmosphereModels/*.hpp IO/*.cpp IO/*.hpp)
else ()
  file(GLOB BASE9_FILES *.cpp *.hpp MsRgbModels/*.cpp MsRgbModels/*.hpp WdCoolingModels/*.cpp WdCoolingModels/*.hpp WdAtmosphereModels/*.cpp WdAtmosphereModels/*.hpp IO/*.cpp IO/*.hpp sleef/*.c sleef/*.h)
endif ()

include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_library(base9 STATIC ${BASE9_FILES})

set_target_properties(base9 PROPERTIES POSITION_INDEPENDENT_CODE true)
target_link_libraries(base9 dl)

install(TARGETS base9
        DESTINATION lib)
