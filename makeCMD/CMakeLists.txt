cmake_minimum_required(VERSION 2.8...3.26.0)

file(GLOB MAKECMD_FILES makeCMD.cpp) # *.cpp *.hpp)
file(GLOB MAKECMDS_FILES makeCMDs.cpp) # *.cpp *.hpp)
file(GLOB MAKECMDWD_FILES makeCMDwd.cpp) # *.cpp *.hpp)

include_directories( ${YAML_CPP_INCLUDE_DIRS} )

add_executable(makeCMD ${MAKECMD_FILES})
add_executable(makeCMDs ${MAKECMDS_FILES})
add_executable(makeCMDwd ${MAKECMDWD_FILES})

target_link_libraries(makeCMD base9 m yaml-cpp)
target_link_libraries(makeCMDs base9 m yaml-cpp)
target_link_libraries(makeCMDwd base9 m yaml-cpp)

install(TARGETS makeCMD makeCMDs makeCMDwd
        DESTINATION bin)
