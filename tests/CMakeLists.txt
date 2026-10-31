cmake_minimum_required(VERSION 3.18)

find_package(GTest REQUIRED)

add_executable(test_noBinaries test_noBinaries.cpp)

target_link_libraries(test_noBinaries
    base9
    yaml-cpp
    ${GSL_LIBRARIES}
    GTest::gtest
    GTest::gtest_main
)

target_include_directories(test_noBinaries PRIVATE
    ${CMAKE_SOURCE_DIR}/base9
    ${CMAKE_BINARY_DIR}/include
    ${GSL_INCLUDE_DIRS}
)

include(GoogleTest)
gtest_discover_tests(test_noBinaries)
