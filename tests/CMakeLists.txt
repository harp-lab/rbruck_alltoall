cmake_minimum_required (VERSION 3.5)

file(GLOB_RECURSE CFILES "${PROJECT_SOURCE_DIR}/test/*.cpp")
SET_SOURCE_FILES_PROPERTIES(${CFILES} PROPERTIES LANGUAGE CXX COMPILE_FLAGS "-std=c++11")

add_executable(pingpong ./pingpongcheck.cpp)
add_executable(inverse ./inverseComm.cpp)

target_link_libraries(pingpong rbruckLib)
target_link_libraries(inverse rbruckLib)