cmake_minimum_required (VERSION 2.8)

file(GLOB_RECURSE CFILES "${PROJECT_SOURCE_DIR}/examples/*.cpp")
SET_SOURCE_FILES_PROPERTIES(${CFILES} PROPERTIES LANGUAGE CXX COMPILE_FLAGS "-std=c++11")

add_executable(rbruck ./radix_r_bruck_example.cpp)
add_executable(gerNum ./generate_ncomm_ndb.cpp)
add_executable(twolayer ./twolayer_example.cpp)

target_link_libraries(rbruck rbruckLib)
target_link_libraries(gerNum rbruckLib)
target_link_libraries(twolayer rbruckLib)