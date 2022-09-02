cmake_minimum_required (VERSION 2.8)

file(GLOB_RECURSE CFILES "${PROJECT_SOURCE_DIR}/*.cpp")
SET_SOURCE_FILES_PROPERTIES(${CFILES} PROPERTIES LANGUAGE CXX )

aux_source_directory(. RBURK_SOURCES)
add_library(rbruckLib ${RBURK_SOURCES})

target_link_libraries(rbruckLib ${MPI_C_LIBRARIES})
