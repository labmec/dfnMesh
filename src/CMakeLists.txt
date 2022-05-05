########### next target ###############

# include_directories("." "util")
add_subdirectory(util)

file(GLOB_RECURSE headers *.h *.tpp *Temp.cpp)
file(GLOB_RECURSE sources *.cpp)

# Using regex filter to clear targets from ${sources}. "|" is the 'or' operator in regex
list(FILTER sources EXCLUDE REGEX ".*Temp.cpp|main.*|quicktests.cpp")
# # We could filter them separately too
# list(FILTER sources EXCLUDE REGEX ".*Temp.cpp")
# list(FILTER sources EXCLUDE REGEX "main.cpp")
# list(FILTER sources EXCLUDE REGEX "quicktests.cpp")

# include_directories(PUBLIC ${PZ_INCLUDE_DIRS})
set(DFN_INCLUDE_DIRS ${PZ_INCLUDE_DIRS} ${GMSH_INC} "." "util")

# Create a library for the project so separate targets won't require rebuilding of project source code
add_library(dfnlib ${headers} ${sources})
target_link_libraries(dfnlib PUBLIC ${GMSH_LIB} NeoPZ::pz)
target_include_directories(dfnlib PUBLIC ${DFN_INCLUDE_DIRS})



# Default Target
add_executable(dfnTest main.cpp)
target_link_libraries(dfnTest PRIVATE dfnlib)
# A cleaner more efficient target for animations
add_executable(animations main_animations.cpp)
target_link_libraries(animations PRIVATE dfnlib)
# add_pz_target(
#   NAME dfnTest
#   SOURCES ${headers} "main.cpp" ${sources}
#   REQUIRED PZ_USING_MKL)
# target_link_libraries(dfnTest PRIVATE ${GMSH_LIB})

# Target for testing rollback
add_executable(dfnGoBack main_goback.cpp)
target_link_libraries(dfnGoBack PRIVATE dfnlib)
# A cleaner more efficient target for animations
# add_executable(animations main_animations.cpp)
# target_link_libraries(animations PRIVATE dfnlib)


# Target for performance tests
# add_executable(perfTest main_perftests.cpp)
# target_link_libraries(perfTest dfnlib)


# A cleaner target to summarize main methods
# add_executable(SimpleDFN main_simple.cpp)
# target_link_libraries(SimpleDFN PRIVATE dfnlib)

# # Quick Tests
# add_executable(quickTests quicktests.cpp)
# target_link_libraries(quickTests dfnlib)


# CTests are in a separate file
include(CTestDFN.cmake)

# Create folder with input files
set(INPUTMESHESDIR ${CMAKE_SOURCE_DIR}/examples)
set (CONFIG_H "${CMAKE_SOURCE_DIR}/src/dfn_config.h")
configure_file ("${CMAKE_SOURCE_DIR}/src/dfn_config.h.in" ${CONFIG_H})

