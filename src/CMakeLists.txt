########### next target ###############

include_directories("." "util")
add_subdirectory(util)

file(GLOB headers *.h *.tpp *Temp.cpp)
file(GLOB sources *.cpp)
# file(GLOB sources DFNNamespace.cpp DFNRib.cpp DFNFace.cpp DFNVolume.cpp DFNPolygon.cpp DFNFracture.cpp DFNMesh.cpp DFNPolyhedron.cpp)

# Using regex filter to clear targets from ${sources}. "|" is the 'or' operator in regex
list(FILTER sources EXCLUDE REGEX ".*Temp.cpp|main.*|quicktests.cpp")
# # We could filter them separately too
# list(FILTER sources EXCLUDE REGEX ".*Temp.cpp")
# list(FILTER sources EXCLUDE REGEX "main.cpp")
# list(FILTER sources EXCLUDE REGEX "quicktests.cpp")

# Create a library for the project so separate targets won't require rebuilding of project source code
add_library(dfnlib ${headers} ${sources})
target_link_libraries(dfnlib ${GMSH_LIB} pz)

# Default Target
add_executable(dfnTest main.cpp)
target_link_libraries(dfnTest dfnlib)


# Target for performance tests
add_executable(perfTest main_perftests.cpp)
target_link_libraries(perfTest dfnlib)


# A cleaner target to summarize main methods
add_executable(SimpleDFN main_simple.cpp)
target_link_libraries(SimpleDFN dfnlib)

# # Quick Tests
# add_executable(quickTests quicktests.cpp)
# target_link_libraries(quickTests dfnlib)


# CTests are in a separate file
include(CTestDFN.cmake)


