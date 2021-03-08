########### next target ###############

include_directories(.)
include_directories(util)
add_subdirectory(util)

file(GLOB headers *.h *.tpp *Temp.cpp)
file(GLOB sources *.cpp)
# file(GLOB sources DFNNamespace.cpp DFNRib.cpp DFNFace.cpp DFNVolume.cpp DFNPolygon.cpp DFNFracture.cpp DFNMesh.cpp DFNPolyhedron.cpp)

list(FILTER sources EXCLUDE REGEX ".*Temp.cpp|main.*|quicktests.cpp")
# list(FILTER sources EXCLUDE REGEX ".*Temp.cpp")
# list(FILTER sources EXCLUDE REGEX "main.cpp")
# list(FILTER sources EXCLUDE REGEX "quicktests.cpp")


# Default Target
add_executable(dfnTest main.cpp  ${sources} ${headers})
target_link_libraries(dfnTest ${GMSH_LIB} pz)


# Target for performance tests
add_executable(perfTest main_perftests.cpp  ${sources} ${headers})
target_link_libraries(perfTest ${GMSH_LIB} pz)

# # # Quick Tests
# # add_executable(quickTests quicktests.cpp ${sources} ${headers})
# # target_link_libraries(quickTests ${GMSH_LIB} pz)

include(CTestDFN.cmake)


