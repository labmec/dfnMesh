########### next target ###############


file(GLOB headers *.h)
file(GLOB sources DFNRib.cpp DFNFace.cpp DFNVolume.cpp DFNPolygon.cpp DFNFracture.cpp DFNMesh.cpp)

# Default Target
add_executable(dfnTest main.cpp  ${sources} ${headers})
target_link_libraries(dfnTest ${GMSH_LIB} pz)

# Quick Tests
add_executable(quickTests quicktests.cpp)
target_link_libraries(quickTests ${GMSH_LIB} pz)