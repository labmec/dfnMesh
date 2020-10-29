########### next target ###############


file(GLOB headers *.h *.tpp *Temp.cpp)
file(GLOB sources DFNNamespace.cpp DFNRib.cpp DFNFace.cpp DFNVolume.cpp DFNPolygon.cpp DFNFracture.cpp DFNMesh.cpp DFNPolyhedron.cpp)

# Default Target
add_executable(dfnTest main.cpp  ${sources} ${headers})
target_link_libraries(dfnTest ${GMSH_LIB} pz)

# # Quick Tests
# add_executable(quickTests quicktests.cpp)
# target_link_libraries(quickTests ${GMSH_LIB} pz)


#test