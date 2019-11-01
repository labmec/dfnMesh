########### next target ###############


file(GLOB header *.h)
file(GLOB sources *.cpp)

add_executable(dfnTest main.cpp  ${sources} ${header})

#target_include_directories(RibTest)

target_link_libraries(dfnTest ${GMSH_LIB} pz)