########### next target ###############


file(GLOB header *.h)
file(GLOB sources *.cpp)

add_executable(RibTest RibFrac.cpp  ${sources} ${header})

#target_include_directories(RibTest)

target_link_libraries(RibTest pz ${OpenCV_LIBS})