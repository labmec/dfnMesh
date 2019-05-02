########### next target ###############


add_executable(LaplacianJorge LaplacianJorge.cpp ${lib_files} ${SOURCE} ${HEADERS})

target_include_directories(LaplacianJorge PUBLIC ${lib_dir})

target_link_libraries(LaplacianJorge pz ${OpenCV_LIBS})