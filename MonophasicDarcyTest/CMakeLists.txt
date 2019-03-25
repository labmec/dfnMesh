########### next target ###############


add_executable(MonophasicDarcyTest MonophasicDarcyTest.cpp ${lib_files})



target_include_directories(MonophasicDarcyTest PUBLIC ${lib_dir})


target_link_libraries(MonophasicDarcyTest pz ${OpenCV_LIBS})



