########### next target ###############


add_executable(LaberintoTest LaberintoTest.cpp ${lib_files})

add_executable(InterpolatorTest InterpolatorTest.cpp ${lib_files})

add_executable(3RelPermTest Interpolator3RelPermTest.cpp ${lib_files})


target_include_directories(LaberintoTest PUBLIC ${lib_dir})
target_include_directories(InterpolatorTest PUBLIC ${lib_dir})
target_include_directories(3RelPermTest PUBLIC ${lib_dir})

target_link_libraries(LaberintoTest pz ${OpenCV_LIBS})

target_link_libraries(InterpolatorTest pz ${OpenCV_LIBS})
target_link_libraries(3RelPermTest pz ${OpenCV_LIBS})

