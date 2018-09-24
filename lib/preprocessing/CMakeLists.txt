# @file ReservoirSimulator/lib/preprocessing/CMakeLists.txt  -- CMake file for the External module

add_subdirectory(interpolators)

set(preprocessing_files ${interpolators_files} PARENT_SCOPE)

set(preprocessing_dir ${interpolators_dir} PARENT_SCOPE)

file(GLOB headers *.h)

file(GLOB sources *.cpp)

set(preprocessing_SRCS ${headers} ${sources} PARENT_SCOPE)

