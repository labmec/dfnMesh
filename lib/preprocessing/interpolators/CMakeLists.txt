########### next target ###############

file(GLOB interpolators_headers *.h)

file(GLOB interpolators_sources *.cpp)

set(interpolators_files ${interpolators_headers} ${interpolators_sources} PARENT_SCOPE)

set(interpolators_dir ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)