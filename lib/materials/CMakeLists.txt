########### next target ###############

file(GLOB materials_headers *.h)

file(GLOB materials_sources *.cpp)

set(materials_files ${materials_headers} ${materials_sources} PARENT_SCOPE)

set(materials_dir ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)