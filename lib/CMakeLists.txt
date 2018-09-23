########### next target ###############



add_subdirectory(preprocessing)
add_subdirectory(materials)

set(lib_files ${materials_files} ${preprocessing_files} PARENT_SCOPE)
set(lib_dir ${materials_dir} ${preprocessing_dir} PARENT_SCOPE)



