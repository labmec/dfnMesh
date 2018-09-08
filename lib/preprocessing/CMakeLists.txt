########### next target ###############

add_subdirectory(interpolators)

set(preprocessing_files ${interpolators_files} PARENT_SCOPE)

set(preprocessing_dir ${interpolators_dir} PARENT_SCOPE)