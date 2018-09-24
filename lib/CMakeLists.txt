

add_subdirectory(preprocessing)
add_subdirectory(materials)

set(lib_files ${materials_files} ${preprocessing_files} PARENT_SCOPE)
set(lib_dir ${materials_dir} ${preprocessing_dir} PARENT_SCOPE)




source_group(PreProcessing FILES ${preprocessing_SRCS})

set(reservoirlibsources
	${preprocessing_SRCS}

)

add_library(reservoirlib ${reservoirlibsources})

