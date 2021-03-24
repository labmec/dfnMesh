########### next target ###############


file(GLOB headers *.h *.hpp *.tpp *Temp.cpp)
file(GLOB sources *.cpp)

include_directories(.)

# target_include_directories(dfnTest PUBLIC 
#                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#                         #    $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/External/FAD>
#                           )