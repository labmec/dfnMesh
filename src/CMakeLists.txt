########### next target ###############

# include_directories("." "util")
add_subdirectory(util)

file(GLOB_RECURSE headers *.h *.tpp *Temp.cpp)
file(GLOB_RECURSE sources *.cpp)

# Using regex filter to clear targets from ${sources}. "|" is the 'or' operator in regex
list(FILTER sources EXCLUDE REGEX ".*Temp.cpp|main.*|quicktests.cpp")
# # We could filter them separately too
# list(FILTER sources EXCLUDE REGEX ".*Temp.cpp")
# list(FILTER sources EXCLUDE REGEX "main.cpp")
# list(FILTER sources EXCLUDE REGEX "quicktests.cpp")

# include_directories(PUBLIC ${PZ_INCLUDE_DIRS})
set(DFN_INCLUDE_DIRS ${PZ_INCLUDE_DIRS} ${GMSH_INC} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/util)

# Create a library for the project so separate targets won't require rebuilding of project source code
add_library(dfnlib ${headers} ${sources})
target_link_libraries(dfnlib PUBLIC ${GMSH_LIB} NeoPZ::pz)
# target_include_directories(dfnlib PUBLIC ${DFN_INCLUDE_DIRS} )

# This package is here because iMRS is the only one that uses it
find_package(libInterpolate REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../libinterpolate_install ${CMAKE_SOURCE_DIR}/libinterpolate_install)
link_libraries(dfnlib libInterpolate::Interpolate)

target_include_directories(dfnlib PUBLIC ${PZ_INCLUDE_DIRS})

target_include_directories(dfnlib PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )

target_include_directories(dfnlib PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/util>
    )

target_include_directories(dfnlib PUBLIC 
    $<BUILD_INTERFACE:${GMSH_INC}>
    $<INSTALL_INTERFACE:${GMSH_INC}>
    )

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
target_sources(dfnlib PRIVATE ${headers} ${sources})



# Default Target
# Main Go back is the default now (may 2022)
# add_executable(dfnTest main.cpp)
# target_link_libraries(dfnTest PRIVATE dfnlib)
# A cleaner more efficient target for animations
# add_executable(animations main_animations.cpp)
# target_link_libraries(animations PRIVATE dfnlib)

# add_pz_target(
#   NAME dfnTest
#   SOURCES ${headers} "main.cpp" ${sources}
#   REQUIRED PZ_USING_MKL)
# target_link_libraries(dfnTest PRIVATE ${GMSH_LIB})

# Target for testing rollback
add_executable(dfnGoBack main_goback.cpp)
target_link_libraries(dfnGoBack PRIVATE dfnlib)
# A cleaner more efficient target for animations
# add_executable(animations main_animations.cpp)
# target_link_libraries(animations PRIVATE dfnlib)


# Target for performance tests
# add_executable(perfTest main_perftests.cpp)
# target_link_libraries(perfTest dfnlib)


# A cleaner target to summarize main methods
# add_executable(SimpleDFN main_simple.cpp)
# target_link_libraries(SimpleDFN PRIVATE dfnlib)

# # Quick Tests
# add_executable(quickTests quicktests.cpp)
# target_link_libraries(quickTests dfnlib)


# CTests are in a separate file
include(CTestDFN.cmake)

