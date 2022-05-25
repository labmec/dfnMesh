# CMakeLists root for DFNMesh project

cmake_minimum_required (VERSION 3.6)
# message("")


set(CMAKE_CONFIGURATION_TYPES Debug Release CACHE STRING  INTERNAL FORCE )

project (DFNMesh)

# Adds definition of macro PROJECT_ROOT, which contains a string of the path to project directory
add_compile_definitions("PROJECT_ROOT=\"${PROJECT_SOURCE_DIR}\"")

# Install paths
include(GNUInstallDirs)

set(DFNLIB_INSTALL_PREFIX "dfnlib")
set(CMAKE_INSTALL_LIBDIR ${DFNLIB_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_INCLUDEDIR ${DFNLIB_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
set(CMAKE_INSTALL_BINDIR ${DFNLIB_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

file(RELATIVE_PATH relDir
${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
)


IF (WIN32)
	# Prevents timespec redefinition problem with Visual Studio 2015
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")
ELSE (WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")#this flag is for preventing errors when casting from double to float
ENDIF()


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Mac OS X specific code
IF(APPLE)
   add_definitions(-DMACOSX)
   set(MACOSX MACOSX)
   FIND_LIBRARY(ACCELERATE_LIB Accelerate)
   link_libraries(${ACCELERATE_LIB})
#   SET(OperatingSystem "Mac OS X")
ENDIF(APPLE)


# NeoPZ
IF (WIN32)
	find_package( NeoPZ PATHS "C:/Arquivos de Programas/PZ" REQUIRED NO_DEFAULT_PATH)
ELSE()
	find_package( NeoPZ REQUIRED HINTS "../neopz_install/" "../lib/neopz-develop")
ENDIF (WIN32)

# include_directories(${PZ_INCLUDE_DIRS})
MESSAGE( STATUS "PZ_INCLUDE_DIRS:         " ${PZ_INCLUDE_DIRS} )
MESSAGE( STATUS "PZ_DIR:                  " ${NeoPZ_DIR} )


# Location to search files or libraries
set(SEARCH_DIRS ENV PATH /usr/local/opt/ /usr/local/opt/qt5/lib/cmake/ /usr/include /usr/local/include /usr/local/include/boost /usr/local/include/log4cxx /usr/local/include/tbb /usr/lib /usr/lib64 /usr/local/lib ./. ./.. ./../libs ./../libs/fad ./../libs/metis ./../libs/metis/include ./../libs/metis/lib ./../externallibs ./../externallibs/fad ./../externallibs/pthread ./../externallibs/metis ./../externallibs/metis/include ./../externallibs/metis/lib ./..metis ./../metis/include ./../metis/lib ./externallibs/lib ./externallibs/include ./../externallibs/lib ./../externallibs/include)

#enabling pthread library
# find_path(PTHREAD_INCLUDE NAMES PThread.h pthread.h PATHS ${SEARCH_DIRS} D:/Projetos/externallibs/include )
# find_library(PTHREAD_LIB NAMES libpthread.dylib pthread.lib pthreadVC2.lib libPThread.so libPThread.a libpthread.so libpthread.a PATHS ${SEARCH_DIRS} D:/Projetos/externallibs/lib  ${PROJECT_SOURCE_DIR}/../pthread/lib)
# if(PTHREAD_INCLUDE-NOTFOUND)
# 	set(PTHREAD_INCLUDE "" CACHE PATH "Directory where pthread.h can be found")
# endif()
# if(PTHREAD_LIB-NOTFOUND)
# 	set(PTHREAD_LIBexit "" CACHE PATH "Directory where the pthread library can be found")
# endif()
# include_directories(${PTHREAD_INCLUDE})

#message(${PTHREAD_INCLUDE})

# Environment variables
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	set(PZDEBUG True)
	add_definitions(-DPZDEBUG)
else() # Debug
	# file(MAKE_DIRECTORY "release")
endif()
if (CMAKE_BUILD_TYPE STREQUAL "Release")
	set(NODEBUG NODEBUG)
	set(ZERO_INTERNAL_RESIDU ZERO_INTERNAL_RESIDU)
        set(MAKEINTERNAL MAKEINTERNAL)
	add_definitions(-DNODEBUG)
endif()

# Windows porting
if(WIN32)
	#add a define saying it is a VC compiler
	set(VC "#Will add a define VC on pz_config.h")
	add_definitions(-DVC)
	#define _USE_MATH_DEFINES for example: M_PI constant definitions.
	add_definitions(-D_USE_MATH_DEFINES)
	#disabling VC warnings
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	#disabling All warnings
	#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W0")
	#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /W0")
	#disabling STL warning for std::_Vector_alloc when LOG4CXX is enabled
	if (USING_LOG4CXX)
	    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /wd4251")
	    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /wd4251")
	endif()
	#define use of pthread static lib.
	add_definitions(-DPTW32_STATIC_LIB)
    #enabling /bigobj
    add_definitions("/bigobj")

    FOREACH(FLAG_TYPE EXE MODULE SHARED)
        STRING (REPLACE "INCREMENTAL:YES" "INCREMENTAL" FLAG_TMP "${CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG}")
        STRING (REPLACE "INCREMENTAL:NO" "INCREMENTAL" FLAG_TMP ${FLAG_TMP})
        STRING (REPLACE "INCREMENTAL" "INCREMENTAL:NO" FLAG_TMP ${FLAG_TMP})
        STRING (REPLACE "/EDITANDCONTINUE" "" FLAG_TMP ${FLAG_TMP})
        SET(CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG "${FLAG_TMP}" CACHE STRING "Overriding default debug ${FLAG_TYPE} linker flags." FORCE)
        MARK_AS_ADVANCED (CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG)
    ENDFOREACH ()

endif()



# linking with GMSH
find_library(GMSH_LIB gmsh HINTS "/home/pedro/projetos/library/gmsh-install/lib" "../lib/gmsh-install/lib" "../gmsh-4.8.4-MacOSX-sdk/lib" NO_DEFAULT_PATH)
find_path(GMSH_INC gmsh.h HINTS "/home/pedro/projetos/library/gmsh-install/include" "../lib/gmsh-install/include" "../gmsh-4.8.4-MacOSX-sdk/include" NO_DEFAULT_PATH)

link_libraries(${GMSH_LIB})
if(NOT GMSH_LIB)
  message(FATAL_ERROR "Could not find libgmsh")
endif()

if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find gmsh.h")
endif()
# include_directories(${GMSH_INC})


# linking std threads
find_package(Threads REQUIRED)

# Copy 'examples' into build directory. 
# @WARNING this is only triggered during configure time.
# I've added it as a request from Phil, but we gotta remember this duplicated file-tree can potentially lead to confusion since now the developer has to mind which example they are editing and which example the CLI is calling.
# Notice that 'build/examples' will only get updated when calling 'cmake' command, and this also means that, if developer was editing 'build/examples' instead of 'source/examples', the former will get overwritten in an acidental/automatic call of 'cmake' command
file(COPY "examples" DESTINATION "./src")

add_subdirectory(src)
# add_subdirectory(testcollapsed)

############# INSTALL #############
message(STATUS "CMAKE_INSTALL_LIBDIR = ${CMAKE_INSTALL_LIBDIR}")
message(STATUS "CMAKE_INSTALL_INCLUDEDIR = ${CMAKE_INSTALL_INCLUDEDIR}")
message(STATUS "CMAKE_INSTALL_BINDIR = ${CMAKE_INSTALL_BINDIR}")


if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(WIN32)
    set(CMAKE_INSTALL_PREFIX "C:/neopz" CACHE PATH "Directory where PZ will be installed" FORCE)
  else()
    set(CMAKE_INSTALL_PREFIX "/opt/neopz" CACHE PATH "Directory where PZ will be installed" FORCE)
    message(FATAL_ERROR "Please set install prefix")
  endif()
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

install(TARGETS dfnlib
        EXPORT DFNTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )
#where cmake config files are going to be set
set(DFN_CMAKE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib/cmake/dfnmesh)

# install(EXPORT DFNTargets
#         FILE DFNTargets.cmake
#         NAMESPACE DFN::
#         DESTINATION ${DFN_CMAKE_INSTALL_DIR}
#         )

# NS: Not sure if we need this      
# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/NeoPZConfig.cmake.in
#   "${CMAKE_CURRENT_BINARY_DIR}/NeoPZConfig.cmake" @ONLY)

# NS: Not sure if we need this      
# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/NeoPZConfigVersion.cmake.in
#   "${CMAKE_CURRENT_BINARY_DIR}/NeoPZConfigVersion.cmake" @ONLY)

# NS: Not sure if we need this      
# install(FILES
#         "${CMAKE_CURRENT_BINARY_DIR}/NeoPZConfig.cmake"
#         "${CMAKE_CURRENT_BINARY_DIR}/NeoPZConfigVersion.cmake"
#         "${CMAKE_CURRENT_SOURCE_DIR}/cmake/add_pz_target.cmake"
#         "${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_pz_opt.cmake"
#         DESTINATION ${DFN_CMAKE_INSTALL_DIR}
#         )

# include(cmake/InstallPZConfigFile.cmake)


############# Testing #############
enable_testing()