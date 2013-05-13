#
# nešlo vygenerovať
#
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

ADD_DEFINITIONS("-std=c++0x -fopenmp")

PROJECT ( Hmira )

#IF(NOT MEDV4D_ROOT_DIRECTORY)
#        SET(MEDV4D_ROOT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../../MedV4D" )
#ENDIF(NOT MEDV4D_ROOT_DIRECTORY)

#IF(NOT MEDV4D_CMAKE_MODULES_DIR )
#        SET(MEDV4D_CMAKE_MODULES_DIR ${MEDV4D_ROOT_DIRECTORY}/../cmake_project/Modules)
#        SET(CMAKE_MODULE_PATH "${MEDV4D_CMAKE_MODULES_DIR}" "${CMAKE_ROOT}/Modules" )
#ENDIF(NOT MEDV4D_CMAKE_MODULES_DIR )

#FIND_PACKAGE( MedV4D REQUIRED )

SET( HMIRA_INCLUDE_DIRECTORY "include" )


# Add FindTBB directory to CMake's module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

find_package(OpenMP)
if(OPENMP_FOUND)
message("damn im good")
endif(OPENMP_FOUND)


# Execute the FindTBB macro
find_package(TBB)
# If you want FindTBB to fail if TBB is not found, use this form:
# find_package(TBB REQUIRED)

if(NOT TBB_FOUND)
        MESSAGE (STATUS "TBB not found. Example project will not be built.")
else(NOT TBB_FOUND)

# Here goes the actual project
include_directories(${TBB_INCLUDE_DIRS})
link_directories(${TBB_LIBRARY_DIRS})
endif(NOT TBB_FOUND)



FIND_LIBRARY( OPENMESH_CORE_DEBUG_LIB OpenMeshCored  PATHS "/usr/local/lib/OpenMesh")
FIND_LIBRARY( OPENMESH_CORE_RELEASE_LIB OpenMeshCore )
SET( OPENMESH_CORE_LIB "debug" ${OPENMESH_CORE_DEBUG_LIB} "debug" ${OPENMESH_CORE_DEBUG_LIB} )

#SET( HMIRA_GRID 
#	"include/grids/RegularGrid.cc" 
#	"include/grids/ScalarGridT.cc" 
#	)

FIND_PATH( OPENMESH_INCLUDE_DIR "OpenMesh/Core/Geometry/VectorT.hh" )
#FIND_PATH( OPENMESH_INCLUDE_DIR "/usr/local/include/OpenMesh" )

INCLUDE_DIRECTORIES(${MEDV4D_INCLUDE_DIR} 
#${MEDV4D_INCLUDE_DIRS_3RD_PARTY} 
${HMIRA_INCLUDE_DIRECTORY} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${OPENMESH_INCLUDE_DIR} )

#INCLUDE( "${MEDV4D_CMAKE_MODULES_DIR}/StandardCompileOptions.cmake" NO_POLICY_SCOPE )

find_package(Boost COMPONENTS program_options REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})                                                                                                                                                                                                                            
endif()

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests )


ADD_EXECUTABLE( main "main.cpp" ${HMIRA_HEADERS} ${EXTENSION_HEADERS} )
add_executable( rasterize_and_triangulate "tests/rasterize_and_triangulate.cpp" ${HMIRA_HEADERS} ${EXTENSION_HEADERS} )

MESSAGE( "${MEDV4D_EXECUTABLE_SOURCES}" )
MESSAGE( "${MEDV4D_EXECUTABLE_HEADERS}" )
MESSAGE( "${MEDV4D_LIBRARIES}" )


#TARGET_LINK_LIBRARIES( main 
	#${MEDV4D_LIBRARIES} 
	#${MEDV4D_LIBRARIES_3RD_PARTY} 
#	)

#TARGET_LINK_LIBRARIES( main
#	OpenMeshTools
#	OpenMeshCore)

FIND_LIBRARY(OPENMESH_CORE_LIBRARY_DEBUG NAMES OpenMeshCored PATHS "/usr/local/lib/OpenMesh")
FIND_LIBRARY(OPENMESH_TOOLS_LIBRARY_DEBUG NAMES OpenMeshToolsd PATHS "/usr/local/lib/OpenMesh")

SUBDIRS( ${HMIRA_INCLUDE_DIRECTORY} )

TARGET_LINK_LIBRARIES( main gomp)
TARGET_LINK_LIBRARIES( main ${TBB_LIBRARIES})
TARGET_LINK_LIBRARIES( main ${OPENMESH_CORE_LIBRARY_DEBUG} )
TARGET_LINK_LIBRARIES( main ${OPENMESH_TOOLS_LIBRARY_DEBUG} )

target_link_libraries( rasterize_and_triangulate gomp ${TBB_LIBRARIES} ${OPENMESH_CORE_LIBRARY_DEBUG} ${OPENMESH_TOOLS_LIBRARY_DEBUG})
target_link_libraries( rasterize_and_triangulate boost_program_options ${Boost_LIBRARIES})