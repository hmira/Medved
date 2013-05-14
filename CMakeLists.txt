#
# nešlo vygenerovať
#
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

ADD_DEFINITIONS("-std=c++0x -fopenmp")

PROJECT ( Hmira )

SET( HMIRA_INCLUDE_DIRECTORY "include/" )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests )

find_package(OpenMP)
if(OPENMP_FOUND)
endif(OPENMP_FOUND)

FIND_LIBRARY(OPENMESH_CORE_LIBRARY_DEBUG NAMES OpenMeshCored PATHS "/usr/local/lib/OpenMesh")
FIND_LIBRARY(OPENMESH_TOOLS_LIBRARY_DEBUG NAMES OpenMeshToolsd PATHS "/usr/local/lib/OpenMesh")

# add_subdirectory( ${HMIRA_INCLUDE_DIRECTORY} )

find_package(TBB)
if(NOT TBB_FOUND)
MESSAGE (STATUS "TBB not found")
else(NOT TBB_FOUND)

include_directories(${TBB_INCLUDE_DIRS})
link_directories(${TBB_LIBRARY_DIRS})
endif(NOT TBB_FOUND)

INCLUDE_DIRECTORIES(
	${HMIRA_INCLUDE_DIRECTORY}
	${CMAKE_CURRENT_BINARY_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${OPENMESH_INCLUDE_DIR} )

find_package(Boost COMPONENTS program_options REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})                                                                                                                                                                                                                            
endif()

ADD_EXECUTABLE( main "main.cpp" ${HMIRA_HEADERS} ${EXTENSION_HEADERS} )
TARGET_LINK_LIBRARIES( main gomp)
TARGET_LINK_LIBRARIES( main ${TBB_LIBRARIES})
TARGET_LINK_LIBRARIES( main ${OPENMESH_CORE_LIBRARY_DEBUG} )
TARGET_LINK_LIBRARIES( main ${OPENMESH_TOOLS_LIBRARY_DEBUG} )

# add_executable( rasterize_and_triangulate "tests/rasterize_and_triangulate.cpp" ${HMIRA_HEADERS} ${EXTENSION_HEADERS} )
# target_link_libraries( rasterize_and_triangulate gomp ${TBB_LIBRARIES} ${OPENMESH_CORE_LIBRARY_DEBUG} ${OPENMESH_TOOLS_LIBRARY_DEBUG})
# target_link_libraries( rasterize_and_triangulate boost_program_options ${Boost_LIBRARIES})
# 
# add_executable( rasterize "tests/rasterize.cpp" ${HMIRA_HEADERS} ${EXTENSION_HEADERS} )
# target_link_libraries( rasterize gomp ${TBB_LIBRARIES} ${OPENMESH_CORE_LIBRARY_DEBUG} ${OPENMESH_TOOLS_LIBRARY_DEBUG})
# target_link_libraries( rasterize boost_program_options ${Boost_LIBRARIES})
# 
# add_executable( triangulate "tests/triangulate.cpp" ${HMIRA_HEADERS} ${EXTENSION_HEADERS} )
# target_link_libraries( triangulate gomp ${TBB_LIBRARIES} ${OPENMESH_CORE_LIBRARY_DEBUG} ${OPENMESH_TOOLS_LIBRARY_DEBUG})
# target_link_libraries( triangulate boost_program_options ${Boost_LIBRARIES})