CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT(OptimalTransport)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/")
find_package(Eigen3 REQUIRED)
find_package(CPLEX REQUIRED)

#FIND_PATH(LEMON_SOURCE_ROOT_DIR CMakeLists.txt
#  PATHS ${CMAKE_SOURCE_DIR}/lemon ${CMAKE_SOURCE_DIR}/deps/lemon 
#  NO_DEFAULT_PATH
#  DOC "Location of LEMON source as a CMAKE subproject")
  
#IF(EXISTS ${LEMON_SOURCE_ROOT_DIR})
#    ADD_SUBDIRECTORY(${LEMON_SOURCE_ROOT_DIR} deps/lemon)
#    SET(LEMON_INCLUDE_DIRS
#      ${LEMON_SOURCE_ROOT_DIR}
#      ${CMAKE_BINARY_DIR}/deps/lemon
#    )
#    SET(LEMON_LIBRARIES lemon)
#    UNSET(LEMON_ROOT_DIR CACHE)
#    UNSET(LEMON_DIR CACHE)
#    UNSET(LEMON_INCLUDE_DIR CACHE)
#    UNSET(LEMON_LIBRARY CACHE)
#ELSE()
#    FIND_PACKAGE(LEMON QUIET NO_MODULE)
#    FIND_PACKAGE(LEMON REQUIRED)
#ENDIF()


INCLUDE_DIRECTORIES( lib
                     external/gmra/lib
                     external/gmra/external/utils/lib
                     external/gmra/external/kmeans/lib
                     external/gmra/external/randomsvd/lib
                     external/gmra/external/metrics/lib 
                     ${EIGEN3_INCLUDE_DIR} 
                     ${CPLEX_INCLUDE_DIRS} 
                     )
ADD_SUBDIRECTORY(commandline)
ADD_SUBDIRECTORY(Rpackage)

