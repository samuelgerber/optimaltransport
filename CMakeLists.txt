CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT(RandomSVD)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/")
find_package(Eigen3 REQUIRED)

FIND_PACKAGE(CPLEX)


INCLUDE_DIRECTORIES( lib
                     external/utils/lib
                     external/kmeans/lib
                     external/flinalg/lib
                     external/randomsvd/lib
                     external/metrics/lib 
                     ${EIGEN3_INCLUDE_DIR} 
                     ${CPLEX_INCLUDE_DIRS} 
                     )
ADD_SUBDIRECTORY(commandline)
ADD_SUBDIRECTORY(Rpackage)

