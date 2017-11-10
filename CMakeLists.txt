CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT(OptimalTransport)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/")
find_package(Eigen3 REQUIRED)

option(USE_CPLEX "Use CPLEX libraries")

IF(USE_CPLEX)
  FIND_PACKAGE(CPLEX)
ENDIF(USE_CPLEX)


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

