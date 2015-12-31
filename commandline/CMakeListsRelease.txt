cmake_minimum_required(VERSION 2.6)

PROJECT(GMRA)

FIND_PACKAGE(LAPACK)
FIND_PACKAGE(BLAS)

INCLUDE_DIRECTORIES(. config lib tclap)

ADD_EXECUTABLE(CreateIPCATree CreateIPCATree.cxx)
TARGET_LINK_LIBRARIES (CreateIPCATree lapack blas)

ADD_EXECUTABLE(IPCATreeStatistics IPCATreeStatistics.cxx)
TARGET_LINK_LIBRARIES (IPCATreeStatistics lapack blas)

ADD_EXECUTABLE(IPCAGWT IPCAGWT.cxx)
TARGET_LINK_LIBRARIES (IPCAGWT lapack blas)




