cmake_minimum_required(VERSION 2.6)

PROJECT(GMRA)


if(BUILD_MCFSimplex)
LINK_DIRECTORIES( ${MCFSIMPLEX_LINK_DIR} )
endif(BUILD_MCFSimplex)

ADD_EXECUTABLE(CreateAATree CreateAATree.cxx)
TARGET_LINK_LIBRARIES (CreateAATree lapack blas)

ADD_EXECUTABLE(CreateIPCATree CreateIPCATree.cxx)
TARGET_LINK_LIBRARIES ( CreateIPCATree )

ADD_EXECUTABLE(CreateIKMTree CreateIKMTree.cxx)
TARGET_LINK_LIBRARIES ( CreateIKMTree )

ADD_EXECUTABLE(IPCATreeStatistics IPCATreeStatistics.cxx)
TARGET_LINK_LIBRARIES (IPCATreeStatistics lapack blas)

ADD_EXECUTABLE(IPCATreeWasserstein IPCATreeWasserstein.cxx)
TARGET_LINK_LIBRARIES (IPCATreeWasserstein lapack blas)


ADD_EXECUTABLE(IPCATreeLeafDistances IPCATreeLeafDistances.cxx)
TARGET_LINK_LIBRARIES (IPCATreeLeafDistances lapack blas)

ADD_EXECUTABLE(IPCATreeCrossWasserstein IPCATreeCrossWasserstein.cxx)
TARGET_LINK_LIBRARIES (IPCATreeCrossWasserstein lapack blas)

ADD_EXECUTABLE(IPCATreePairwiseWasserstein IPCATreePairwiseWasserstein.cxx)
TARGET_LINK_LIBRARIES (IPCATreePairwiseWasserstein lapack blas)

ADD_EXECUTABLE(IPCAGWT IPCAGWT.cxx)
TARGET_LINK_LIBRARIES (IPCAGWT lapack blas)
  
ADD_EXECUTABLE(IPCADensityEstimator IPCADensityEstimator.cxx)
TARGET_LINK_LIBRARIES (IPCADensityEstimator lapack blas)

ADD_EXECUTABLE(IPCAGWTSquared IPCAGWTSquared.cxx)
TARGET_LINK_LIBRARIES (IPCAGWTSquared lapack blas)
 
ADD_EXECUTABLE(IPCADists IPCADists.cxx)
TARGET_LINK_LIBRARIES (IPCADists lapack blas)

ADD_EXECUTABLE(IPCAGWTCompress IPCAGWTCompress.cxx)
TARGET_LINK_LIBRARIES (IPCAGWTCompress lapack blas)

ADD_EXECUTABLE(NodeSimilarityIPCA NodeSimilarityIPCA.cxx)
TARGET_LINK_LIBRARIES (NodeSimilarityIPCA lapack blas ${ITK_LIBRARIES})


ADD_EXECUTABLE(CrossWaveletSimilarityIPCA  CrossWaveletSimilarityIPCA.cxx)
TARGET_LINK_LIBRARIES (CrossWaveletSimilarityIPCA lapack blas)

ADD_EXECUTABLE(LPDistanceIPCA LPDistanceIPCA.cxx)
TARGET_LINK_LIBRARIES (LPDistanceIPCA lapack blas glpk)


option(BUILD_GMRA_TRANSPORT "Build GMRA optimal transport" ON)
if(BUILD_GMRA_TRANSPORT)

 FIND_PACKAGE(CPLEX)
INCLUDE_DIRECTORIES( ${CPLEX_INCLUDE_DIRS} )

ADD_EXECUTABLE(IPCATransport IPCATransport.cxx)
TARGET_LINK_LIBRARIES (IPCATransport lapack blas ${CPLEX_LIBRARIES} pthread )

endif(BUILD_GMRA_TRANSPORT)



ADD_CUSTOM_TARGET(GMRA_release 
    COMMAND rsync -aC ${GMRA_SOURCE_DIR}/ ./gmra_src/
    COMMAND mkdir -p ./gmra_src/lib 
    COMMAND rsync -aC ${DIAT_LIB_DIR}/gmra/ ./gmra_src/lib/
    COMMAND rsync -aC ${DIAT_SOURCE_DIR}/tools/config/ ./gmra_src/config/
    COMMAND rsync -aC ${DIAT_SOURCE_DIR}/lib/tclap/ ./gmra_src/tclap/
    COMMAND cp 
${DIAT_SOURCE_DIR}/cmake_modules/*
${DIAT_LIB_DIR}/linalg/DenseMatrix.h
${DIAT_LIB_DIR}/linalg/Linalg.h
${DIAT_LIB_DIR}/linalg/LapackDefs.h
${DIAT_LIB_DIR}/linalg/Matrix.h 
${DIAT_LIB_DIR}/linalg/Vector.h 
${DIAT_LIB_DIR}/linalg/DenseVector.h
${DIAT_LIB_DIR}/linalg/SVD.h
${DIAT_LIB_DIR}/linalg/RandomSVD.h
${DIAT_LIB_DIR}/linalg/RandomRange.h
${DIAT_LIB_DIR}/statistics/Random.h
${DIAT_LIB_DIR}/statistics/Wasserstein.h
${DIAT_LIB_DIR}/metrics/EuclideanMetric.h
${DIAT_LIB_DIR}/metrics/SquaredEuclideanMetric.h
${DIAT_LIB_DIR}/metrics/Metric.h
${DIAT_LIB_DIR}/io/LinalgIO.h
./gmra_src/lib/   
    COMMAND mv ./gmra_src/CMakeListsRelease.txt ./gmra_src/CMakeLists.txt
    COMMAND tar -c gmra_src > gmra_src.tar
)

