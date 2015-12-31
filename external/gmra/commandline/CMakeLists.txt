


ADD_EXECUTABLE(CreateIPCATree CreateIPCATree.cxx)
TARGET_LINK_LIBRARIES ( CreateIPCATree )

ADD_EXECUTABLE(CreateIKMTree CreateIKMTree.cxx)
TARGET_LINK_LIBRARIES ( CreateIKMTree )

ADD_EXECUTABLE(IPCADensityEstimator IPCADensityEstimator.cxx)
TARGET_LINK_LIBRARIES (IPCADensityEstimator lapack blas)


##old code
#ADD_EXECUTABLE(CreateAATree CreateAATree.cxx)
#TARGET_LINK_LIBRARIES (CreateAATree lapack blas)

#ADD_EXECUTABLE(IPCATreeStatistics IPCATreeStatistics.cxx)
#TARGET_LINK_LIBRARIES (IPCATreeStatistics lapack blas)

#ADD_EXECUTABLE(IPCATreeWasserstein IPCATreeWasserstein.cxx)
#TARGET_LINK_LIBRARIES (IPCATreeWasserstein lapack blas)


#ADD_EXECUTABLE(IPCATreeLeafDistances IPCATreeLeafDistances.cxx)
#TARGET_LINK_LIBRARIES (IPCATreeLeafDistances lapack blas)

#ADD_EXECUTABLE(IPCATreeCrossWasserstein IPCATreeCrossWasserstein.cxx)
#TARGET_LINK_LIBRARIES (IPCATreeCrossWasserstein lapack blas)

#ADD_EXECUTABLE(IPCATreePairwiseWasserstein IPCATreePairwiseWasserstein.cxx)
#TARGET_LINK_LIBRARIES (IPCATreePairwiseWasserstein lapack blas)

#ADD_EXECUTABLE(IPCAGWT IPCAGWT.cxx)
#TARGET_LINK_LIBRARIES (IPCAGWT lapack blas)
  

#ADD_EXECUTABLE(IPCAGWTSquared IPCAGWTSquared.cxx)
#TARGET_LINK_LIBRARIES (IPCAGWTSquared lapack blas)
 
#ADD_EXECUTABLE(IPCADists IPCADists.cxx)
#TARGET_LINK_LIBRARIES (IPCADists lapack blas)

#ADD_EXECUTABLE(IPCAGWTCompress IPCAGWTCompress.cxx)
#TARGET_LINK_LIBRARIES (IPCAGWTCompress lapack blas)

#ADD_EXECUTABLE(NodeSimilarityIPCA NodeSimilarityIPCA.cxx)
#TARGET_LINK_LIBRARIES (NodeSimilarityIPCA lapack blas ${ITK_LIBRARIES})


#ADD_EXECUTABLE(CrossWaveletSimilarityIPCA  CrossWaveletSimilarityIPCA.cxx)
#TARGET_LINK_LIBRARIES (CrossWaveletSimilarityIPCA lapack blas)

#ADD_EXECUTABLE(LPDistanceIPCA LPDistanceIPCA.cxx)
#TARGET_LINK_LIBRARIES (LPDistanceIPCA lapack blas glpk)



