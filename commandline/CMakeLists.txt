ADD_EXECUTABLE(IPCATransport IPCATransport.cxx)
TARGET_LINK_LIBRARIES (IPCATransport lapack blas pthread )

ADD_EXECUTABLE(TransportLP TransportLP.cxx)

ADD_EXECUTABLE(MatrixMinFlow MatrixMinFlow.cxx)

ADD_EXECUTABLE(Grid3dMinFlow Grid3dMinFlow.cxx)

IF(USE_CPLEX)
TARGET_LINK_LIBRARIES (IPCATransport ${CPLEX_LIBRARIES} )
TARGET_LINK_LIBRARIES (TransportLP ${CPLEX_LIBRARIES} )
TARGET_LINK_LIBRARIES (Grid3dMinFlow ${CPLEX_LIBRARIES} )
ENDIF(USE_CPLEX)


ADD_EXECUTABLE(SinkhornTransport SinkhornTransport.cxx)
TARGET_LINK_LIBRARIES (SinkhornTransport)
  
ADD_EXECUTABLE(SparseSinkhornTransport SparseSinkhornTransport.cxx)
TARGET_LINK_LIBRARIES (SparseSinkhornTransport)


