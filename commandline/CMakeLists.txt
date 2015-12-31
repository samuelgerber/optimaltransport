ADD_EXECUTABLE(IPCATransport IPCATransport.cxx)
TARGET_LINK_LIBRARIES (IPCATransport lapack blas ${CPLEX_LIBRARIES} pthread )

ADD_EXECUTABLE(TransportLP TransportLP.cxx)
TARGET_LINK_LIBRARIES (TransportLP gfortran lapack blas ${CPLEX_LIBRARIES} pthread)

ADD_EXECUTABLE(SinkhornTransport SinkhornTransport.cxx)
TARGET_LINK_LIBRARIES (SinkhornTransport)
  
ADD_EXECUTABLE(SparseSinkhornTransport SparseSinkhornTransport.cxx)
TARGET_LINK_LIBRARIES (SparseSinkhornTransport)


