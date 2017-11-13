CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT(R)


############## multiscale optimal transport package ##############
#IF( USE_CPLEX )
#  IF( USE_MOSEK )
#    IF( USE_GLPK )
#      COMMAND cp ${R_SOURCE_DIR}/mop/MakevarsALL ./mop/src/Makevars
#      COMMAND cp ${R_SOURCE_DIR}/mop/mop_config_ALL.h ./mop/src/mop_config.h
#    ELSE( USE_GLPK )
#      COMMAND cp ${R_SOURCE_DIR}/mop/MakevarsMOSEK_CPLEX ./mop/src/Makevars
#      COMMAND cp ${R_SOURCE_DIR}/mop/mop_config_MOSEK_CPLEX.h ./mop/src/mop_config.h
#    ENDIF( USE_GLPK )
#  ELSEIF( USE_GLPK )
#    COMMAND cp ${R_SOURCE_DIR}/mop/MakevarsGLPK_CPLEX ./mop/src/Makevars
#    COMMAND cp ${R_SOURCE_DIR}/mop/mop_config_GLPK_CPLEX.h ./mop/src/mop_config.h
#  ELSE( USE_MOSEK )
#    COMMAND cp ${R_SOURCE_DIR}/mop/MakevarsCPLEX ./mop/src/Makevars
#    COMMAND cp ${R_SOURCE_DIR}/mop/mop_config_CPLEX.h ./mop/src/mop_config.h
#  ENDIF( USE_MOSEK )
#ELSEIF( USE_CPLEX )
#  IF( USE_MOSEK )
#    IF( USE_GLPK )
#      COMMAND cp ${R_SOURCE_DIR}/mop/MakevarsMOSEK_GLPK ./mop/src/Makevars
#      COMMAND cp ${R_SOURCE_DIR}/mop/mop_config_MOSEK_GLPK.h ./mop/src/mop_config.h
#    ELSE( USE_GLPK )
#      COMMAND cp ${R_SOURCE_DIR}/mop/MakevarsMOSEK ./mop/src/Makevars
#      COMMAND cp ${R_SOURCE_DIR}/mop/mop_config_MOSEK.h ./mop/src/mop_config.h
#    ENDIF( USE_GLPK )
#  ELSE( USE_MOSEK )
#    COMMAND cp ${R_SOURCE_DIR}/mop/MakevarsGLPK ./mop/src/Makevars
#    COMMAND cp ${R_SOURCE_DIR}/mop/mop_config_GLPK.h ./mop/src/mop_config.h
#  ENDIF( USE_MOSEK)
#ENDIF( USE_CPLEX )



ADD_CUSTOM_TARGET(R_mop
    COMMAND rsync -aC ${R_SOURCE_DIR}/mop . 
    COMMAND rsync -aC ${CMAKE_SOURCE_DIR}/external/gmra/lib/ ./mop/src 
    COMMAND rsync -aC ${R_SOURCE_DIR}/Eigen/ ./mop/src/Eigen 
    COMMAND rsync -aC ${R_SOURCE_DIR}/unsupported/ ./mop/src/unsupported 
    COMMAND cp 
    ${CMAKE_SOURCE_DIR}/lib/SinkhornTransport.h
    ${CMAKE_SOURCE_DIR}/lib/SparseSinkhornTransport.h
    ${CMAKE_SOURCE_DIR}/lib/MultiscaleTransport.h
    ${CMAKE_SOURCE_DIR}/lib/MultiscaleTransportLP.h
    ${CMAKE_SOURCE_DIR}/lib/NeighborhoodStrategy.h
    ${CMAKE_SOURCE_DIR}/lib/ExpandNeighborhoodStrategy.h
    ${CMAKE_SOURCE_DIR}/lib/RefineNeighborhoodStrategy.h
    ${CMAKE_SOURCE_DIR}/lib/PotentialNeighborhoodStrategy.h
    ${CMAKE_SOURCE_DIR}/lib/PropagationStrategy.h
    ${CMAKE_SOURCE_DIR}/lib/NeighborhoodPropagationStrategy.h
    ${CMAKE_SOURCE_DIR}/lib/CapacityPropagationStrategy.h
    ${CMAKE_SOURCE_DIR}/lib/IteratedCapacityPropagationStrategy.h
    ${CMAKE_SOURCE_DIR}/lib/RandomizedNeighborhoodPropagationStrategy.h
    ${CMAKE_SOURCE_DIR}/lib/MaxEntropyPropagationStrategy.h
    ${CMAKE_SOURCE_DIR}/lib/SinkhornPropagationStrategy.h
    ${CMAKE_SOURCE_DIR}/lib/TransportLP.h
    ${CMAKE_SOURCE_DIR}/lib/TransportLPSolver.h
    ${CMAKE_SOURCE_DIR}/lib/GLPKSolver.h
    ${CMAKE_SOURCE_DIR}/lib/MOSEKSolver.h
    ${CMAKE_SOURCE_DIR}/lib/CPLEXSolver.h
    ${CMAKE_SOURCE_DIR}/lib/CPLEXNetworkSolver.h
    ${CMAKE_SOURCE_DIR}/lib/LPSolver.h
    ${CMAKE_SOURCE_DIR}/lib/MatrixMinFlow.h
    ${CMAKE_SOURCE_DIR}/lib/Grid3dMinFlow.h
    ${CMAKE_SOURCE_DIR}/lib/mosek.h
    ${CMAKE_SOURCE_DIR}/external/gmra/external/metrics/lib/EigenMetric.h
    ${CMAKE_SOURCE_DIR}/external/gmra/external/metrics/lib/EigenEuclideanMetric.h 
    ${CMAKE_SOURCE_DIR}/external/gmra/external/metrics/lib/EigenSquaredEuclideanMetric.h 
    ${CMAKE_SOURCE_DIR}/external/gmra/external/metrics/lib/EigenL1Metric.h 
    ${CMAKE_SOURCE_DIR}/external/gmra/external/metrics/lib/EigenWasserstein.h
    ${CMAKE_SOURCE_DIR}/external/gmra/external/utils/lib/Random.h
    ${CMAKE_SOURCE_DIR}/external/gmra/external/randomsvd/lib/EigenRandomSVD.h
    ${CMAKE_SOURCE_DIR}/external/gmra/external/randomsvd/lib/EigenRandomRange.h
    ./mop/src
)

ADD_CUSTOM_TARGET(R_mop_build  
    COMMAND R CMD build --resave-data mop)
ADD_DEPENDENCIES(R_mop_build R_mop)

ADD_CUSTOM_TARGET(R_mop_install
    COMMAND R CMD INSTALL  --no-multiarch   mop)
ADD_DEPENDENCIES(R_mop_install R_mop)

ADD_CUSTOM_TARGET(R_mop_check 
    COMMAND R CMD check --no-multiarch mop)
ADD_DEPENDENCIES(R_mop_check R_mop)



########### Entropy regularized optimal transport #############

ADD_CUSTOM_TARGET(R_sinkhorn
    COMMAND rsync -aC ${R_SOURCE_DIR}/sinkhorn . 
    COMMAND rsync -aC ${CMAKE_SOURCE_DIR}/external/gmra/lib/ ./sinkhorn/src 
    COMMAND rsync -aC ${R_SOURCE_DIR}/Eigen/ ./sinkhorn/src/Eigen 
    COMMAND cp 
    ${CMAKE_SOURCE_DIR}/lib/SinkhornTransport.h
    ${CMAKE_SOURCE_DIR}/lib/SparseSinkhornTransport.h
    ${CMAKE_SOURCE_DIR}/lib/MultiscaleTransport.h
    ${CMAKE_SOURCE_DIR}/lib/MultiscaleSinkhornTransport.h
    ${CMAKE_SOURCE_DIR}/external/gmra/external/metrics/lib/EigenMetric.h
    ${CMAKE_SOURCE_DIR}/external/gmra/external/metrics/lib/EigenEuclideanMetric.h 
    ${CMAKE_SOURCE_DIR}/external/gmra/external/metrics/lib/EigenSquaredEuclideanMetric.h 
    ${CMAKE_SOURCE_DIR}/external/gmra/external/metrics/lib/EigenL1Metric.h 
    ${CMAKE_SOURCE_DIR}/external/gmra/external/utils/lib/Random.h
    ${CMAKE_SOURCE_DIR}/external/gmra/external/metrics/lib/EigenWasserstein.h
    ${CMAKE_SOURCE_DIR}/external/gmra/external/randomsvd/lib/EigenRandomSVD.h
    ${CMAKE_SOURCE_DIR}/external/gmra/external/randomsvd/lib/EigenRandomRange.h
    ./sinkhorn/src
)

ADD_CUSTOM_TARGET(R_sinkhorn_build  
    COMMAND R CMD build --resave-data sinkhorn)
ADD_DEPENDENCIES(R_sinkhorn_build R_sinkhorn)

ADD_CUSTOM_TARGET(R_sinkhorn_install
    COMMAND R CMD INSTALL  --no-multiarch   sinkhorn)
ADD_DEPENDENCIES(R_sinkhorn_install R_sinkhorn)

ADD_CUSTOM_TARGET(R_sinkhorn_check 
    COMMAND R CMD check --no-multiarch sinkhorn)
ADD_DEPENDENCIES(R_sinkhorn_check R_sinkhorn)









