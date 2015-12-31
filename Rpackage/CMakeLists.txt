CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT(R)




############## gmra package ##############


ADD_CUSTOM_TARGET(R_gmra 
    COMMAND rsync -aC ${R_SOURCE_DIR}/gmra . 
    COMMAND rsync -aC ${CMAKE_SOURCE_DIR}/lib/ ./gmra/src 
    COMMAND cp 
    ${CMAKE_SOURCE_DIR}/external/metrics/lib/EigenMetric.h
    ${CMAKE_SOURCE_DIR}/external/metrics/lib/EigenEuclideanMetric.h 
    ${CMAKE_SOURCE_DIR}/external/metrics/lib/EigenL1Metric.h 
    ${CMAKE_SOURCE_DIR}/external/metrics/lib/EigenWasserstein.h
    ${CMAKE_SOURCE_DIR}/external/randomsvd/lib/EigenSquaredEuclideanMetric.h 
    ${CMAKE_SOURCE_DIR}/external/randomsvd/lib/EigenRandomSVD.h
    ${CMAKE_SOURCE_DIR}/external/randomsvd/lib/EigenRandomRange.h
    ${CMAKE_SOURCE_DIR}/external/utils/lib/Random.h
    ${CMAKE_SOURCE_DIR}/external/kmeans/lib/Kmeans.h
    ./gmra/src
)

ADD_CUSTOM_TARGET(R_gmra_build  
    COMMAND R CMD build --resave-data gmra)
ADD_DEPENDENCIES(R_gmra_build R_gmra)

ADD_CUSTOM_TARGET(R_gmra_install
    COMMAND R CMD INSTALL  --no-multiarch  gmra)
ADD_DEPENDENCIES(R_gmra_install R_gmra)

ADD_CUSTOM_TARGET(R_gmra_check 
    COMMAND R CMD check --no-multiarch gmra)
ADD_DEPENDENCIES(R_gmra_check R_gmra)



