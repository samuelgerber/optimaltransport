CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT(Rcode)


######## msr package #########

ADD_CUSTOM_TARGET(R_msr 
    COMMAND rsync -aC ${Rcode_SOURCE_DIR}/msr . 
    COMMAND rsync -aC ${DIAT_LIB_DIR}/ann/ ./msr/src 
    COMMAND cp 
    ${DIAT_LIB_DIR}/linalg/DenseMatrix.h
    ${DIAT_LIB_DIR}/metrics/EuclideanMetric.h 
    ${DIAT_LIB_DIR}/utils/Geometry.h
    ${DIAT_LIB_DIR}/statistics/KernelDensity.h 
    ${DIAT_LIB_DIR}/linalg/Linalg.h
    ${DIAT_LIB_DIR}/linalg/LapackDefs.h
    ${DIAT_LIB_DIR}/linalg/Matrix.h
    ${DIAT_LIB_DIR}/metrics/Metric.h
    ${DIAT_LIB_DIR}/metrics/SquaredEuclideanMetric.h
    ${DIAT_LIB_DIR}/linalg/Vector.h 
    ${DIAT_LIB_DIR}/linalg/DenseVector.h
    ${DIAT_LIB_DIR}/statistics/GaussianKernel.h
    ${DIAT_LIB_DIR}/datastructures/heap/Heap.h
    ${DIAT_LIB_DIR}/statistics/Kernel.h
    ${DIAT_LIB_DIR}/datastructures/heap/MaxHeap.h
    ${DIAT_LIB_DIR}/datastructures/heap/MinHeap.h
    ${DIAT_LIB_DIR}/topology/NNMSComplex.h
    ${DIAT_LIB_DIR}/topology/NNMSComplex2old.h
    ${DIAT_LIB_DIR}/io/LinalgIO.h
    ./msr/src
    COMMAND rm ./msr/src/CMakeLists.txt
    COMMAND rm ./msr/src/perf.cpp
    COMMAND rm ./msr/src/kd_dump.cpp
    COMMAND mv ./msr/src/ReadMe_ANN.txt .
)

ADD_CUSTOM_TARGET(R_msr_build  
    COMMAND R CMD build --resave-data msr)
ADD_DEPENDENCIES(R_msr_build R_msr)

ADD_CUSTOM_TARGET(R_msr_install
    COMMAND R CMD INSTALL msr)
ADD_DEPENDENCIES(R_msr_install R_msr)

ADD_CUSTOM_TARGET(R_msr_check 
    COMMAND R CMD check --as-cran msr)
ADD_DEPENDENCIES(R_msr_check R_msr)







############## cem package ##############

ADD_CUSTOM_TARGET(R_cems 
    COMMAND rsync -aC ${Rcode_SOURCE_DIR}/cems . 
    COMMAND rsync -aC ${DIAT_LIB_DIR}/ann/ ./cems/src 
    COMMAND cp 
    ${DIAT_LIB_DIR}/statistics/Random.h
    ${DIAT_LIB_DIR}/linalg/Matrix.h
    ${DIAT_LIB_DIR}/linalg/DenseMatrix.h
    ${DIAT_LIB_DIR}/linalg/Vector.h 
    ${DIAT_LIB_DIR}/linalg/DenseVector.h
    ${DIAT_LIB_DIR}/metrics/Metric.h
    ${DIAT_LIB_DIR}/metrics/EuclideanMetric.h 
    ${DIAT_LIB_DIR}/metrics/SquaredEuclideanMetric.h 
    ${DIAT_LIB_DIR}/utils/Geometry.h
    ${DIAT_LIB_DIR}/linalg/Linalg.h
    ${DIAT_LIB_DIR}/linalg/LapackDefs.h
    ${DIAT_LIB_DIR}/linalg/SymmetricEigensystem.h
    ${DIAT_LIB_DIR}/linalg/SymmetricDenseMatrix.h
    ${DIAT_LIB_DIR}/statistics/Kernel.h
    ${DIAT_LIB_DIR}/statistics/GaussianKernel.h
    ${DIAT_LIB_DIR}/statistics/MahalanobisKernel.h
    ${DIAT_LIB_DIR}/statistics/KernelDensity.h
    ${DIAT_LIB_DIR}/statistics/CEM.h
    ${DIAT_LIB_DIR}/statistics/FastCEM.h
    ${DIAT_LIB_DIR}/statistics/BlockCEM.h
    ${DIAT_LIB_DIR}/statistics/CEMRegression.h
    ${DIAT_LIB_DIR}/datastructures/heap/Heap.h
    ${DIAT_LIB_DIR}/datastructures/heap/MaxHeap.h
    ${DIAT_LIB_DIR}/datastructures/heap/MinHeap.h
    ./cems/src
    COMMAND rm ./cems/src/CMakeLists.txt
    COMMAND rm ./cems/src/perf.cpp
    COMMAND rm ./cems/src/kd_dump.cpp
    COMMAND rm ./cems/R/cemr* 
    COMMAND rm ./cems/R/cem.old* 
    COMMAND rm ./cems/src/cems-* 
    COMMAND mv ./cems/src/ReadMe_ANN.txt .
)

ADD_CUSTOM_TARGET(R_cems_build  
    COMMAND R CMD build --resave-data cems)
ADD_DEPENDENCIES(R_cems_build R_cems)

ADD_CUSTOM_TARGET(R_cems_install
    COMMAND R CMD INSTALL  --no-multiarch  cems)
ADD_DEPENDENCIES(R_cems_install R_cems)

ADD_CUSTOM_TARGET(R_cems_check 
    COMMAND R CMD check --no-multiarch cems)
ADD_DEPENDENCIES(R_cems_check R_cems)









############## MDM project #####################

ADD_CUSTOM_TARGET(R_mdm 
    COMMAND rsync -aC ${Rcode_SOURCE_DIR}/mdm . 
    COMMAND rsync -aC ${DIAT_LIB_DIR}/ann/ ./mdm/src 
    COMMAND cp 
    ${DIAT_LIB_DIR}/statistics/Random.h
    ${DIAT_LIB_DIR}/linalg/Matrix.h
    ${DIAT_LIB_DIR}/linalg/DenseMatrix.h
    ${DIAT_LIB_DIR}/linalg/Vector.h 
    ${DIAT_LIB_DIR}/linalg/DenseVector.h
    ${DIAT_LIB_DIR}/metrics/Metric.h
    ${DIAT_LIB_DIR}/metrics/EuclideanMetric.h 
    ${DIAT_LIB_DIR}/metrics/SquaredEuclideanMetric.h 
    ${DIAT_LIB_DIR}/utils/Geometry.h
    ${DIAT_LIB_DIR}/linalg/Linalg.h
    ${DIAT_LIB_DIR}/linalg/LapackDefs.h
    ./mdm/src
    COMMAND rm ./mdm/src/CMakeLists.txt
)

ADD_CUSTOM_TARGET(R_mdm_build  
    COMMAND R CMD build --resave-data mdm)
ADD_DEPENDENCIES(R_mdm_build R_mdm)

ADD_CUSTOM_TARGET(R_mdm_install
    COMMAND R CMD INSTALL --no-multiarch  mdm)
ADD_DEPENDENCIES(R_mdm_install R_mdm)

ADD_CUSTOM_TARGET(R_mdm_check 
    COMMAND R CMD check --no-multiarch mdm)
ADD_DEPENDENCIES(R_mdm_check R_mdm)









################  AUCLOSS project  ######################

ADD_CUSTOM_TARGET(R_aucloss 
    COMMAND rsync -aC ${Rcode_SOURCE_DIR}/aucloss . 
)

ADD_CUSTOM_TARGET(R_aucloss_tar
    COMMAND tar -c aucloss > aucloss.tar)
ADD_DEPENDENCIES(R_aucloss_tar R_aucloss)

ADD_CUSTOM_TARGET(R_aucloss_install  
    COMMAND R CMD INSTALL aucloss)
ADD_DEPENDENCIES(R_aucloss_install R_aucloss)









################# mio package ###############

ADD_CUSTOM_TARGET(R_mio 
    COMMAND rsync -aC ${Rcode_SOURCE_DIR}/mio . 
    COMMAND cp 
    ${DIAT_LIB_DIR}/io/LinalgIO.h
    ${DIAT_LIB_DIR}/linalg/DenseMatrix.h
    ${DIAT_LIB_DIR}/linalg/Matrix.h
    ${DIAT_LIB_DIR}/linalg/Vector.h 
    ${DIAT_LIB_DIR}/linalg/DenseVector.h
    ./mio/src
)

ADD_CUSTOM_TARGET(R_mio_build  
    COMMAND R CMD build --resave-data mio)
ADD_DEPENDENCIES(R_mio_build R_mio)

ADD_CUSTOM_TARGET(R_mio_install
    COMMAND R CMD INSTALL mio)
ADD_DEPENDENCIES(R_mio_install R_mio)

ADD_CUSTOM_TARGET(R_mio_check 
    COMMAND R CMD check mio)
ADD_DEPENDENCIES(R_mio_check R_mio)




################# itkImage package ###############

ADD_CUSTOM_TARGET(R_itkImage 
    COMMAND rsync -aC ${Rcode_SOURCE_DIR}/itkImage . 
    COMMAND cp 
    ${DIAT_LIB_DIR}/linalg/DenseMatrix.h
    ${DIAT_LIB_DIR}/linalg/Matrix.h
    ./itkImage/src
)

ADD_CUSTOM_TARGET(R_itkImage_build  
    COMMAND R CMD build --resave-data itkImage)
ADD_DEPENDENCIES(R_itkImage_build R_itkImage)

ADD_CUSTOM_TARGET(R_itkImage_install
    COMMAND R CMD INSTALL itkImage)
ADD_DEPENDENCIES(R_itkImage_install R_itkImage)

ADD_CUSTOM_TARGET(R_itkImage_check 
    COMMAND R CMD check itkImage)
ADD_DEPENDENCIES(R_itkImage_check R_itkImage)






################ gyroscope #################

ADD_CUSTOM_TARGET(R_gyroscope
    COMMAND rsync -aC ${Rcode_SOURCE_DIR}/gyroscope .
    COMMAND rsync -aC ${DIAT_LIB_DIR}/ann/ ./gyroscope/src
    COMMAND rsync -aC ${DIAT_LIB_DIR}/gmra/ ./gyroscope/src
    COMMAND cp 
    ${DIAT_LIB_DIR}/statistics/PCA.h
    ${DIAT_LIB_DIR}/statistics/Random.h
    ${DIAT_LIB_DIR}/metrics/Metric.h
    ${DIAT_LIB_DIR}/metrics/EuclideanMetric.h 
    ${DIAT_LIB_DIR}/metrics/SquaredEuclideanMetric.h 
    ${DIAT_LIB_DIR}/utils/Geometry.h
    ${DIAT_LIB_DIR}/io/IO.h
    ${DIAT_LIB_DIR}/datastructures/heap/Heap.h
    ${DIAT_LIB_DIR}/datastructures/heap/MaxHeap.h
    ${DIAT_LIB_DIR}/datastructures/heap/MinHeap.h
    ${DIAT_LIB_DIR}/linalg/Linalg.h
    ${DIAT_LIB_DIR}/linalg/SVD.h
    ${DIAT_LIB_DIR}/linalg/RandomSVD.h
    ${DIAT_LIB_DIR}/linalg/RandomRange.h
    ${DIAT_LIB_DIR}/linalg/LapackDefs.h
    ${DIAT_LIB_DIR}/linalg/SymmetricEigensystem.h
    ${DIAT_LIB_DIR}/linalg/Matrix.h
    ${DIAT_LIB_DIR}/linalg/Vector.h 
    ${DIAT_LIB_DIR}/linalg/DenseVector.h
    ${DIAT_LIB_DIR}/linalg/DenseMatrix.h
    ${DIAT_LIB_DIR}/linalg/SymmetricDenseMatrix.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/Gyroscope.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/ShmGyroscope.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/Font.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/GlutStrokeFont.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/Animation.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/Animator.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/Data.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/Display.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/DisplayElement.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/DummyAnimation.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/LabelsDEL.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/PatchDEL.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/PCADEL.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/ProfileDEL.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/ProjectionDEL.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/DensityProjectionDEL.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/GMRADensityProjectionDEL.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/Rotation.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/SearchDEL.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/TextDEL.h
    ${DIAT_SOURCE_DIR}/tools/gyroscope/colormapper.h
    ./gyroscope/src
)

ADD_CUSTOM_TARGET(R_gyroscope_build  
    COMMAND R CMD build --resave-data gyroscope)
ADD_DEPENDENCIES(R_gyroscope_build R_gyroscope)

ADD_CUSTOM_TARGET(R_gyroscope_install
    COMMAND R CMD INSTALL --no-multiarch gyroscope)
ADD_DEPENDENCIES(R_gyroscope_install R_gyroscope)

ADD_CUSTOM_TARGET(R_gyroscope_check 
    COMMAND R CMD check --no-multiarch gyroscope)
ADD_DEPENDENCIES(R_gyroscope_check R_gyroscope)




################ gyroscopeV2 #################

ADD_CUSTOM_TARGET(R_gyroscopeV2
    COMMAND rsync -aC ${Rcode_SOURCE_DIR}/gyroscopeV2 .
    COMMAND rsync -aC ${DIAT_LIB_DIR}/ann/ ./gyroscopeV2/src
    COMMAND rsync -aC ${DIAT_LIB_DIR}/Eigen/ ./gyroscopeV2/src/Eigen
    COMMAND rsync -aC ${DIAT_LIB_DIR}/gmra/ ./gyroscopeV2/src 
    COMMAND cp 
    ${DIAT_LIB_DIR}/metrics/EigenMetric.h
    ${DIAT_LIB_DIR}/metrics/EigenEuclideanMetric.h 
    ${DIAT_LIB_DIR}/metrics/EigenL1Metric.h 
    ${DIAT_LIB_DIR}/metrics/EigenSquaredEuclideanMetric.h 
    ${DIAT_LIB_DIR}/linalg/EigenRandomSVD.h
    ${DIAT_LIB_DIR}/linalg/EigenRandomRange.h
    ${DIAT_LIB_DIR}/statistics/EigenWasserstein.h
    ${DIAT_LIB_DIR}/statistics/Random.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/ShmGyroscope.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/Gyroscope.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/Font.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/GlutStrokeFont.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/Animation.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/Animator.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/DataInterface.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/StandardData.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/GMRAData.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/SingleGMRAData.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/MultiGMRAData.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/ClusteredData.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/Display.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/DisplayElement.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/DummyAnimation.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/ImageDEL.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/PatchDEL.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/PCADEL.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/ProfileDEL.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/ProjectionDEL.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/Rotation.h
    ${DIAT_SOURCE_DIR}/tools/gyroscopeV2/colormapper.h
    ./gyroscopeV2/src
)

ADD_CUSTOM_TARGET(R_gyroscopeV2_build  
    COMMAND R CMD build --resave-data gyroscopeV2)
ADD_DEPENDENCIES(R_gyroscopeV2_build R_gyroscopeV2)

ADD_CUSTOM_TARGET(R_gyroscopeV2_install
    COMMAND R CMD INSTALL --no-multiarch gyroscopeV2)
ADD_DEPENDENCIES(R_gyroscopeV2_install R_gyroscopeV2)

ADD_CUSTOM_TARGET(R_gyroscopeV2_check 
    COMMAND R CMD check --no-multiarch gyroscopeV2)
ADD_DEPENDENCIES(R_gyroscopeV2_check R_gyroscopeV2)




################ focus #################

ADD_CUSTOM_TARGET(R_focus
    COMMAND rm -rf ./focus
    COMMAND rsync -aC ${Rcode_SOURCE_DIR}/focus .
    COMMAND cp -a ${DIAT_LIB_DIR}/ann/ANN ./focus/src/ANN
    COMMAND cp -a ${DIAT_LIB_DIR}/Eigen ./focus/src/Eigen
    COMMAND cp -a ${DIAT_LIB_DIR}/gmra/* ./focus/src/. 
    COMMAND cp  ${DIAT_SOURCE_DIR}/tools/focus/*.h ./focus/src/. 
    COMMAND cp  ${DIAT_SOURCE_DIR}/tools/focus/projectors/* ./focus/src/.
    COMMAND cp  ${DIAT_SOURCE_DIR}/tools/focus/dels/* ./focus/src/. 
    COMMAND cp  
    ${DIAT_LIB_DIR}/metrics/EigenMetric.h
    ${DIAT_LIB_DIR}/metrics/EigenEuclideanMetric.h 
    ${DIAT_LIB_DIR}/metrics/EigenL1Metric.h 
    ${DIAT_LIB_DIR}/metrics/EigenSquaredEuclideanMetric.h 
    ${DIAT_LIB_DIR}/linalg/EigenRandomSVD.h
    ${DIAT_LIB_DIR}/linalg/EigenRandomRange.h
    ${DIAT_LIB_DIR}/statistics/EigenWasserstein.h
    ${DIAT_LIB_DIR}/statistics/Random.h
    ./focus/src
)

ADD_CUSTOM_TARGET(R_focus_build  
    COMMAND R CMD build --resave-data focus)
ADD_DEPENDENCIES(R_focus_build R_focus)

ADD_CUSTOM_TARGET(R_focus_install
    COMMAND R CMD INSTALL --no-multiarch focus)
ADD_DEPENDENCIES(R_focus_install R_focus)

ADD_CUSTOM_TARGET(R_focus_check 
    COMMAND R CMD check --no-multiarch focus)
ADD_DEPENDENCIES(R_focus_check R_focus)







############## gmra package ##############


ADD_CUSTOM_TARGET(R_gmra 
    COMMAND rsync -aC ${Rcode_SOURCE_DIR}/gmra . 
    COMMAND rsync -aC ${DIAT_LIB_DIR}/gmra/ ./gmra/src 
    COMMAND rsync -aC ${DIAT_LIB_DIR}/Eigen/ ./gmra/src/Eigen 
    COMMAND cp 
    ${DIAT_LIB_DIR}/metrics/EigenMetric.h
    ${DIAT_LIB_DIR}/metrics/EigenEuclideanMetric.h 
    ${DIAT_LIB_DIR}/metrics/EigenL1Metric.h 
    ${DIAT_LIB_DIR}/metrics/EigenSquaredEuclideanMetric.h 
    ${DIAT_LIB_DIR}/linalg/EigenRandomSVD.h
    ${DIAT_LIB_DIR}/linalg/EigenRandomRange.h
    ${DIAT_LIB_DIR}/statistics/Random.h
    ${DIAT_LIB_DIR}/statistics/EigenWasserstein.h
    ${DIAT_LIB_DIR}/statistics/Kmeans.h
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


######## proceduraldata package #########

ADD_CUSTOM_TARGET(R_proceduraldata 
    COMMAND rsync -aC ${Rcode_SOURCE_DIR}/proceduraldata . 
)

ADD_CUSTOM_TARGET(R_proceduraldata_build  
    COMMAND R CMD build --resave-data proceduraldata)
ADD_DEPENDENCIES(R_proceduraldata_build R_proceduraldata)

ADD_CUSTOM_TARGET(R_proceduraldata_install
    COMMAND R CMD INSTALL proceduraldata)
ADD_DEPENDENCIES(R_proceduraldata_install R_proceduraldata)

ADD_CUSTOM_TARGET(R_proceduraldata_check 
    COMMAND R CMD check proceduraldata)
ADD_DEPENDENCIES(R_proceduraldata_check R_proceduraldata)




  


############## multiscale persistent homology package ##############

ADD_CUSTOM_TARGET(R_mph
    COMMAND rsync -aC ${Rcode_SOURCE_DIR}/mph . 
    COMMAND rsync -aC ${DIAT_LIB_DIR}/gmra/ ./mph/src 
    COMMAND rsync -aC ${DIAT_LIB_DIR}/Eigen/ ./mph/src/Eigen 
    COMMAND cp 
    ${DIAT_LIB_DIR}/topology/Simplex.h
    ${DIAT_LIB_DIR}/topology/MappingRips.h
    ${DIAT_LIB_DIR}/topology/MappingRips2.h
    ${DIAT_LIB_DIR}/topology/MultiscaleMappingRips.h
    ${DIAT_LIB_DIR}/topology/MappingPersistentHomology.h
    ${DIAT_LIB_DIR}/metrics/EigenMetric.h
    ${DIAT_LIB_DIR}/metrics/EigenL1Metric.h
    ${DIAT_LIB_DIR}/metrics/EigenEuclideanMetric.h 
    ${DIAT_LIB_DIR}/metrics/EigenSquaredEuclideanMetric.h 
    ${DIAT_LIB_DIR}/linalg/EigenRandomSVD.h
    ${DIAT_LIB_DIR}/linalg/EigenRandomRange.h
    ${DIAT_LIB_DIR}/statistics/EigenWasserstein.h
    ${DIAT_LIB_DIR}/statistics/Random.h
    ${DIAT_LIB_DIR}/io/LinalgIO.h
    ./mph/src
)

ADD_CUSTOM_TARGET(R_mph_build  
    COMMAND R CMD build --resave-data mph)
ADD_DEPENDENCIES(R_mph_build R_mph)

ADD_CUSTOM_TARGET(R_mph_install
    COMMAND R CMD INSTALL  --no-multiarch   mph)
ADD_DEPENDENCIES(R_mph_install R_mph)

ADD_CUSTOM_TARGET(R_mph_check 
    COMMAND R CMD check --no-multiarch mph)
ADD_DEPENDENCIES(R_mph_check R_mph)



############## multiscale optimal tranpsort package ##############

ADD_CUSTOM_TARGET(R_mop
    COMMAND rsync -aC ${Rcode_SOURCE_DIR}/mop . 
    COMMAND rsync -aC ${DIAT_LIB_DIR}/gmra/ ./mop/src 
    COMMAND rsync -aC ${DIAT_LIB_DIR}/Eigen/ ./mop/src/Eigen 
    COMMAND cp 
    ${DIAT_LIB_DIR}/transport/SinkhornTransport.h
    ${DIAT_LIB_DIR}/transport/SparseSinkhornTransport.h
    ${DIAT_LIB_DIR}/transport/MultiscaleTransport.h
    ${DIAT_LIB_DIR}/transport/MultiscaleTransportLP.h
    ${DIAT_LIB_DIR}/transport/NeighborhoodStrategy.h
    ${DIAT_LIB_DIR}/transport/ExpandNeighborhoodStrategy.h
    ${DIAT_LIB_DIR}/transport/RefineNeighborhoodStrategy.h
    ${DIAT_LIB_DIR}/transport/PotentialNeighborhoodStrategy.h
    ${DIAT_LIB_DIR}/transport/PropagationStrategy.h
    ${DIAT_LIB_DIR}/transport/NeighborhoodPropagationStrategy.h
    ${DIAT_LIB_DIR}/transport/CapacityPropagationStrategy.h
    ${DIAT_LIB_DIR}/transport/IteratedCapacityPropagationStrategy.h
    ${DIAT_LIB_DIR}/transport/RandomizedNeighborhoodPropagationStrategy.h
    ${DIAT_LIB_DIR}/transport/MaxEntropyPropagationStrategy.h
    ${DIAT_LIB_DIR}/transport/SinkhornPropagationStrategy.h
    ${DIAT_LIB_DIR}/transport/TransportLP.h
    ${DIAT_LIB_DIR}/transport/TransportLPSolver.h
    ${DIAT_LIB_DIR}/transport/GLPKSolver.h
    ${DIAT_LIB_DIR}/transport/MOSEKSolver.h
    ${DIAT_LIB_DIR}/transport/CPLEXSolver.h
    ${DIAT_LIB_DIR}/transport/CPLEXNetworkSolver.h
    ${DIAT_LIB_DIR}/transport/LPSolver.h
    ${DIAT_LIB_DIR}/metrics/EigenMetric.h
    ${DIAT_LIB_DIR}/metrics/EigenEuclideanMetric.h 
    ${DIAT_LIB_DIR}/metrics/EigenSquaredEuclideanMetric.h 
    ${DIAT_LIB_DIR}/metrics/EigenL1Metric.h 
    ${DIAT_LIB_DIR}/statistics/Random.h
    ${DIAT_LIB_DIR}/statistics/EigenWasserstein.h
    ${DIAT_LIB_DIR}/linalg/EigenRandomSVD.h
    ${DIAT_LIB_DIR}/linalg/EigenRandomRange.h
    ${DIAT_LIB_DIR}/mosek/mosek.h
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
    COMMAND rsync -aC ${Rcode_SOURCE_DIR}/sinkhorn . 
    COMMAND rsync -aC ${DIAT_LIB_DIR}/gmra/ ./sinkhorn/src 
    COMMAND rsync -aC ${DIAT_LIB_DIR}/Eigen/ ./sinkhorn/src/Eigen 
    COMMAND cp 
    ${DIAT_LIB_DIR}/transport/SinkhornTransport.h
    ${DIAT_LIB_DIR}/transport/SparseSinkhornTransport.h
    ${DIAT_LIB_DIR}/transport/MultiscaleTransport.h
    ${DIAT_LIB_DIR}/transport/MultiscaleSinkhornTransport.h
    ${DIAT_LIB_DIR}/metrics/EigenMetric.h
    ${DIAT_LIB_DIR}/metrics/EigenEuclideanMetric.h 
    ${DIAT_LIB_DIR}/metrics/EigenSquaredEuclideanMetric.h 
    ${DIAT_LIB_DIR}/metrics/EigenL1Metric.h 
    ${DIAT_LIB_DIR}/statistics/Random.h
    ${DIAT_LIB_DIR}/statistics/EigenWasserstein.h
    ${DIAT_LIB_DIR}/linalg/EigenRandomSVD.h
    ${DIAT_LIB_DIR}/linalg/EigenRandomRange.h
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





################ install all gmra related packages ###############

ADD_CUSTOM_TARGET(R_gmra_suite_build) 
ADD_DEPENDENCIES(R_gmra_suite_build R_gmra_build R_mop_build R_mph_build R_sinkhorn_build)

ADD_CUSTOM_TARGET(R_gmra_suite_install) 
ADD_DEPENDENCIES(R_gmra_suite_install R_gmra_install R_mop_install R_mph_install R_sinkhorn_install)











############## multisacle entropy package ##############


ADD_CUSTOM_TARGET(R_mmsentropy 
    COMMAND rsync -aC ${Rcode_SOURCE_DIR}/mmsentropy . 
    COMMAND rsync -aC ${DIAT_LIB_DIR}/Eigen/ ./mmsentropy/src/Eigen 
    COMMAND cp 
    ${DIAT_LIB_DIR}/statistics/MultiscaleEntropy.h
    ./mmsentropy/src
)

ADD_CUSTOM_TARGET(R_mmsentropy_build  
    COMMAND R CMD build --resave-data mmsentropy)
ADD_DEPENDENCIES(R_mmsentropy_build R_mmsentropy)

ADD_CUSTOM_TARGET(R_mmsentropy_install
    COMMAND R CMD INSTALL  --no-multiarch mmsentropy)
ADD_DEPENDENCIES(R_mmsentropy_install R_mmsentropy)

ADD_CUSTOM_TARGET(R_mmsentropy_check 
    COMMAND R CMD check --no-multiarch mmsentropy)
ADD_DEPENDENCIES(R_mmsentropy_check R_mmsentropy)



