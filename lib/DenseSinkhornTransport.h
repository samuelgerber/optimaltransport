

#ifndef SINKHORNTRANSPORT_H
#define SINKHORNTRANSPORT_H

#include <Eigen/Dense>
#include <Eigen/Sparse>


template <typename MatrixType>
class SinkhornDistances{
   public:
     MatrixXd lScaling;
     MatrixXd rScaling;
     MatrixType K;
     MatrixType U;
     

     VectorXd getDistances(){
       MatrixXd d = lScaling .* (U * rScaling)
       return d.colSums();
     };

     MatrixType getTransportPlan(){
     };

};

template <typename MatrixType>
class SinkhornTransport {

  public:
    
    SinkhornTransport(){
    };


    SinkhornDistances<MatrixType> transport(VectorXd &mu, MatrixXd &nu,
        MatrixType &D, double lambda, int maxIter=1000){
      SinkhornDistances<MatrixType> dists;
      dists.K = exp(-lambda*D);
      dists.U = K .* M;
      dists.lScaling = ones(mu.length(), nu.ncol()) / u.length();
      
      MatrixType muK = K rowwise/ mu;
      VectorXd u = ones(mu.length())/mu.length();

      MatrixType Kt = K.transpose();
      VectorXd dOld(-1, nu.ncol());

      for(int i=0; i< maxIter; i++){
        u = 1 ./  ( muK * (nu./ (Kt * u) ) );

        dists.rScaling = b./(Kt * dists.lScaling);
        dists.lScaling = 1./(rK * dist.rScaling );
        VectorXd d = dists.getDistances();
        if( fabs( 1-(d ./ dOld).max() ) < tol){
          break;
        }
        dOld = d; 
      }

      return dists;

    };  





};


#endif
