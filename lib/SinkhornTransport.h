

#ifndef SINKHORNTRANSPORT_H
#define SINKHORNTRANSPORT_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iostream>



class SinkhornTransport {


  private:
  
    Eigen::MatrixXd lScaling;
    Eigen::MatrixXd rScaling;
    Eigen::MatrixXd K;
    Eigen::VectorXd d;

  
  public:

    SinkhornTransport(){
    };


    void transport(Eigen::VectorXd &mu, Eigen::MatrixXd &nu, Eigen::MatrixXd &Kin,
        Eigen::MatrixXd &U, double tol = 0.00001, int maxIter=1000){ 
      using namespace Eigen;
      K = Kin;

      //K = ( -lambda * D ).array().exp();
      //U = K.cwiseProduct( D );
      
      MatrixXd muK = K.array().colwise() / mu.array();
      lScaling = MatrixXd::Constant(mu.size(), nu.cols(), 1.0 /  mu.size() );


      //MatrixType Kt = K.transpose();
      d = VectorXd::Constant(nu.cols(), -1);

      for(int i=0; i< maxIter; i++){
        //lScaling = 1 /  ( muK * (nu.cwiseQuotient( Kt * lScaling ) ) ).array();
        rScaling = nu.array() / (K.transpose() * lScaling).array();
        lScaling = 1.0 / ( muK * rScaling ).array();
        VectorXd dTmp = ( lScaling.array() * (U * rScaling).array()).colwise().sum().transpose();
        if( ( 1.0 - dTmp.array() / d.array() ).abs().maxCoeff()  < tol){
          d = dTmp;
          break;
        }
        d = dTmp; 
      }

    }; 


    Eigen::VectorXd getDistances(){
      return d; 
    };


    Eigen::MatrixXd getTransportPlan(int i=0){
      using namespace Eigen;
      return ( K.array().colwise() * lScaling.col(i).array() ).rowwise() *
        rScaling.col(i).transpose().array(); 
    };

 


};


#endif
