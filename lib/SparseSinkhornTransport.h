

#ifndef SPARSESINKHORNTRANSPORT_H
#define SPARSESINKHORNTRANSPORT_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iostream>



class SparseSinkhornTransport {


  private:
  
    Eigen::MatrixXd lScaling;
    Eigen::MatrixXd rScaling;
    Eigen::SparseMatrix<double> K;
    Eigen::VectorXd d;

  
  public:

    SparseSinkhornTransport(){
    };


    void initalizeLeftScaling(Eigen::MatrixXd &left){
      lScaling = left;
    };

    void transport(Eigen::VectorXd &mu, Eigen::MatrixXd &nu,
        Eigen::SparseMatrix<double> &Kin, Eigen::SparseMatrix<double> &U, double
        tol = 0.00001, int maxIter=1000){ 
      using namespace Eigen;
      
      K = Kin;

      SparseMatrix<double> muK = K;
      SparseMatrix<double> Kt = K.transpose();
      //= K.array().colwise() / mu.array();
      for (int k=0; k<muK.outerSize(); ++k){
        for (SparseMatrix<double>::InnerIterator it(muK,k); it; ++it) {
          it.valueRef() = it.value() / mu( it.row() );
        }
      }

      if( lScaling.rows() != mu.size() ){
        lScaling = MatrixXd::Constant(mu.size(), nu.cols(), 1.0 /  mu.size() );
      }


      d = VectorXd::Constant(nu.cols(), -1);

      for(int i=0; i< maxIter; i++){
        rScaling = nu.array() / ( Kt * lScaling ).array();
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



    Eigen::SparseMatrix<double> getTransportPlan(int i=0){
      using namespace Eigen;
      SparseMatrix<double> T = K;
      //= K.array().colwise() / mu.array();
      for (int k=0; k<T.outerSize(); ++k){
        for (SparseMatrix<double>::InnerIterator it(T,k); it; ++it) {
          it.valueRef() *=  lScaling( it.row(), i ) * rScaling( it.col(), i );
        }
      }
      T.makeCompressed();
      return T; 
    };


    Eigen::MatrixXd getLeftScaling(){
      return lScaling;
    };
 


};


#endif
