
#ifndef EIGENWASSERSTEIN_H
#define EIGENWASSERSTEIN_H

#include <Eigen/Dense>
#include <Eigen/SVD>
#include "EigenSquaredEuclideanMetric.h"

#include <cmath>

template <typename TPrecision>
class Wasserstein{

  private:
    SquaredEuclideanMetric<TPrecision> metric;
  
  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
  
 
    Wasserstein(){}; 
    ~Wasserstein(){};

    TPrecision distance(const MatrixXp &U1, const VectorXp &S1, const VectorXp &C1, 
                        const MatrixXp &U2, const VectorXp &S2, const VectorXp &C2){ 
      TPrecision md = metric.distance(C1, C2);
      TPrecision cd = covarianceDistSquared(U1, S1, U2, S2);
      return sqrt(md + cd);
    
    };
   

    //U eigenvectors, S variances of covariance matrices
    TPrecision covarianceDistSquared(const MatrixXp &U1, const VectorXp &S1, 
                                     const MatrixXp &U2, const VectorXp &S2){
      using namespace Eigen; 
      
      TPrecision ts1 = S1.sum();
      TPrecision ts2 = S2.sum();

      //std::cout << U1.M() << " x " << U1.N() << std::endl;
      //std::cout << U2.M() << " x " << U2.N() << std::endl;
      TPrecision tc = 0;
      if(U1.cols() == 0 || U2.cols() == 0 || U1.rows() == 0 || U2.rows() == 0){
        //nothing todo tc = 0
      }
      else{
        MatrixXp U1tU2 = U1.transpose() * U2;
        MatrixXp U2tU1 = U1tU2.transpose();
        for(int i=0; i< S2.size(); i++){
          U1tU2.col(i) *= S2(i);
        }

        //std::cout << U1tU2.M() << " x " << U1tU2.N() << std::endl;
        //std::cout << U2tU1.M() << " x " << U2tU1.N() << std::endl << std::endl;

        MatrixXp A = U1tU2 * U2tU1;
        for(int i=0; i< S1.size(); i++){
          TPrecision s= sqrt( S1(i) );
          A.col(i) *= s;
          A.row(i) *= s;
        }

        JacobiSVD<MatrixXp> svd(A, ComputeThinU | ComputeThinV);
        for(int i=0; i<svd.singularValues().size(); i++){
          tc += sqrt(svd.singularValues()(i));
        }
        
      }

      return ts1 + ts2 - 2*tc;
        
    };








    TPrecision covarianceDistSquared(const MatrixXp &Cov1, const MatrixXp &Cov2){ 
      using namespace Eigen; 

      JacobiSVD<MatrixXp> svd1(Cov1, ComputeThinU | ComputeThinV);
      JacobiSVD<MatrixXp> svd2(Cov2, ComputeThinU | ComputeThinV);

      MatrixXp U1;
      MatrixXp U2;
      if( svd1.U.rows() < svd2.matrixU().rows() ){
        U2 = svd2.matrixU();
        U1 = MatrixXp( U2.rows(), svd1.matrixU().cols() );
        U1.topLeftCorner( svd1.matrixU().rows(),  svd1.matrixU().cols() ) = svd1.matrixU();
      }
      else if( svd2.U.M() < svd1.U.M() ){
        U1 = svd1.matrixU();
        U2 = MatrixXp( U1.rows(), svd2.matrixU().cols() );
        U2.topLeftCorner( svd2.matrixU().rows(),  svd2.matrixU().cols() ) = svd2.matrixU();
      }
  
      TPrecision d = distance2(U1, svd1.singularValues(), U2, svd2.singularValues() );
     
      return d;
    };







    //U eigenvectors, S variances of covariance matrices
    MatrixXp covarianceMap(const MatrixXp &U1, const VectorXp &S1, 
                           const MatrixXp &U2, const VectorXp &S2){
      
      TPrecision ts1 = S1.sum();
      TPrecision ts2 = S2.sum();

      //std::cout << U1.M() << " x " << U1.N() << std::endl;
      //std::cout << U2.M() << " x " << U2.N() << std::endl;
      TPrecision tc = 0;
      if(U1.cols() == 0 || U2.cols() == 0 || U1.rows() == 0 || U2.rows() == 0){
        //use identity map
        MatrixXp T = MatrixXp::Zero(S2.size(), S1.size());
        for(int i=0; i < std::min(T.cols(), T.rows()); i++){
          T(i,i) = 1;
        }
        return T;
      }
      else{
        MatrixXp U1tU2 = U1.transpose() * U2;
        MatrixXp U2tU1 = U1tU2.transpose();
        for(int i=0; i< S2.size(); i++){
          U1tU2.col(i ) *=  S2(i);
        }

        //std::cout << U1tU2.M() << " x " << U1tU2.N() << std::endl;
        //std::cout << U2tU1.M() << " x " << U2tU1.N() << std::endl << std::endl;

        MatrixXp A = U1tU2 * U2tU1;
        for(int i=0; i< S1.size(); i++){
          TPrecision s= sqrt(S1(i));
          A.col(i) *= s;
          A.row(i) *= s;
        }


        MatrixXp U1t = U1.transpose();
        for(int i=0; i< S1.size(); i++){
          U1t.row(i) /= sqrt(S1(i));
        }
        MatrixXp C1 = U1 * U1t;
        MatrixXp T = (C1 * A) * C1;

        return T;
      }
    };



    
};
#endif
