#ifndef TRANSPORTBREGMAN_H
#define TRANSPORTBREGMAN_H

#include "LPSolver.h"

#include <map>



//Split-Bregman iterations for solving optimal transport
template <typename TPrecision>
class TransportBregman {

  

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
    
    
    TransportBregman() {
    };

    virtual ~TransportBregman(){
    };





    MatrixXp solve(const MatrixXp &C, const VectorXp &from, const VectorXp &to, TPrecision lambda = 1, maxIter = 10){

        MatrixXp B = MatrixXp::Zero( from.length(),  to.length() );
        MatrixXp Z = B;
        MatrixXp X = Z;

        TPrecision lambdaInv = 1.0/lambda;

        for(int iter=0; iter<maxIter; iter++){

          //solve step 1, quadratic problem with diagonal "regularization"
          for(int i=0; i<X.M(); i++){
            for(int j=0; j<X.N(); j++){
              X(i, j) = ( ( from(i) + to(j) ) + lambda * ( B(i, j) - Z(i, j) ) )/( 1.0 - lambda*C(i, j) );
            }
          }

          //solve step 2, shrinkage
          for(int i=0; i<X.M(); i++){
            for(int j=0; j<X.N(); j++){
              TPrecision tmp = C(i, j) * X(i, j) + B(i, j)
              Z(i, j) = std::copysign(1.0, tmp) * styd::max( fabs( tmp ) - lambdaInv, 0.0);
            }
          }

          //solve step 3
          for(int i=0; i<X.M(); i++){
            for(int j=0; j<X.N(); j++){
              B(i, j) = B(i,j) + C(i, j)*X(i, j) - Z(i, j);
            }
          }
        }

        return X;

      };




};


#endif


