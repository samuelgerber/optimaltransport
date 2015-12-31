#ifndef EIGENRANDOMRANGE_H
#define EIGENRANDOMRANGE_H

#include "Random.h"
#include <Eigen/QR>

namespace EigenLinalg{

template <typename TPrecision>
class RandomRange{
  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
  
    virtual ~RandomRange(){};

    virtual MatrixXp &GetRange() = 0;

};


template <typename TPrecision>
class StandardRandomRange : public RandomRange<TPrecision>{

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
  
  private:
    MatrixXp Q;

  public:


    StandardRandomRange(Eigen::Ref<MatrixXp> X, int
        d, int nPowerIt = 0){

      using namespace Eigen;



      static Random<double> rand;
      
      MatrixXp N(X.cols(), d);
      for(unsigned int i=0; i< N.rows(); i++){
        for(unsigned int j=0; j< N.cols(); j++){
          N(i, j) = rand.Normal();
        }
      }

      
      Q = X *N;
      HouseholderQR<MatrixXp> qr(Q);

      if(nPowerIt > 0){
        MatrixXp Z;
        for( int i=0; i<nPowerIt; i++){
          Z.noalias() = X.transpose() * qr.householderQ();

          qr.compute(Z);

          Q.noalias() = X * qr.householderQ();
          qr.compute(Q);
        }
      }
      
      Q = qr.householderQ() * MatrixXp::Identity( X.rows(), d );
    };
   
    MatrixXp &GetRange(){
      return Q;
    };

};




}

#endif 
