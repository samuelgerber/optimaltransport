#ifndef EIGENRANDOMSVD_H
#define EIGENRANDOMSVD_H

#include "EigenRandomRange.h"
#include <Eigen/SVD>

namespace EigenLinalg{

  
  
template <typename TPrecision>
class RandomSVD{    
  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
    
  private:
    MatrixXp U;
    MatrixXp B;
    VectorXp S;

  public:

    RandomSVD(){
  
    };
    
   RandomSVD(Eigen::Ref<MatrixXp> Xin, int d, int nPowerIt = 0){

      using namespace Eigen;
      using namespace EigenLinalg;

      StandardRandomRange<TPrecision> range(Xin, d, nPowerIt);
      MatrixXp &Q = range.GetRange();
      B = Q.transpose() * Xin;

      JacobiSVD<MatrixXp> svd(B, ComputeThinU | ComputeThinV);
      S = svd.singularValues();
      U = Q * svd.matrixU();

    };
  

    MatrixXp &GetU(){
      return U;
    };
    
    VectorXp &GetS(){
      return S;
    };

    MatrixXp &GetProjected(){
      return B;
    };
      
    
};


}


#endif 
