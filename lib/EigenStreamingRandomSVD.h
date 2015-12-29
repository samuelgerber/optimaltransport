#ifndef EIGENSTREAMINGRANDOMSVD_H
#define EIGENSTREAMINGRANDOMSVD_H

#include "EigenStreamingRandomRange.h"
#include <Eigen/SVD>

namespace EigenLinalg{

  
class<typename TPrecision>  
class RowStreamingRandomSVD{
  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
    
  private:
    
    MatrixXp Q;
    MatrixXp B;
  
    MatrixXp U;
    VectorXp S;

    int index;

  public:
    


    RowStreamingRandomSVD(RandomRange<TPrecision> &range, int n) {

      using namespace Eigen;

      Q = range.GetRange();

      B = MatrixXp(Q.cols(), n);
      B.setConstant(0);
      
      index = 0;
    };


    void Add(VectorXp &x){
      using namespace Eigen;

      VectorXp q = Q.row(index);
      for(int i=0; i<x.size(); i++){
        B.col(i) += q * x(i);
      }

      index++;
    };


    void Compute(){
      using namespace Eigen;

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
