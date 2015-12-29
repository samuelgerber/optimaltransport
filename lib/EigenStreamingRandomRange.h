#ifndef EIGENSTREAMINGRANDOMRANGE_H
#define EIGENSTREAMINGRANDOMRANGE_H

#include "Random.h"
#include "EigenRandomRange.h"
#include <Eigen/QR>

namespace EigenLinalg{
  
class<typename TPrecision>
class RowStreamingRandomRange : public RandomRange{
  
  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
    

  private:
    MatrixXp N;
    MatrixXp Y;
    MatrixXp Q;
    int index;

  public:

   

    RowStreamingRandomRange(int m, int n, int d){
      using namespace Eigen;

      N = MatrixXp(n, d);      
      static Random<double> rand;
      for(unsigned int i=0; i< N.rows(); i++){
        for(unsigned int j=0; j< N.cols(); j++){
          N(i, j) = rand.Normal();
        }
      }

      Y = MatrixXp(m, d);
      index = 0;
    };


    void Add(VectorXp &x){
     using namespace Eigen;
     Y.row(index) =   x.transpose() * N;
     index++;
    }; 


    void Compute(){
     using namespace Eigen;
      HouseholderQR<MatrixXp> qr(Y);
      Q = qr.householderQ() * MatrixXp::Identity( Y.rows(), N.cols() );   
    };


    MatrixXp &GetRange(){
      return Q;
    };

    MatrixXp &GetRandomProjector(){
      return N;
    };
    
    MatrixXp &GetProjected(){
      return Y;
    };

};

}

#endif 
