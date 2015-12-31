#ifndef GMRASPATIALVARIANCEREFINE_H
#define GMRASPATIALVARIANCEREFINE_H

#include "GMRAVisitor.h"




template <typename TPrecision>
class GMRAJointVariancePrune : public Visitor<TPrecision>{
  
  public:    
  
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
  
  private:
  
    MatrixXp &X;  
    int nCor;
    int nSpatial;
    TPrecision radius;


  public:
    
    GMRASpatialVarianceRefine(MatrixXp &data, int nC, int nS, TPrecision
        maxRadius) : X(data), nCor(nC), nSpatial(nS){
   
    };


    virtual void visit(GMRANode<TPrecision> *node){
          
    };

};






#endif
