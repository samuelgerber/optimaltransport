#ifndef GMRAVARIANCEPRUNE_H
#define GMRAVARIANCEPRUNE_H

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


  public:
    
    GMRAJointVariancePruneVisitor(MatrixXp &data, int nC, int nS) : X(data), nCor(nC), nSpatial(nS){
   

    };

    virtual void visit(GMRANode<TPrecision> *node){
      std::vector<int> pts;
      
      TPrecision varCor = 0;
      TPrecision varSpatial = 0;
      VectorXp mean = VectorXp::Zero(X.rows());

      std::vector<GMRANode<TPrecision> *> chitlums = node->getChildren();

      if(chitlums.size() == 0){
        pts = node->getPoints();
      }
      else{
        for(int i=0; i<chitlums.size(); i++){
          std::vector<int> &tmp = chitlums[i]->getPoints();
          pts.add(pts.end(), tmp.begin(), tmp.end() );
        }
      }

      for(int i = 0; i<  pts.size(); i++){
        mean += X.col(pits[i]);
      }
       
      mean.head( nCor ).normalize();
      mean.tail( nSpatial ) /= pts.size();
      
      VectorXp cMean = mean.head(nCor);
      VectorXp sMean = mean.tail(nSpatial);
    
      for(int i = 0; i<  pts.size(); i++){
        TPrecision tmp = cMean.dot( X.col(pts[i]).head(nCor) );
        varCor += tmp*tmp;

        varSpatial += (sMean - X.col(pts[i]).tail(nSpatial)).squaredNorm(); 
      }

      if( (varSpatial * lambda) > (varCor * factor) ){
      }
      else{
        node->setCenter(center);
        node->setPoints(
      } 
    };

};






#endif
