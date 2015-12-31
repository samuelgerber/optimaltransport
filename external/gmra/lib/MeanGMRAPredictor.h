#ifndef MEANGMRAPREDICTOR_H
#define MEANGMRAPREDICTOR_H

#include "GMRATree.h"

#include <vector>

#include <cmath>

template <typename TPrecision>
class MeanGMRAPredictor  {

  private:
    GMRATree<TPrecision> &tree;
    std::vector< double > &signal;
    int minPoints;
  
  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;



    MeanGMRAPredictor(GMRATree<TPrecision> &t, std::vector<double> &s, int
        nStop) : tree(t), signal(s), minPoints(nStop) { 
    };




    std::vector<double> getPredictions(MatrixXp &X){
      std::vector<double> predictions( X.cols(), 0 );
      for(int i=0; i<X.cols(); i++){
        std::vector<GMRANode<TPrecision> *> path = tree.getLeafPath( X.col(i) );

        GMRANode<TPrecision> *node = path[0];
        for(int k=1; k<path.size(); k++){
          //TODO: check for number of labeled instances
          if(path[k]->getPoints().size() < minPoints){
            break;
          }
          node =path[k];
        }
        predictions[i] = getPrediction(node);
      }

      return predictions;
    };



  private:
    
    
    double getPrediction(GMRANode<TPrecision> *node){
      std::vector<int> &pts = node->getPoints();
      int nL=0;
      double sum=0;
      for(int i=0; i<pts.size(); i++){
        int index = pts[i];
        if( !isnan( signal[index] ) ){
          nL++;
          sum += signal[index];
        }
      }
      return sum/nL;
    };

};


#endif
