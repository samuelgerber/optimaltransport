#ifndef GMRADENSITYESTIMATOR_H
#define GMRADENSITYESTIMATOR_H

#include <Eigen/Dense>
#include "GMRATree.h"


#include <vector>
#include <set>

#include <cmath>
#include <limits>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

template <typename TPrecision>
class NodeDensityEstimator {
  
  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
 
    std::vector<TPrecision> estimateDensity(GMRATree<TPrecision> &gmra, const MatrixXp &X){
      std::vector<TPrecision> p( X.cols() );
      for(int i=0;i<X.cols(); i++){
        p[i] = density(gmra, X.col(i) );
      }

      return p;
    };

    virtual TPrecision density(GMRATree<TPrecision> &gmra, const VectorXp &x) = 0;

};

template <typename TPrecision>
class KernelNodeDensityEstimator : public NodeDensityEstimator<TPrecision>{
  
  private:
    int minP, maxP;
    TPrecision sigma;
  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;

    KernelNodeDensityEstimator(int minPoints, int maxPoints, TPrecision s ){
      minP = minPoints;
      maxP = maxPoints;
      sigma = s;
    };





    TPrecision density(GMRATree<TPrecision> &gmra, const VectorXp &x){

      //TPrecision nPoints = gmra.getRoot()->getPoints().size();

      std::vector<GMRANode<TPrecision> *> path = gmra.getLeafPath(x);
      TPrecision p = -1;
      for(int i=0; i<path.size(); i++){
        GMRANode<TPrecision> *node = path[i]; 
        if(node->getPoints().size() < maxP){
          GMRANode<TPrecision> *parent = node->getParent();
          if( parent == NULL || parent->getPoints().size() > maxP ){
            std::vector<int> ind = node->getPoints();
            if(ind.size() > minP){
              p = 0;
              int dim = node->getIntrinsicDimension();
              TPrecision trace = pow(sigma, dim);
              TPrecision norm = 1.0 /( sqrt(trace) * pow(2*M_PI, dim/2));
              for(int i=0; i<ind.size(); i++){
                const VectorXp &x2 = gmra.getPoint( ind[i] );
                TPrecision d = (x - x2).squaredNorm();
                p +=  norm * exp(- d /(2* pow(sigma, dim) ) );
              }
              p /= ind.size();
            }
            
            break;
          }
        }
      }

      return p;
    };


};



template <typename TPrecision>
class GMRADensityEstimator{


private:
  
  std::vector< GMRATree<TPrecision> *> trees;


public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;

   
    virtual ~GMRADensityEstimator(){};
    
    virtual void addTree(GMRATree<TPrecision> *tree){
      trees.push_back(tree);
    };

    int getNumberOfTrees(){
      return trees.size();
    };

    GMRATree<TPrecision> *getTree(int i){
      return trees[i];
    };

    virtual std::vector<TPrecision> estimate(const MatrixXp &X,
        NodeDensityEstimator<TPrecision> &estimator ) = 0;
};

template <typename TPrecision>
class MeanGMRADensityEstimator : public GMRADensityEstimator<TPrecision>{

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;

    std::vector<TPrecision> estimate(const MatrixXp &X,
        NodeDensityEstimator<TPrecision> &estimator ) { 

      std::vector<TPrecision> p(X.cols(), 0);
      std::vector<TPrecision> counts(X.cols(), 0);

      for(int i=0; i<this->getNumberOfTrees(); i++){
        std::vector<TPrecision> tmp =  estimator.estimateDensity(*this->getTree(i), X);
        for(int j=0; j<p.size(); j++){
          if(tmp[j] >= 0){
            p[j] += tmp[j];
            counts[j] += 1;
          }
        }
      }

      for(int j=0; j<p.size(); j++){
        p[j] /= counts[j];
      }

      return p;

    };

};


template <typename TPrecision>
class MaxGMRADensityEstimator : public GMRADensityEstimator<TPrecision>{

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;

    std::vector<TPrecision> estimate(const MatrixXp &X,
        NodeDensityEstimator<TPrecision> &estimator ) { 

      std::vector<TPrecision> p(X.cols(), 0);

      for(int i=0; i<this->getNumberOfTrees(); i++){
        std::vector<TPrecision> tmp =  estimator.estimateDensity(*this->getTree(i), X);
        for(int j=0; j<p.size(); j++){
          if(tmp[j] >= 0){
            p[j] = std::max(p[j], tmp[j]);
          }
        }
      }

      return p;

    };

};







#endif
