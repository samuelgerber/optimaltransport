#ifndef NODEDISTANCE_H
#define NODEDISTANCE_H

#include <Eigen/Dense>


#include "GMRATree.h"
#include "EigenMetric.h"

template <typename TPrecision>
class GMRANode;


template < typename TPrecision >
class NodeDistance{

  public:
    virtual ~NodeDistance(){};

    virtual TPrecision distance(GMRANode<TPrecision> *n1, GMRANode<TPrecision> *n2) = 0;

};




template < typename TPrecision >
class CenterNodeDistance : public NodeDistance<TPrecision> {
  private:
    Metric<TPrecision> *metric;

  public: 

    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
    CenterNodeDistance(Metric<TPrecision> *m):metric(m){};

    ~CenterNodeDistance(){
    };

    TPrecision distance(GMRANode<TPrecision> *n1, GMRANode<TPrecision> *n2){
      VectorXp &x1 = n1->getCenter();
      VectorXp &x2 = n2->getCenter();
      return metric->distance(x1, x2);
    };
};




#endif
