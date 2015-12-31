#ifndef GMRAPREDICTOR_H
#define GMRAPREDICTOR_H

#include "GMRATree.h"
#include "IKMNode.h"
#include "IKMKmeansData.h"
#include "Kmeans.h"


#include <map>
#include <iostream>
#include <fstream>

#include <set>


template <typename TPrecision, typename TLabel>
class GMRAPredictor<TPrecision>{

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
    



    virtual void predict(GMRATree<TPrecision> &tree, std::vector< Label<TLabel> > &labels) = 0; 

};


#endif
