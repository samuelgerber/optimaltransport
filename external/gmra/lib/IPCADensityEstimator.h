#ifndef IPCADENSITYESTIMATOR_H
#define IPCADENSITYESTIMATOR_H

#include "GMRADensityEstimator.h"
#include "IPCATree.h"

#include <vector>

#include <iostream>




//Tree data structure
template <typename TPrecision>
class IPCADensityEstimator : public GMRADensityEstimator<TPrecision>{


public:


    IPCADensityEstimator(IPCATreeFactory<TPrecision> &c):config(c){};

    //X an m * n matrix of n  m-dimensional data points stored in column major
    //order
    virtual void construct(TPrecision *X, int m, int n, int nRuns){
      for(int i=0; i<nRuns; i++){
        IPCATree<TPrecision> *tree = new IPCATree<TPrecision>(config);
        tree->construct(X, m, n); 
        this->addTree(tree);
      }
      
    };
  


};




#endif
