#ifndef BINARYGMRAPREDICTOR_H
#define BINARYGMRAPREDICTOR_H

#include "GMRATree.h"

#include <vector>


template <typename TPrecision>
class BinaryGMRAPredictor<TPrecision>{

  private:
    int nLabeled;

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;




    //Expect 0, 1 labels and -1 for missing labels
    std::vector<double> predict(GMRATree<TPrecision> &tree, std::vector< int > &labels){

    }; 

    MatrixXp getPredictions(){

    };

    MatrixXp getLOO(){

    };


};


#endif
