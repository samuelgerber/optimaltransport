#ifndef GMRADATAOBJECT_H
#define GMRADATAOBJECT_H


#include <Eigen/Dense>




template <typename TPrecision>
class GMRADataObject{
  public:    
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;

    virtual ~GMRADataObject(){};

    virtual VectorXp getPoint(int i) = 0;

    virtual int numberOfPoints() = 0;

    virtual int dimension() = 0;

    virtual TPrecision getMass(int i){
      return 1.0/this->numberOfPoints();
    };


};



template <typename TPrecision>
class MatrixGMRADataObject : public GMRADataObject<TPrecision>{
  public:    
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;

  private:
    MatrixXp X;

  public:    

    MatrixGMRADataObject(const MatrixXp &data):X(data){
    };

    virtual VectorXp getPoint(int i){
      return X.col(i);
    };

    virtual int numberOfPoints(){
      return X.cols();
    };

    virtual int dimension(){
      return X.rows();
    };

};

template <typename TPrecision>
class WeightedMatrixGMRADataObject : public MatrixGMRADataObject<TPrecision>{
  public:    
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;

  private:
    VectorXp weights;

  public:    

    WeightedMatrixGMRADataObject(MatrixXp &data, VectorXp &w):MatrixGMRADataObject<TPrecision>(data), weights(w){
    };

    virtual TPrecision getMass(int i){
      return weights(i);
    };

};






#endif
