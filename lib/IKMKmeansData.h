#ifndef IKMKMEANSDATA_H
#define IKMKMEANSDATA_H

#include "GMRATree.h"
#include "Kmeans.h"


#include <map>
#include <iostream>
#include <fstream>

#include <set>


template <typename TPrecision>
class IKMKmeansDataFactory{
  public:    
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;

    virtual ~IKMKmeansDataFactory(){};

    virtual KmeansData<TPrecision> *createKmeansData(GMRADataObject<TPrecision> *d, std::vector<int> &sub) = 0;

};



template <typename TPrecision>
class L2GMRAKmeansData : public KmeansData<TPrecision>{
  private:
    GMRADataObject<TPrecision> *data;
    std::vector<int> subset;

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
    
    L2GMRAKmeansData(GMRADataObject<TPrecision> *d, std::vector<int> &sub) : data(d), subset(sub){
    };

    virtual VectorXp getPoint(int index){
      return data->getPoint( subset[index] );
    };

    virtual int getNumberOfPoints(){
      return subset.size();
    };
    
    VectorXp getMean( std::vector<int> &pts ){
      VectorXp mean = VectorXp::Zero(data->dimension());
      for(int i=0; i<pts.size(); i++){
        mean += getPoint( pts[i] );
      }
      mean /= pts.size();
      return mean;
    };

    
    virtual TPrecision getSquaredSimilarity(int index, const VectorXp &p){
      return (p -getPoint(index)).squaredNorm();
    };
      

};


template <typename TPrecision>
class L2GMRAKmeansDataFactory : public IKMKmeansDataFactory<TPrecision>{
  public:
    virtual KmeansData<TPrecision> *createKmeansData(GMRADataObject<TPrecision> *d, std::vector<int> &sub){
      return new L2GMRAKmeansData<TPrecision>(d, sub);
    };
};














template <typename TPrecision>
class CorrelationGMRAKmeansData : public KmeansData<TPrecision>{
  private:
    GMRADataObject<TPrecision> *data;
    std::vector<int> subset;

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
    
    CorrelationGMRAKmeansData(GMRADataObject<TPrecision> *d, std::vector<int> &sub) : data(d), subset(sub){
    };

    virtual VectorXp getPoint(int index){
      return data->getPoint( subset[index] );
    };

    virtual int getNumberOfPoints(){
      return subset.size();
    };
    
    VectorXp getMean( std::vector<int> &pts ){
      VectorXp mean = VectorXp::Zero(data->dimension());
      for(int i=0; i<pts.size(); i++){
        mean += getPoint( pts[i] );
      }
      mean.normalize();
      return mean;
    };

    
    virtual TPrecision getSquaredSimilarity(int index, const VectorXp &p){
      TPrecision a = 1 - p.dot( getPoint(index) );
      return a*a;
    };
      

};


template <typename TPrecision>
class CorrelationGMRAKmeansDataFactory : public IKMKmeansDataFactory<TPrecision>{
  public:
    virtual KmeansData<TPrecision> *createKmeansData(GMRADataObject<TPrecision> *d, std::vector<int> &sub){
      return new CorrelationGMRAKmeansData<TPrecision>(d, sub);
    };
};





















template <typename TPrecision>
class JointCorrelationAndSpatialGMRAKmeansData : public KmeansData<TPrecision>{
  private:
    GMRADataObject<TPrecision> *data;
    int nCor;
    int nSpatial;
    std::vector<int> subset;

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
    
    JointCorrelationAndSpatialGMRAKmeansData(GMRADataObject<TPrecision> *d, int
        nC, int nS, std::vector<int> &sub) : data(d), nCor(nC),
        nSpatial(nS), subset(sub){ 
    };

    virtual VectorXp getPoint(int index){
      return data->getPoint( subset[index] ) ;
    };

    virtual int getNumberOfPoints(){
      return subset.size();
    };
    
    VectorXp getMean( std::vector<int> &pts ){
      VectorXp mean = VectorXp::Zero(data->dimension());
      for(int i=0; i<pts.size(); i++){
        mean += getPoint( pts[i] );
      }
      mean.head( nCor ).normalize();
      mean.tail( nSpatial ) /= pts.size();
      return mean;
    };

    
    virtual TPrecision getSquaredSimilarity(int index, const VectorXp &p){
      VectorXp x = getPoint(index);
      VectorXp x1 = x.head( nCor );
      VectorXp x2 = x.tail( nSpatial );
      VectorXp y1 = p.head( nCor );
      VectorXp y2 = p.tail( nSpatial );

      TPrecision a = 1 - y1.dot( x1 );
      TPrecision b = (x2 - y2).squaredNorm();
      return a*a + b;
    };
      

};


template <typename TPrecision>
class JointCorrelationAndSpatialGMRAKmeansDataFactory : public IKMKmeansDataFactory<TPrecision>{
  private:
    
    int nCor;
    int nSpatial;

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;


    JointCorrelationAndSpatialGMRAKmeansDataFactory(int nC, int nS) : nCor(nC), nSpatial(nS){};

    virtual KmeansData<TPrecision> *createKmeansData(GMRADataObject<TPrecision> *d, std::vector<int> &sub){
      return new JointCorrelationAndSpatialGMRAKmeansData<TPrecision>(d, nCor, nSpatial, sub);
    };


};







/*

template <typename TPrecision>
class JointL2GMRAKmeansData : public KmeansData<TPrecision>{
  private:
    GMRADataObject<TPrecision> *data;
    int n;
    int nSpatial;
    TPrecision lambda;
    std::vector<int> subset;

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
    
    JointCorrelationAndSpatialGMRAKmeansDataFactory(GMRADataObject<TPrecision>
        *d, int nD, int nS, TPrecision l, std::vector<int> &sub) : data(d), nData(nD), nSignal(nS),
    lambda(l), subset(sub){};

    virtual VectorXp getPoint(int index){
      return data->getPoint( subset[index] ) ;
    };

    virtual int getNumberOfPoints(){
      return subset.size();
    };
    
    VectorXp getMean( std::vector<int> &pts ){
      VectorXp mean = VectorXp::Zero(data->dimension());
      for(int i=0; i<pts.size(); i++){
        mean += getPoint( pts[i] );
      }
      mean.head( nCor ).normalize();
      mean.tail( nSpatial ) /= pts.size();
      return mean;
    };

    
    virtual TPrecision getSquaredSimilarity(int index, const VectorXp &p){
      VectorXp x = getPoint(index);
      VectorXp x1 = x.head( nCor );
      VectorXp x2 = x.tail( nSpatial );
      VectorXp y1 = p.head( nCor );
      VectorXp y2 = p.tail( nSpatial );

      TPrecision a = 1 - y1.dot( x1 );
      TPrecision b = (x2 - y2).squaredNorm();
      return a*a + b;
    };
      

};


template <typename TPrecision>
class JointL2GMRAKmeansDataFactory : public IKMKmeansDataFactory<TPrecision>{
  private:
    
    int nData;
    int nSignal;
    TPrecision lambda,

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;


    JointCorrelationAndSpatialGMRAKmeansDataFactory(int nD, int nS, TPrecision l) : nData(nD), nSignal(nS), lambda(l){};

    virtual KmeansData<TPrecision> *createKmeansData(GMRADataObject<TPrecision> *d, std::vector<int> &sub){
      return new JointL2GMRAKmeansData<TPrecision>(d, nCor, nSpatial, sub);
    };

};

*/









#endif
