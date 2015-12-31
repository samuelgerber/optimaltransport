#ifndef IPCANODE_H
#define IPCANODE_H

#include "GMRATree.h"

#include <queue>
#include <map>
#include <iostream>
#include <fstream>



//GMRANode subclass
template <typename TPrecision> 
class IPCANode : public GMRANodeBase<TPrecision>{
  public:
    typedef typename GMRANode<TPrecision>::NodeVector NodeVector;
    typedef typename NodeVector::iterator NodeVectorIterator;
    
    typedef typename GMRADataObject<TPrecision>::MatrixXp MatrixXp;
    typedef typename GMRADataObject<TPrecision>::VectorXp VectorXp;


  private:
    std::vector<int> indices;
    MatrixXp phi;
    VectorXp sigma;
    TPrecision l2Radius;
    VectorXp center;
    
    MatrixXp dir;
    VectorXp a; 
    
    TPrecision totalVar;
    
    NodeVector children;
    std::map<int, int> childmap;



  public:    





    IPCANode(){
    };

    IPCANode(VectorXp &mean, std::vector<int> &pts, MatrixXp &phiIn, VectorXp &sigmaIn,
        TPrecision radius, MatrixXp &splitDir, VectorXp &split, TPrecision tV): indices(pts),
         phi(phiIn), sigma(sigmaIn), l2Radius(radius), center(mean), dir(splitDir),
         a(split), totalVar(tV){ 
    };



    virtual ~IPCANode<TPrecision>(){
    };





    int getIntrinsicDimension(){
      return phi.cols();
    };


    void addChild(GMRANode<TPrecision> *node, int childIndex){
      children.push_back(node);
      childmap[childIndex] = children.size()-1;
    };


    int getMaxKids(){
      return pow(2.0, dir.cols() );
    };


    GMRANode<TPrecision> *getChild(int childIndex){
      std::map<int, int>::iterator it = childmap.find(childIndex);
      if(it == childmap.end()){
        return NULL;
      }
      return children[it->second];
    };




    int getChildIndex(const VectorXp &x){
      VectorXp s =  dir.transpose() *x;

      int childIndex = 0;
      int factor=1;
      for(int i=0; i<s.size(); i++){
        if( s(i) > a(i) ){
          childIndex += factor;
        }
        factor *= 2;
      }
      return childIndex;

    };



    virtual GMRANode<TPrecision> *findDescendant(const VectorXp &x ){
       int index = getChildIndex(x);
       return getChild(index);
    };



    virtual NodeVector &getChildren(){
       return children;
    };


    std::vector<int> &getPoints(){
      return indices;
    };

    VectorXp &getCenter(){
      return center;
    };


    TPrecision getL2Radius(){
      return l2Radius;
    };

    TPrecision getTotalVariance(){
      return totalVar;
    };

    VectorXp &getSigma(){
      return sigma;
    };

    MatrixXp &getPhi(){
      return phi;
    };



    virtual void translate(VectorXp &x){
      center += x;
    };


    virtual void affine(MatrixXp &A){
      center = A*center; 
      phi = A * phi;
    };

};



#endif
