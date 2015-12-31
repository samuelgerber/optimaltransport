#ifndef IKMNODE_H
#define IKMNODE_H

#include "GMRATree.h"

#include <queue>
#include <map>
#include <iostream>
#include <fstream>



//GMRANode subclass
template <typename TPrecision> 
class IKMNode : public GMRANodeBase<TPrecision>{
  public:
    typedef typename GMRANode<TPrecision>::NodeVector NodeVector;
    typedef typename NodeVector::iterator NodeVectorIterator;
    
    typedef typename GMRADataObject<TPrecision>::MatrixXp MatrixXp;
    typedef typename GMRADataObject<TPrecision>::VectorXp VectorXp;


  private:
    std::vector<int> indices;
    TPrecision kmRadius;
    VectorXp center;
    
    TPrecision mse;
    
    NodeVector children;




  public:    





    IKMNode(){
    };

    IKMNode(VectorXp &mean, std::vector<int> &pts, TPrecision radius, TPrecision meanSE): indices(pts),
         kmRadius(radius), center(mean), mse(meanSE) { 
    };



    virtual ~IKMNode<TPrecision>(){
    };





    int getIntrinsicDimension(){
      return center.size();
    };


    void addChild(GMRANode<TPrecision> *node){
      children.push_back(node);
    };




    int getChildIndex(const VectorXp &x){
      TPrecision dist = std::numeric_limits<TPrecision>::max();
      int index =-1;
      for(int i=0; i<children.size(); i++){
        TPrecision tmp = (children[i]->getCenter() - x).squaredNorm();
        if(tmp < dist){
          dist = tmp;
          index = i;
        }
      }
       
      return index;
    };



    virtual GMRANode<TPrecision> *findDescendant(const VectorXp &x ){
       int index = getChildIndex(x);
       if(index == -1){
         return NULL;
       }
       return children[index];
    };



    virtual NodeVector &getChildren(){
       return children;
    };


    std::vector<int> &getPoints(){
      return indices;
    };
   

    void setPoints(const std::vector<int> &pts){
       indices = pts;
    };


    VectorXp &getCenter(){
      return center;
    };


    void setCenter(const VectorXp &x){
      center = x;
    };


    TPrecision getKMRadius(){
      return kmRadius;
    };
    

    void setKMRadius(TPrecision r){
      kmRadius = r;
    };

    TPrecision getMSE(){
      return mse;
    };


    void setMSE(TPrecision x){
      mse = x;
    };



    virtual void translate(VectorXp &x){
      center += x;
    };


    virtual void affine(MatrixXp &A){
      center = A*center; 
    };

};



#endif
