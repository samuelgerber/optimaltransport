#ifndef GMRANODE_H
#define GMRANODE_H


#include "NodeDistance.h"

#include <vector>




//Node class
template <typename TPrecision>
class GMRANode{
  
  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
    
    typedef typename std::vector< GMRANode<TPrecision> *> NodeVector;
    typedef typename NodeVector::iterator NodeVectorIterator;


  private:

    TPrecision radius;
    TPrecision localRadius;

    class CalculateRadius{
      private:
        NodeDistance<TPrecision> *d;
        GMRANode<TPrecision> *source;
        TPrecision radius;

        void visit(GMRANode<TPrecision> *node){
          if(node->getChildren().size() == 0){
            radius = std::max( d->distance(source, node), radius );
          }

          NodeVector kids = node->getChildren();
          for(unsigned int i=0; i<kids.size(); i++){
            visit(kids[i]);
          }
        };

      public:
        CalculateRadius(GMRANode<TPrecision> *n, NodeDistance<TPrecision> *dist) :
          d(dist), source(n), radius(0){
        };

        TPrecision calculateRadius(){
          radius = 0;
          visit(source);
          return radius;
        };

    };




  public:
    
    GMRANode(){
      radius = -1;
      localRadius = -1;
    };

    virtual ~GMRANode(){};

    //Get a list of children nodes for this node
    virtual NodeVector &getChildren() = 0;

    virtual GMRANode<TPrecision> *getParent() = 0;    
    virtual void setParent(GMRANode<TPrecision> *p) = 0;

    virtual int getScale() = 0;
    virtual void setScale(int s) = 0;

    virtual bool isStop() = 0;
    virtual void setStop(bool s) = 0;

    virtual GMRANode<TPrecision> *findDescendant(const VectorXp &x) = 0;

    //Get the points as indicies into the data matrix X used to construct the
    //tree
    virtual std::vector<int> &getPoints() = 0;
   
    virtual int getIntrinsicDimension() = 0;
    
    //Partition representative and radius 
    virtual VectorXp &getCenter() = 0;
    
    virtual TPrecision getRadius(){
      return radius;
    };

    //Radius with respect to children represntatives
    virtual TPrecision getLocalRadius(){
      return localRadius;
    };


    void removeChild(GMRANode<TPrecision> *kid){
      NodeVector &kids = getChildren();
      for(unsigned int i=0; i<kids.size(); i++){
        if(kids[i] == kid){
          kids.erase( kids.begin()+i );
        }
      }
    };

    virtual void computeRadius(NodeDistance<TPrecision> *dist){
      CalculateRadius rc(this, dist);
      radius = rc.calculateRadius();
    };



    virtual void computeLocalRadius(NodeDistance<TPrecision> *dist){
      NodeVector &kids = getChildren();
      localRadius = 0;
      for(unsigned int i=0; i<kids.size(); i++){
        localRadius = std::max(localRadius, dist->distance(this, kids[i]) );
      }
    };


    virtual void translate(VectorXp &x) = 0;
    virtual void affine(MatrixXp &A) = 0;
   

};





//Node class
template <typename TPrecision>
class GMRANodeBase : public GMRANode<TPrecision>{
  
  public:

    typedef typename GMRANode<TPrecision>::NodeVector NodeVector;
    typedef typename NodeVector::iterator NodeVectorIterator;


  private:

    GMRANode<TPrecision> *parent;
    int scale;
    bool stop;


  public:
    
    GMRANodeBase(){
      parent = NULL;
      stop = false;
      scale = -1;
    };

    virtual ~GMRANodeBase(){
    };

    virtual GMRANode<TPrecision> *getParent(){
      return parent;
    };
    
    virtual void setParent(GMRANode<TPrecision> *p){
      parent = p;
    };

    virtual int getScale(){
      return scale;
    };

    virtual void setScale(int s){
      scale = s;
    };

    virtual bool isStop(){
      return stop;
    };

    virtual void setStop(bool s){
      stop = s;
    };


};



template <typename TPrecision>
class QueryGMRANode : public GMRANodeBase<TPrecision>{
  
  public:
    
    typedef typename std::vector< GMRANode<TPrecision> *> NodeVector;
    typedef typename NodeVector::iterator NodeVectorIterator;


    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;


  private: 
    VectorXp center;

  public:
 
  
    QueryGMRANode(const VectorXp &c) : center(c){
    };

    VectorXp &getCenter(){
      return center;
    };

    NodeVector &getChildren(){
      static NodeVector children;
      return children;
    };

    TPrecision getRadius(){
      return -1;
    };

    std::vector<int> &getPoints(){
      static std::vector<int> empty;
      return empty;
    };

    int getIntrinsicDimension(){
      return center.size();
    };

    virtual GMRANode<TPrecision> *findDescendant(const VectorXp &x){
      return NULL;
    };
    
    virtual void translate(VectorXp &x){
      center += x;
    }
    
    virtual void affine(MatrixXp &A){
    };
}; 






template <typename TPrecision>
class GMRANodeDecorator : public GMRANode<TPrecision>{

  private:
    GMRANode<TPrecision> *node;

  public:
    typedef typename GMRANode<TPrecision>::NodeVector NodeVector;
    typedef typename NodeVector::iterator NodeVectorIterator;
    
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
    
    
    GMRANodeDecorator(GMRANode<TPrecision> *n) : node(n){
    };
    
    virtual ~ GMRANodeDecorator(){
      delete node;
    };


    //Get a list of children nodes for this node
    virtual NodeVector &getChildren() {
      return node->getChildren();
    };

    virtual GMRANode<TPrecision> *getParent(){
      return node->getParent();
    };

    virtual void setParent(GMRANode<TPrecision> *p) {
      node->setParent(p);
    };

    virtual int getScale(){
      return node->getScale();
    };

    virtual void setScale(int s){
      node->setScale(s);
    };

    virtual bool isStop(){
      return node->isStop();
    }

    virtual void setStop(bool s) {
      node->setStop(s);
    };

    virtual GMRANode<TPrecision> * findDescendant(const VectorXp &x){
      return node->findDescendant(x);
    };

    //Get the points as indicies into the data matrix X used to construct the
    //tree
    virtual std::vector<int> &getPoints(){
      return node->getPoints();
    };
   

    virtual int getIntrinsicDimension(){
      return node->getIntrinsicDimension();
    };

    //Partition representative and radius 
    virtual VectorXp &getCenter(){
      return node->getCenter();
    };

    virtual TPrecision getRadius(){
      return node->getRadius();
    };


    virtual TPrecision getLocalRadius(){
      return node->getLocalRadius();
    };


    virtual GMRANode<TPrecision> *getDecoratedNode(){
      return node;
    };

   
    virtual void translate(VectorXp &x){
      node->translate(x);
    };


    virtual void affine(MatrixXp &A){
      node->affine(A);
    }; 
    
    virtual void computeRadius(NodeDistance<TPrecision> *dist){
      node->computeRadius(dist);
    };

    virtual void computeLocalRadius(NodeDistance<TPrecision> *dist){
      node->computeLocalRadius(dist);
    };

};

#endif
