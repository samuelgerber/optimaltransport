#ifndef GMRATREE_H
#define GMRATREE_H

#include "GMRANode.h"
#include "GMRAVisitor.h"
#include "GMRADecorator.h"
#include "GMRADataObject.h"
#include "NodeDistance.h"

#include <Eigen/Dense>

#include <list>
#include <vector>
#include <iostream>







//Tree data structure
template <typename TPrecision>
class GMRATree{


private:

  GMRANode<TPrecision> *root;


protected:

  GMRADataObject<TPrecision> *data;

public:
    typedef typename GMRADataObject<TPrecision>::MatrixXp MatrixXp;
    typedef typename GMRADataObject<TPrecision>::VectorXp VectorXp;

    
    GMRATree(GMRADataObject<TPrecision> *D):data(D){
      root = NULL;
    };

   
    virtual ~GMRATree(){
      DeleteVisitor<TPrecision> del;
      depthFirstVisitor(&del);
    };


    //Add points to the GMRA Tree from data object
    virtual void addPoints(std::vector<int> &points) = 0;
   
    //Get Leaf node for data point x
    virtual std::vector<GMRANode<TPrecision> *> getLeafPath(const VectorXp &x) = 0;

    VectorXp getPoint(int i){
      return data->getPoint(i);
    };

    //Get root node of the tree
    GMRANode<TPrecision> *getRoot(){ 
      return root;
    };

    GMRADataObject<TPrecision> *getDataObject(){
      return data;
    };


    void setupParents(){
      
      //Set up parent pointers
      class SetParent : public Visitor<TPrecision>{
        public:
          //
          virtual void visit(GMRANode<TPrecision> *node){
            int scale = 1;
            GMRANode<TPrecision> *p = node->getParent();
            if(p != NULL){
              scale = node->getScale()+1;
            }
            else{
              node->setScale(0);
            }
            std::vector< GMRANode<TPrecision> * > &children = node->getChildren();
            for(unsigned int i=0; i<children.size(); i++){
              children[i]->setParent( node );
              children[i]->setScale( scale );
            }
          };
      
      };

      SetParent parenter;
      //Breadth first required  for setting scale correctly
      breadthFirstVisitor(&parenter);

    };


    void setRoot(GMRANode<TPrecision> *r){
      root= r;
    };



    //Pass each node in depth first order to the Visitor v
    void depthFirstVisitor(Visitor<TPrecision> *v){
      GMRATree<TPrecision>::depthFirst( v, this->getRoot() );
    };

    static void depthFirst(Visitor<TPrecision> *v, GMRANode<TPrecision> *node){
      std::vector<GMRANode<TPrecision> *> children = node->getChildren();
      for(typename std::vector< GMRANode<TPrecision>* >::iterator it = children.begin(); it !=
          children.end(); ++it){ 
        GMRATree<TPrecision>::depthFirst( v, *it );
      }
      v->visit( node );
    };





    //Pass each node in breadth first order to the Visitor v
    void breadthFirstVisitor(Visitor<TPrecision> *v){
      std::list<GMRANode<TPrecision> *> nodes;
      nodes.push_back( getRoot() );
      while(!nodes.empty()){
        GMRANode<TPrecision> *node = nodes.front();
        nodes.pop_front();

        v->visit(node);
        std::vector<GMRANode<TPrecision> *> children = node->getChildren();
        for(typename std::vector<GMRANode<TPrecision>*>::iterator it = children.begin(); it !=
           children.end(); ++it){ 
          nodes.push_back(*it);
        }
      }
    };




    void decorate(Decorator<TPrecision> &decorator){
     this->setRoot( decorate(this->getRoot(), NULL, decorator) );
    };




    void computeRadii(NodeDistance<TPrecision> *dist){
      class RadiusVisitor : public Visitor<TPrecision>{
        private:
          NodeDistance<TPrecision> *dist;
        public:

          RadiusVisitor(NodeDistance<TPrecision> *d):dist(d){
          };

          void visit(GMRANode<TPrecision> *node){
              node->computeRadius(dist);  
          };
      };

      RadiusVisitor vis(dist);
      this->breadthFirstVisitor(&vis);

    };




    void computeLocalRadii(NodeDistance<TPrecision> *dist){
      class RadiusVisitor : public Visitor<TPrecision>{
        private:
          NodeDistance<TPrecision> *dist;
        public:

          RadiusVisitor(NodeDistance<TPrecision> *d):dist(d){
          };

          void visit(GMRANode<TPrecision> *node){
              node->computeLocalRadius(dist);  
          };
      };

      RadiusVisitor vis(dist);
      this->breadthFirstVisitor(&vis);

    };






  private:


    GMRANode<TPrecision> *decorate(GMRANode<TPrecision> *node,
        GMRANode<TPrecision> *parent, Decorator<TPrecision> &decorator){
      
      GMRANode<TPrecision> *dnode = decorator.decorate(node, parent);
      dnode->setParent(parent);

      std::vector< GMRANode<TPrecision>* > &children = node->getChildren();
      for( int i=0; i < children.size(); i++ ){ 
        children[i] =  decorate( children[i], dnode, decorator );
      }

      return dnode;

    };


};




#endif
