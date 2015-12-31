#ifndef GMRAVISITOR_H
#define GMRAVISITOR_H

#include "GMRANode.h"

#include <vector>





//Visitor for implementing actions on the GMRATree structure
template <typename TPrecision>
class Visitor{
  public:
    //
    virtual void visit(GMRANode<TPrecision> *node) = 0;

};




template <typename TPrecision>
class DeleteVisitor : public Visitor<TPrecision>{
  public:
    //
    virtual void visit(GMRANode<TPrecision> *node){
      delete node;
    };

};


template <typename TPrecision>
class CollectLeavesVisitor : public Visitor<TPrecision>{
  private:
    std::vector<GMRANode<TPrecision>* > leaves;
  public:
    //
    virtual void visit(GMRANode<TPrecision> *node){
      if(node->getChildren().size() == 0){
        leaves.push_back(node);
      }
    };

    void clearLeaves(){
      leaves.clear();
    };

    std::vector<GMRANode<TPrecision> *> &getLeaves(){
      return leaves;
    };

};



#endif
