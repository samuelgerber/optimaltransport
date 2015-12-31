#ifndef GMRAPRUNE_H
#define GMRAPRUNE_H

#include "GMRAVisitor.h"




template <typename TPrecision>
class MinPointsPruneVisitor : public Visitor<TPrecision>{
  private:
    int scale;
    int nPoints;

  public:
    MinPointsPruneVisitor(int s, int n) :scale(s), nPoints(n){};

    virtual void visit(GMRANode<TPrecision> *node){
      if(node->getScale() <= scale){
        std::vector<GMRANode<TPrecision> *> chitlums = node->getChildren();
        for(int i=0; i<chitlums.size(); i++){
          GMRANode<TPrecision> *kid = chitlums[i];
          if( kid->getPoints().size()  < nPoints){
#ifdef VERBOSE
            std::cout << "Pruned node at scale " << kid->getScale() << "with " << kid->getPoints().size() <<" points." << std::endl << std::endl;
#endif
            node->removeChild( kid );
            DeleteVisitor<TPrecision> del;
            GMRATree<TPrecision>::depthFirst(&del, kid);
          }
        }
      }
    };

};






#endif
