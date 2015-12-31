#ifndef NEIGHBORHOODSTRATEGY_H
#define NEIGHBORHOODSTRATEGY_H

#include "MultiscaleTransport.h"
#include "LPSolver.h"



template <typename TPrecision>
class NeighborhoodStrategy  {

  public:
    
    
    NeighborhoodStrategy() {
    };

    virtual ~NeighborhoodStrategy(){
    };

    virtual TransportPlan<TPrecision> *solveNeighborhoodLP(MultiscaleTransportLevel<TPrecision> *source,
        MultiscaleTransportLevel<TPrecision> *target, TransportPlan<TPrecision>
        *sol, TransportPlan<TPrecision> *nhood, LPSolver *solver, TPrecision p) = 0;


  protected:


    TPrecision getLocalNodeRadius(TransportNode<TPrecision> *node){
      
      TransportNode<TPrecision> *parent = node->getParent();
      if(parent == NULL){
        return node->getLocalNodeRadius();
      }
      return parent->getLocalNodeRadius();

    };








};


#endif


