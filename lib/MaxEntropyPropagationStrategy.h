#ifndef MAXENTROPYPROPAGATIONSTRATEGY_H
#define MAXENTROPYPROPAGATIONSTRATEGY_H

#include "PropagationStrategy.h"

#include <ctime>


template <typename TPrecision>
class MaxEntropyPropagationStrategy : public PropagationStrategy<TPrecision> {


  public:


    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;

    typedef typename TransportPlan<TPrecision>::Path Path; 



  public:


    virtual TransportPlanSolutions<TPrecision> *propagate(LPSolver *solver,
        MultiscaleTransportLevel<TPrecision> *source,
        MultiscaleTransportLevel<TPrecision> *target, TransportPlanSolutions<TPrecision>
        *pSol, double p, bool lastScale){
    
      TransportPlanSolutions<TPrecision> *sols = new TransportPlanSolutions<TPrecision>(source, target); 
      TransportPlan<TPrecision> *sol = sols->getPrimarySolution(); 
      
      TransportPlan<TPrecision> *prevSol = pSol->getPrimarySolution(); 
      

      sol->cost = 0;  
      //Store solution
      for( prevSol->pathIteratorBegin(); !prevSol->pathIteratorIsAtEnd();
          prevSol->pathIteratorNext() ){
        const Path &path = prevSol->pathIteratorCurrent();


        const TransportNodeVector &fKids = path.from->getChildren();
        const TransportNodeVector &tKids = path.to->getChildren(); 

        TPrecision sumw =  0;
        for(TransportNodeVectorCIterator fkIt = fKids.begin(); fkIt !=
            fKids.end(); ++fkIt){
          TransportNode<TPrecision> *f2 = *fkIt;
          for(TransportNodeVectorCIterator tkIt = tKids.begin(); tkIt !=
              tKids.end(); ++tkIt){
            TransportNode<TPrecision> *t2 = *tkIt;
            sumw += f2->getMass() * t2->getMass();
          }
        }
        for(TransportNodeVectorCIterator fkIt = fKids.begin(); fkIt !=
            fKids.end(); ++fkIt){
          TransportNode<TPrecision> *f2 = *fkIt;
          for(TransportNodeVectorCIterator tkIt = tKids.begin(); tkIt !=
              tKids.end(); ++tkIt){
            TransportNode<TPrecision> *t2 = *tkIt;

            Path path(f2, t2);
            path.cost = f2->getTransportCost(t2, p) ;
            path.w = f2->getMass() * t2->getMass() * path.w / sumw;
            sol->addPath(path);
            sol->cost += path.w * path.cost;
          }
        }

      } 

      //Potential of from nodes
      for(TransportNodeVectorCIterator it = prevSol->source->getNodes().begin(); it !=
          prevSol->source->getNodes().end(); ++it){
        TransportNode<TPrecision> *node = *it;
        const TransportNodeVector &kids = node->getChildren();
        TPrecision pi = node->getPotential();
        for(int i=0; i < kids.size(); i++){
          kids[i]->setPotential(pi);
        }
      }


      //Potential of to nodes
      for(TransportNodeVectorCIterator it = prevSol->target->getNodes().begin(); it !=
          prevSol->target->getNodes().end(); ++it){
        TransportNode<TPrecision> *node = *it;
        const TransportNodeVector &kids = node->getChildren();
        TPrecision pi = node->getPotential();
        for(int i=0; i < kids.size(); i++){
          kids[i]->setPotential(pi);
        }
      }


      return sols;

    };


};


#endif


