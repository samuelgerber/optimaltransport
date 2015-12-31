#ifndef NEIGHBORHOODPROPAGATIONSTRATEGY_H
#define NEIGHBORHOODPROPAGATIONSTRATEGY_H

#include "PropagationStrategy.h"

#include <ctime>
#include "TransportLPSolver.h"

template <typename TPrecision>
class NeighborhoodPropagationStrategy : public PropagationStrategy<TPrecision>{



  public:


    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;

    typedef typename TransportPlan<TPrecision>::Path Path; 


  private:
    TPrecision expansionFactor;

  public:

    NeighborhoodPropagationStrategy(TPrecision eFactor=0) :  expansionFactor(eFactor) {
    };

    virtual ~NeighborhoodPropagationStrategy(){};

    virtual TransportPlanSolutions<TPrecision> *propagate( LPSolver *solver,
        MultiscaleTransportLevel<TPrecision> *source,
        MultiscaleTransportLevel<TPrecision> *target,
        TransportPlanSolutions<TPrecision> *pSol, double p, bool lastScale){
      

      if(pSol == NULL){
        TransportPlanSolutions<TPrecision> *sols = 
              new TransportPlanSolutions<TPrecision>(source, target);
        TransportPlan<TPrecision> *sol = sols->getPrimarySolution();
        //Add all variables
        const TransportNodeVector &sourceNodes = source->getNodes();
        for(TransportNodeVectorCIterator sIt = sourceNodes.begin(); sIt !=
            sourceNodes.end(); ++sIt){

          const TransportNodeVector &targetNodes = target->getNodes();
          for(TransportNodeVectorCIterator tIt = targetNodes.begin(); tIt !=
              targetNodes.end(); ++tIt){
            Path path(*sIt, *tIt);
            path.cost = (*sIt)->getTransportCost( *tIt, p );
            sol->addPath(path);
          }
        }
        
        TransportLPSolver<TPrecision>::createLP(sol, solver); 
        solver->solveLP(); 
        TransportLPSolver<TPrecision>::storeLP(sol, solver, p);     
        
        return sols;
      }
      
     
      
      this->computeAlternateSolutions(solver, pSol, p, lastScale);


      TransportPlanSolutions<TPrecision> *sols = 
           new TransportPlanSolutions<TPrecision>(source, target);
      TransportPlan<TPrecision> *sol = sols->getPrimarySolution();
     
      clock_t t1 = clock(); 
      this->ballNeighborhood(sol, pSol->getPrimarySolution(), 0, p, solver, false);
      TransportLPSolver<TPrecision>::createLP(sol, solver); 
      clock_t t2 = clock(); 
      solver->solveLP(); 
      clock_t t3 = clock(); 
      this->ballNeighborhood(sol, pSol->getCombinedPaths(), expansionFactor, p, solver, true);
      clock_t t4 = clock(); 
      solver->solveLP(); 
      clock_t t5 = clock(); 

      
      sol->timeSolve += t3 - t2 + t5 - t4;
      sol->timePropagate += t2 - t1 + t4 - t3;
      TransportLPSolver<TPrecision>::storeLP(sol, solver, p);     

      return sols;
    
    };







  protected:

    virtual void computeAlternateSolutions( LPSolver *solver,
      TransportPlanSolutions<TPrecision> *pSol, double p, bool lastScale ){
    };






  private:


    void addAllCombinations(const TransportNodeVector &fKids, const
        TransportNodeVector &tKids, TransportPlan<TPrecision> *sol, double p, std::vector<Path> &toAdd){

      for(TransportNodeVectorCIterator fkIt = fKids.begin(); fkIt != fKids.end();
          ++fkIt){
          
        TransportNode<TPrecision> *f2 = *fkIt;
        
        for(TransportNodeVectorCIterator tkIt = tKids.begin(); tkIt !=
            tKids.end(); ++tkIt){
        
          TransportNode<TPrecision> *t2 = *tkIt;

          Path path(f2, t2);

          if( !sol->hasPath(path) ){
            path.cost = f2->getTransportCost(t2, p);
            toAdd.push_back(path);
          }
        
        }
      }

    };



        
    void ballNeighborhood( TransportPlan<TPrecision> *sol,
        TransportPlan<TPrecision>  *prevSol, TPrecision rFactor, double p, LPSolver *solver, bool add){

      //TPrecision rFrom = prevSol->source->getMaximalRadius();
      //TPrecision rTo = prevSol->target->getMaximalRadius();
      //TPrecision r = (rTo + rFrom) * rFactor;

      std::vector<Path> toAdd;
      double sumw = 0;
      for(prevSol->pathIteratorBegin(); !prevSol->pathIteratorIsAtEnd();
          prevSol->pathIteratorNext()){
     
        Path &path = prevSol->pathIteratorCurrent();
        //No mass moved at previous solution 
        if(path.w == 0){
          continue;
        }
        sumw += path.w;

        TransportNode<TPrecision> *f1  = path.from;
        TransportNode<TPrecision> *t1  = path.to;

        if(rFactor > 0){
          TPrecision rFrom = f1->getLocalNodeRadius() * rFactor;
          TPrecision rTo   = t1->getLocalNodeRadius() * rFactor;
          //TPrecision r = (rTo + rFrom) * rFactor;

          const TransportNodeVector &f1n = prevSol->source->getNeighborhood(f1, rFrom);
          const TransportNodeVector &t1n = prevSol->target->getNeighborhood(t1, rTo);


          for(TransportNodeVectorCIterator fIt = f1n.begin(); fIt != f1n.end();
              ++fIt){
            for(TransportNodeVectorCIterator tIt = t1n.begin(); tIt != t1n.end();
                ++tIt){

              const TransportNodeVector &fKids = (*fIt)->getChildren();
              const TransportNodeVector &tKids = (*tIt)->getChildren();

              addAllCombinations(fKids, tKids, sol, p, toAdd);

            } 
          }

        }
        else{
          const TransportNodeVector &fKids = f1->getChildren();
          const TransportNodeVector &tKids = t1->getChildren();

          addAllCombinations(fKids, tKids, sol, p, toAdd);
        }

      }


      int offset = sol->source->getNodes().size();
      if(add){
        solver->addCols(toAdd.size() );
      }
      for(typename std::vector<Path>::iterator it = toAdd.begin(); it != toAdd.end(); ++it){
        Path &path = *it; 
        int index = sol->addPath(path);

        if(add){
          solver->setCoefficent(index, path.cost);
          solver->setColBounds(index );
          solver->setColStatus(index, LPSolver::LOWER);

          solver->setColConstraints(index, path.from->getID(), offset + path.to->getID() );
        } 
      }
      
      std::cout << "Neighborhood weight sum: " << sumw << std::endl;

    
    };



};


#endif


