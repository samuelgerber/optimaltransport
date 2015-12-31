#ifndef CAPACITYPROPAGATIONSTRATEGY_H
#define CAPACITYPROPAGATIONSTRATEGY_H

#include "NeighborhoodPropagationStrategy.h"
#include "Random.h"

#include <ctime>


template <typename TPrecision>
class CapacityPropagationStrategy : public NeighborhoodPropagationStrategy<TPrecision>{



  public:


    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;

    typedef typename TransportPlan<TPrecision>::Path Path; 


  private:
    int k;

  public:

    CapacityPropagationStrategy(int kFactor, TPrecision eFactor=0) :
      NeighborhoodPropagationStrategy<TPrecision>(eFactor), k(kFactor){ 
    };

    virtual ~CapacityPropagationStrategy(){};




  protected:

    virtual void computeAlternateSolutions( LPSolver *solver,
                        TransportPlanSolutions<TPrecision> *pSol, double p, bool lastScale ){

      TransportPlan<TPrecision> *sol = pSol->getPrimarySolution();

      clock_t t1 = clock();
      for( sol->pathIteratorBegin(); !sol->pathIteratorIsAtEnd();
          sol->pathIteratorNext() ){
        Path &path = sol->pathIteratorCurrent();            

        
        int nFrom = sol->getNumberOfToPaths( path.from->getID() );
        int nTo = sol->getNumberOfFromPaths( path.to->getID() );
        TPrecision ub =  std::max( path.from->getMass() / ( std::min(k, nFrom) - 0.01 ) , 
                                   path.to->getMass()   / ( std::min(k, nTo)   - 0.01 )  );
        solver->setColBounds(path.index, 0, ub);
      }
      clock_t t2 = clock();

      solver->solveLP();
      clock_t t3 = clock();

      TransportPlan<TPrecision> *res = sol->createCopy();
      res->timePropagate += t2 - t1;
      res->timeSolve += t3 - t2;
      TransportLPSolver<TPrecision>::storeLP(res, solver, p);
      pSol->addAlternativeSolution(res);


    };




};


#endif


