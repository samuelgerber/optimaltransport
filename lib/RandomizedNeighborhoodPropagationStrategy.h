#ifndef RANDOMIZEDNEIGHBORHOODPROPAGATIONSTRATEGY_H
#define RANDOMIZEDNEIGHBORHOODPROPAGATIONSTRATEGY_H

#include "NeighborhoodPropagationStrategy.h"
#include "Random.h"

#include <ctime>


template <typename TPrecision>
class RandomizedNeighborhoodPropagationStrategy : public NeighborhoodPropagationStrategy<TPrecision>{



  public:


    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;

    typedef typename TransportPlan<TPrecision>::Path Path; 


  private:
    int nRandom;

  public:

    RandomizedNeighborhoodPropagationStrategy(int n, TPrecision eFactor=0) :
      NeighborhoodPropagationStrategy<TPrecision>(eFactor), nRandom(n){ 
    };

    virtual ~RandomizedNeighborhoodPropagationStrategy(){};




  protected:

    virtual void computeAlternateSolutions(LPSolver *solver,
                      TransportPlanSolutions<TPrecision> *pSol, double p, bool lastScale ){

      TransportPlan<TPrecision> *sol = pSol->getPrimarySolution();
      for(int i=0; i<nRandom; i++){
         
        for( sol->pathIteratorBegin(); !sol->pathIteratorIsAtEnd();
             sol->pathIteratorNext() ){
          Path &path = sol->pathIteratorCurrent();            
          
          TransportNode<TPrecision> *from = path.from;
          TransportNode<TPrecision> *to = path.to;
          TPrecision r = from->getNodeRadius() + to->getNodeRadius();
          TPrecision dist = pow(path.cost, 1.0/p);
          TPrecision delta = pow(dist+r, p) - pow(dist-r, p);

          //TPrecision change = (random.Uniform()-0.5) * delta;
          //TPrecision change = (random.Uniform()-0.5) * delta/2.0;
          static Random<TPrecision> random;
          TPrecision change = random.Normal() * delta/5.0;
          solver->setCoefficent(path.index, path.cost + change );
        }

        solver->solveLP();

        TransportPlan<TPrecision> *res = sol->createCopy();
        TransportLPSolver<TPrecision>::storeLP(res, solver, p);
        pSol->addAlternativeSolution(res);

      }

    };




};


#endif


