#ifndef ITERATEDCAPACITYPROPAGATIONSTRATEGY_H
#define ITERATEDCAPACITYPROPAGATIONSTRATEGY_H

#include "NeighborhoodPropagationStrategy.h"
#include "Random.h"

#include <ctime>


template <typename TPrecision>
class IteratedCapacityPropagationStrategy : public NeighborhoodPropagationStrategy<TPrecision>{



  public:


    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;

    typedef typename TransportPlan<TPrecision>::Path Path; 


  private:
    int nIterations;

  public:

    IteratedCapacityPropagationStrategy(int nIter, TPrecision eFactor=0) :
      NeighborhoodPropagationStrategy<TPrecision>(eFactor), nIterations(nIter){ 
    };

    virtual ~IteratedCapacityPropagationStrategy(){};




  protected:

    virtual void computeAlternateSolutions( LPSolver *solver,
                        TransportPlanSolutions<TPrecision> *pSol, double p, bool lastScale ){

      TransportPlan<TPrecision> *sol = pSol->getPrimarySolution();
    
/*
      for(int i=0; i<nIterations; i++){
      
        TransportPlan<TPrecision> *res = 
            new TransportPlan<TPrecision>( sol->source, sol->target );
        for( sol->pathIteratorBegin(); !sol->pathIteratorIsAtEnd();
          sol->pathIteratorNext() ){
          Path &path = sol->pathIteratorCurrent();            

          if(path.w == 0){
            res->addPath(path);
          }
        }


        TransportLPSolver<TPrecision>::createLP(res, solver);
        solver->solveLP();       
	TransportLPSolver<TPrecision>::storeLP(res, solver, p);
        pSol->addAlternativeSolution(res);
        
        if( !solver->isOptimal() ){
	  std::cout << "suboptimal" << std::endl;
          break;
	}
        sol = res;
      }  

*/
      TransportPlan<TPrecision> *res = sol->createCopy();
      TransportLPSolver<TPrecision>::createLP(res, solver);

      for(int i=0; i<nIterations; i++){
        clock_t t1 = clock();
        for( res->pathIteratorBegin(); !res->pathIteratorIsAtEnd();
          res->pathIteratorNext() ){
          Path &path = res->pathIteratorCurrent();            

          if(path.w > 0){
            static Random<TPrecision> random;
            TPrecision ub =  path.w * ( 0.9 + 0.05*random.Uniform() );
            solver->setColBounds(path.index, 0, ub);
            //path.cost *= 10;
            //solver->setCoefficent(path.index, path.cost);
          }
        }
        clock_t t2 = clock();
        solver->solveLP();
        clock_t t3 = clock();
        

        TransportLPSolver<TPrecision>::storeLP(res, solver, p);
        pSol->addAlternativeSolution(res);
        if( !solver->isOptimal() ){
          std::cout << "suboptimal" << std::endl;
      
          //reset bounds and recompute solution
          /*
           for( res->pathIteratorBegin(); !res->pathIteratorIsAtEnd();
		      res->pathIteratorNext() ){
	      Path &path = res->pathIteratorCurrent();            
	      solver->setColBounds(path.index);
          }
          solver->solveLP();
          */
          break;
        }
        res = res->createCopy();

      }
      

      //pSol->addAlternativeSolution(res);


    };




};


#endif


