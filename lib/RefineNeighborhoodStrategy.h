#ifndef REFINENEIGHBORHOODSTRATEGY_H
#define REFINENEIGHBORHOODSTRATEGY_H

#include "NeighborhoodStrategy.h"


template <typename TPrecision>
class RefineNeighborhoodStrategy : public NeighborhoodStrategy<TPrecision> {


  public:


    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;

    typedef typename TransportPlan<TPrecision>::Path Path; 

    typedef LPSolver::Status Status;

  private:
    TPrecision expansionFactor;
    TPrecision expansionTolerance;
    int nRefinementIterations;
    int nExpansionAdd;


  public:

    
    RefineNeighborhoodStrategy( TPrecision rFactor, TPrecision eTolerance, int nIters, int nAdd=1000000) :
      expansionFactor(rFactor), expansionTolerance(eTolerance), nRefinementIterations(nIters), nExpansionAdd(nAdd)  { };

    virtual ~RefineNeighborhoodStrategy(){
    };



    virtual TransportPlan<TPrecision> *solveNeighborhoodLP(MultiscaleTransportLevel<TPrecision> *source,
        MultiscaleTransportLevel<TPrecision> *target, TransportPlan<TPrecision>
        *sol, TransportPlan<TPrecision> *nhood, LPSolver *solver, TPrecision p){
      int nIter = nRefinementIterations;

      TPrecision prevCost = -1;
      while( prevCost != sol->cost && nIter != 0){
        nIter--;

        std::cout << std::endl << "---- Refine strategy ----" << std::endl;

        clock_t t1 = clock();
        prevCost = sol->cost;


        TransportLPSolver<TPrecision>::setPotentials( source->getNodes(), target->getNodes(), solver );

        TransportPlan<TPrecision> *newSol = new
          TransportPlan<TPrecision>(source, target);
        refine(sol, newSol, p, expansionFactor);
        
        //Setup basis
        int nrow = solver->getNumRows();

        std::vector<Status> rowStatus( nrow );

        for( int i=0; i < nrow; i++ ){
          Status status = solver->getRowStatus(i);
          rowStatus[i] = status;
        }

        std::vector< Status > colStatus( newSol->getNumberOfPaths() );

        for(int i=0; i<colStatus.size(); i++){
          colStatus[i] = LPSolver::LOWER;
        }
        for( sol->pathIteratorBegin(); !sol->pathIteratorIsAtEnd();
            sol->pathIteratorNext() ){

          Path &path = sol->pathIteratorCurrent();
          if(path.w > 0){
            int index = newSol->getPathIndex(path);
            colStatus[index] = LPSolver::BASIC;
          }
        }

        newSol->timeRefine = sol->timeRefine; 
        newSol->timeSolve = sol->timeSolve;
        newSol->timePropagate = sol->timePropagate;
        delete sol;
        sol = newSol;

        TransportLPSolver<TPrecision>::createLP(sol, solver); 
        TransportLPSolver<TPrecision>::setupBasis(solver, colStatus, rowStatus);

        clock_t t2 = clock();
        sol->timeRefine += t2 -t1;

        std::cout << "ncols: " << solver->getNumCols() << std::endl;

        double prevCost = sol->cost;
        solver->solveLP();
        clock_t t3 = clock();
        sol->timeSolve += t3 -t2;

        if( (prevCost - sol->cost) <= (expansionTolerance * prevCost) ){
          break;
        }
        
        TransportLPSolver<TPrecision>::storeLP(sol, solver, p);     

      }

      return sol;

    };




  private:

    void refine(TransportPlan<TPrecision> *sol, TransportPlan<TPrecision>
        *newSol, double p, TPrecision rFactor){



      for( sol->pathIteratorBegin(); !sol->pathIteratorIsAtEnd();
          sol->pathIteratorNext() ){

        Path &path = sol->pathIteratorCurrent();

        if( path.w > 0 ){ 
          TransportNode<TPrecision> *f1 = path.from;
          TransportNode<TPrecision> *t1 = path.to;

          if( !newSol->hasPath(path) ){
            newSol->addPath(path);
          }

          TPrecision rTo = this->getLocalNodeRadius(t1);
          //TPrecision rTo = t1->getLocalNodeRadius();// * rFactor;
          TPrecision rFrom = this->getLocalNodeRadius(f1);
          //TPrecision rFrom = f1->getLocalNodeRadius();// * rFactor;

          TPrecision r = (rTo + rFrom) * rFactor;
          //r = powf(r, p);

          TransportNodeVector f1N = sol->source->getNeighborhood(f1, r);
          TransportNodeVector t1N = sol->target->getNeighborhood(t1, r);

          for(TransportNodeVectorCIterator fIt = f1N.begin(); fIt !=
              f1N.end(); ++fIt){

            for(TransportNodeVectorCIterator tIt = t1N.begin(); tIt !=
                t1N.end(); ++tIt){


              TransportNode<TPrecision> *f2 = *fIt;
              TransportNode<TPrecision> *t2 = *tIt;
              Path p2(f2, t2);

              if( !newSol->hasPath(p2) ){
                //PathMapIterator find = allPaths.find(p2);
                //if( find == allPaths.end() )
                p2.cost = f2->getTransportCost(t2, p);
                TPrecision rc = p2.cost - f2->getPotential() + t2->getPotential();
                if( rc <= 0 ){
                  newSol->addPath(p2);
                }
              } 
            }
          }
        }
      }


    };


};


#endif


