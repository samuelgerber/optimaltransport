#ifndef EXPANDNEIGHBORHOODSTRATEGY_H
#define EXPANDNEIGHBORHOODSTRATEGY_H

#include "MultiscaleTransport.h"
#include "TransportLPSolver.h"
#include "LPSolver.h"

#include <list>
#include <map>
#include <set>
#include <vector>


template <typename TPrecision>
class ExpandNeighborhoodStrategy : public NeighborhoodStrategy<TPrecision> {

  
  public:

    
    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;

    typedef typename TransportPlan<TPrecision>::Path Path; 


  private:
    TPrecision expansionFactor;
    TPrecision expansionTolerance;
    int nRefinementIterations;
    int nExpansionAdd;

  public:
    
    
    ExpandNeighborhoodStrategy( TPrecision rFactor, TPrecision eTolerance, int
        nIters, int nAdd=1000000) : expansionFactor(rFactor),
    expansionTolerance(eTolerance), nRefinementIterations(nIters),
    nExpansionAdd(nAdd)  { };

    virtual ~ExpandNeighborhoodStrategy(){
    };




    TransportPlan<TPrecision> *solveNeighborhoodLP(MultiscaleTransportLevel<TPrecision> *source,
        MultiscaleTransportLevel<TPrecision> *target, TransportPlan<TPrecision>
        *sol, TransportPlan<TPrecision> *expandPaths, LPSolver *solver, TPrecision p){

      //Add all arcs within epsilon of current source and target of the active
      //arcs
      int outerAdded = 1;
      int nIter = nRefinementIterations;

      TPrecision prevCost = 0;

      while(outerAdded != 0 && nIter != 0){
        nIter--;
        std::cout << std::endl << "--- Expand strategy ---" << std::endl; 

        //store dual variables for each node
        clock_t t1 = clock();
        TransportLPSolver<TPrecision>::setPotentials( source->getNodes(), target->getNodes(), solver );

        TransportPlan<TPrecision> neighborhoodPaths(source, target);
        getNeighborhodArcs(sol, expandPaths, neighborhoodPaths, p,
            expansionFactor);
        outerAdded = neighborhoodPaths.getNumberOfPaths();
        neighborhoodPaths.pathIteratorBegin();

        clock_t t2 = clock();
        sol->timeRefine += t2 - t1;
        while( !neighborhoodPaths.pathIteratorIsAtEnd() ){ 
          std::cout << "R-factor: " << expansionFactor << std::endl; 

          prevCost = sol->cost;

          clock_t t3 = clock();
          addColumns(solver, sol, neighborhoodPaths, nExpansionAdd );
          clock_t t4 = clock();
          sol->timeRefine += t4 -t3;

          std::cout << "#ncols: " << solver->getNumCols() << std::endl; 
          solver->solveLP();
          sol->cost = solver->getObjectiveValue();
          clock_t t5 = clock();
          sol->timeSolve += t5 - t4;

          if(prevCost - sol->cost <= expansionTolerance * prevCost ){
            break;
          }
        }
      
        TransportLPSolver<TPrecision>::storeLP(sol, solver, p);     
      
      }

      return sol;

    };




  private:



    //Compute all the neighboring arcs of the current optimal solution
    void getNeighborhodArcs(TransportPlan<TPrecision> *sol,
        TransportPlan<TPrecision> *expand, TransportPlan<TPrecision>
        &neighborhoodPaths, double p, TPrecision rFactor){

      typename std::list< std::pair<Path, TPrecision> > toAdd;
      for(expand->pathIteratorBegin(); !expand->pathIteratorIsAtEnd();
          expand->pathIteratorNext()){

        Path &path = expand->pathIteratorCurrent();
        if(path.w > 0){
          TransportNode<TPrecision> *f1 = path.from;
          TransportNode<TPrecision> *t1 = path.to;

          TPrecision rTo = this->getLocalNodeRadius(t1);
          //TPrecision rTo = t1->getLocalNodeRadius();//*rFactor;
          TPrecision rFrom = this->getLocalNodeRadius(f1);
          //TPrecision rFrom = f1->getLocalNodeRadius();//*rFactor;

          TPrecision r = (rTo + rFrom) * rFactor;
          //r = powf(r, p);

          TransportNodeVector f1N = expand->source->getNeighborhood(f1, r);
          //TransportNodeVector f1N = sol->source->getNeighborhood(f1, powf(rFrom*rFactor, p));
          TransportNodeVector t1N = expand->target->getNeighborhood(t1, r);
          //TransportNodeVector t1N = sol->target->getNeighborhood(t1, powf(rTo *
          //      rFactor, p) );

          for(TransportNodeVectorCIterator fIt = f1N.begin(); fIt !=
              f1N.end(); ++fIt){

            for(TransportNodeVectorCIterator tIt = t1N.begin(); tIt !=
                t1N.end(); ++tIt){
              TransportNode<TPrecision> *f2 = *fIt;
              TransportNode<TPrecision> *t2 = *tIt;
              Path p2(f2, t2);

              if( !sol->hasPath(p2) && !neighborhoodPaths.hasPath(p2) ){
                p2.cost = p2.from->getTransportCost(p2.to, p);
                TPrecision rc = p2.cost - p2.from->getPotential() + p2.to->getPotential();
                if(rc <= 0){
                  neighborhoodPaths.addPath( p2 );
                }
              }
            }
          }
        }
      }
      
    };



    int addColumns(LPSolver *solver, TransportPlan<TPrecision> *sol, TransportPlan<TPrecision>
        &neighborhoodPaths, int maxToAdd){ 
      
      int nToAdd = std::min(neighborhoodPaths.getNumberOfPaths(), maxToAdd);
      solver->addCols(nToAdd);

      //add columns to lp
      int nAdded = 0;
      int offset = sol->source->getNodes().size();
      for( ; (!neighborhoodPaths.pathIteratorIsAtEnd() ) && (nAdded < maxToAdd);
          neighborhoodPaths.pathIteratorNext(true) ){
       
        Path &path = neighborhoodPaths.pathIteratorCurrent();
        int index = sol->addPath(path);
      
        //set coefficient
        solver->setCoefficent(index, path.cost);
        solver->setColBounds(index );
        solver->setColStatus(index, LPSolver::LOWER);

        solver->setColConstraints(index, path.from->getID(), offset + path.to->getID() );
        nAdded++;
      }
      return nAdded;

    };

  

};


#endif


