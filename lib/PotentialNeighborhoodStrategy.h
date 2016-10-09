#ifndef POTENTIALNEIGHBORHOODSTRATEGY_H
#define POTENTIALNEIGHBORHOODSTRATEGY_H

#include "NeighborhoodStrategy.h"

#include <list>

template <typename TPrecision>
class ReducedCostPath{
  public:
    typedef typename TransportPlan<TPrecision>::Path Path; 

    Path path;
    TPrecision reducedCost;

    ReducedCostPath(Path &p, TPrecision rc) : path(p),
      reducedCost(rc){
    };

  
    bool operator == (const ReducedCostPath& other) const{
      return this->reducedCost == other.reducedCost;
    };

    bool operator < (const ReducedCostPath& other) const{
      return this->reducedCost < other.reducedCost;
    };

    bool operator > (const ReducedCostPath& other) const{
      return this->reducedCost > other.reducedCost;
    };
  

};



template <typename TPrecision>
class PotentialNeighborhoodStrategy : public NeighborhoodStrategy<TPrecision> {


  public:

    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;

    typedef typename TransportPlan<TPrecision>::Path Path; 

    typedef LPSolver::Status Status;


  private:
    
    typedef typename std::list< ReducedCostPath<TPrecision> > RCList;
    typedef typename RCList::iterator RCListIterator;
    
    TPrecision reducedCostThresholdFactor;
    TPrecision expansionTolerance;
    bool sortReducedCost;
    bool expandPotential;
    int nRefinementIterations;
    int nExpansionAdd;



  public:


    PotentialNeighborhoodStrategy(TPrecision threshold, TPrecision eTolerance,
        bool sort=false, bool expand=false,  int nIters=1, int nAdd=1000000) :
      reducedCostThresholdFactor(threshold), expansionTolerance(eTolerance),
      sortReducedCost(sort), expandPotential(expand),
      nRefinementIterations(nIters), nExpansionAdd(nAdd){ 
    };


    virtual ~PotentialNeighborhoodStrategy(){
    };




    TransportPlan<TPrecision> *solveNeighborhoodLP(MultiscaleTransportLevel<TPrecision> *source,
        MultiscaleTransportLevel<TPrecision> *target, TransportPlan<TPrecision>
        *sol, TransportPlan<TPrecision> *nhood, LPSolver *solver, TPrecision p){ 





      for(int nIter =nRefinementIterations; nIter !=0; --nIter){
        std::cout << std::endl << "---- Potential strategy ----" << std::endl;

        clock_t t1 = clock(); 
        
        MultiscaleTransportLevel<TPrecision> *rootS = source->getRootLevel();
        MultiscaleTransportLevel<TPrecision> *rootT = target->getRootLevel();

        int tScale = target->getScale();
        int sScale = source->getScale();
        TransportNodeVector &targetRootNodes = rootT->getNodes();
        TransportNodeVector &sourceRootNodes = rootS->getNodes();

        //setup potential bounds in each node of the multiscale hierarchy
        TransportNodeVector &targetNodes = target->getNodes();
        TransportNodeVector &sourceNodes = source->getNodes();



        //store dual variables for each node
        TransportLPSolver<TPrecision>::setPotentials(sourceNodes, targetNodes, solver);

        //propagate bounds to top of target transport hierarchy
        for(TransportNodeVectorIterator tIt = targetRootNodes.begin(); tIt
            != targetRootNodes.end(); ++tIt){ 
          potentialBounds( *tIt, rootT->getScale(), tScale );
        }

        //propagate bounds to top of source transport hierarchy
        for(TransportNodeVectorIterator sIt = sourceRootNodes.begin(); sIt
            != sourceRootNodes.end(); ++sIt){
          potentialBounds( *sIt, rootS->getScale(), sScale );
        }




        //Find edges that are possibly included in optimal transport plan
        int nrow = solver->getNumRows();

        RCList rcArcs;
        int nComparisons = 0;

        for(TransportNodeVectorIterator sIt = sourceNodes.begin(); sIt != sourceNodes.end(); ++sIt){
          TransportNode<TPrecision> *nFrom = *sIt;

          std::list< TransportNode<TPrecision>* > queue;
          queue.insert( queue.end(), targetRootNodes.begin(), targetRootNodes.end() );
          std::list<int> scales;
          scales.insert( scales.end(), targetRootNodes.size(), rootT->getScale() );

          while( !queue.empty() ){
            nComparisons++;

            TransportNode<TPrecision> *nTo = queue.front();
            queue.pop_front();
            int s = scales.front();
            scales.pop_front();

            TPrecision cost = -1;
            cost = nFrom->getTransportCost(nTo, 1);


            TPrecision d = cost;
            if(s < tScale){
              d = d - nTo->getNodeRadius();
            }
            if(d > 0 ){
              d = pow(d, p);
            }
            TPrecision delta = nFrom->getPiMax() - nTo->getPiMin();
            //TPrecision delta = nFrom->getPiMax() + nTo->getPiMax();

            TPrecision rc = d-delta;

            if( rc <= reducedCostThresholdFactor * pow(cost,p) ){
              if(s == tScale){              
                Path path(nFrom, nTo);
                path.cost = pow(cost, p);
                rcArcs.push_back( ReducedCostPath<TPrecision>(path, rc) );
              }
              else{
                const TransportNodeVector &kids = nTo->getChildren();
                queue.insert( queue.end(), kids.begin(), kids.end() );
                scales.insert( scales.end(), kids.size(), s+1 );
              }

            }
          }
          if(rcArcs.size() > nExpansionAdd){
            break;
          }
        }


        std::cout << "#Comparisons: " << nComparisons << std::endl;

        //int cutoff = rcArcs.size(); //nExpansionAdd;

        //if(sortReducedCost){ 
        //  rcArcs.sort();
        //}

        TransportPlan<TPrecision> *newSol = new
          TransportPlan<TPrecision>(source, target);

        //int nAdded = 0;
        for(RCListIterator it = rcArcs.begin(); it != rcArcs.end(); ++it){
          ReducedCostPath<TPrecision> &rcPath = *it;
          newSol->addPath( rcPath.path );
        }

        //copy current solution from previous solution
        //col stats
        std::vector<Status> colStatus( newSol->getNumberOfPaths(), LPSolver::LOWER );

        for( sol->pathIteratorBegin(); !sol->pathIteratorIsAtEnd();
             sol->pathIteratorNext() ){


          Path &path = sol->pathIteratorCurrent();
          Status status = solver->getColStatus(path.index);
          if(status == LPSolver::BASIC){

            int index = newSol->getPathIndex(path);
            if( index == -1 ){
              newSol->addPath( path );
              colStatus.push_back( LPSolver::BASIC );
            }
            else{ 
              colStatus[ index ] = LPSolver::BASIC;
            }
          }
          else if( expandPotential ){
            int index =  newSol->addPath( path );
            if( index == colStatus.size() ){
              colStatus.push_back( LPSolver::LOWER );
            }
          }



        }


        //row stats
        std::vector<Status>  rowStatus( nrow );
        for( int i=0; i < nrow; i++ ){
          rowStatus[ i ] = solver->getRowStatus(i);
        }


        //setup new problem
        newSol->timeRefine = sol->timeRefine; 
        newSol->timePropagate = sol->timePropagate; 
        newSol->timeSolve = sol->timeSolve;

        delete sol;
        sol = newSol;

        std::cout << "Potential strategy #Paths: " << newSol->getNumberOfPaths() <<
          std::endl;
        std::cout << "Potential strategy #Reduced Cost Paths: " << rcArcs.size() <<
          std::endl;


        TransportLPSolver<TPrecision>::createLP(sol, solver); 
        TransportLPSolver<TPrecision>::setupBasis(solver, colStatus, rowStatus);

        std::cout << "ncols: " << solver->getNumCols() << std::endl;

        double prevCost = sol->cost;

        clock_t t2 = clock();
        sol->timeRefine += t2 - t1;

        solver->solveLP();
        sol->cost = solver->getObjectiveValue();

        clock_t t3 = clock();
        sol->timeSolve += t3 - t2;

        TransportLPSolver<TPrecision>::storeLP(sol, solver, p);

        if( solver->getIterationCount() == 0 ){//|| sol->cost  - prevCost < 10e-15)
          break;
        }

        if( (prevCost - sol->cost) <= (expansionTolerance * prevCost) ){
          break;
        }

      }
      return sol;
    };




  private:


      void potentialBounds(TransportNode<TPrecision> *node, int currentScale, int stopScale){
        if(currentScale == stopScale){
          return;
        }
        node->resetPi();

        const TransportNodeVector &kids = node->getChildren();
        for(int i=0; i< kids.size(); i++){
          potentialBounds( kids[i], currentScale+1, stopScale );
        }

        for(int i=0; i< kids.size(); i++){
          node->setPiMax( kids[i]->getPiMax() );
          node->setPiMin( kids[i]->getPiMin() );
        }

      };

    };


#endif


