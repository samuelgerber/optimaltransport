#ifndef MULTISCALETRANSPORTLP_H
#define MULTISCALETRANSPORTLP_H


#include "MultiscaleTransport.h"
#include "LPSolver.h"
#include "PropagationStrategy.h" 
#include "NeighborhoodStrategy.h" 
#include "NeighborhoodPropagationStrategy.h" 

#include <list>


template <typename TPrecision>
class MultiscaleTransportLP : public MultiscaleTransport<TPrecision> {

  
  private:


    std::list< NeighborhoodStrategy<TPrecision>*  > neighborhood;
    typedef typename std::list< NeighborhoodStrategy<TPrecision>*  >::iterator NSIterator;
    NeighborhoodStrategy<TPrecision> *lastScaleNeighborhood;

    //use neighborhood strategy until 
    int maxNeighborhoodSize;

    //use propagation strategy one until  maxNeighborhoodSize points afterwards
    //use propagation strategy two
    PropagationStrategy<TPrecision> *propagation1;
    PropagationStrategy<TPrecision> *propagation2;

    LPSolver *solver;




  public:
    
    
    MultiscaleTransportLP(LPSolver *lps) : solver(lps) {
      propagation1 = NULL;
      propagation2 = NULL;
      maxNeighborhoodSize = 10000000;
      lastScaleNeighborhood = NULL;
      propagation1 =  new NeighborhoodPropagationStrategy<TPrecision>(0);
    };

    virtual ~MultiscaleTransportLP(){
      delete solver;
    };



    void setPropagationStrategy1(PropagationStrategy<TPrecision> *ps){
      //if(propagation1 != NULL){
      //  delete propagation1;
      //}
      propagation1 = ps;
    };
    

    void setPropagationStrategy2(PropagationStrategy<TPrecision> *ps){
      //if(propagation2 != NULL){
      //  delete propagation2;
      //}
      propagation2 = ps;
    };


    void setMaxNeighborhoodSize(int n){
      maxNeighborhoodSize = n;
    };


    void addNeighborhodStrategy(NeighborhoodStrategy<TPrecision> *ns){
      neighborhood.push_back(ns);
    };
    

    void setLastScaleNeighborhodStrategy(NeighborhoodStrategy<TPrecision> *ns){
      //if(lastScaleNeighborhood != NULL){
      //  delete lastScaleNeighborhood;
      //}
      lastScaleNeighborhood = ns;
    };



  protected:


    TransportPlanSolutions<TPrecision> *solveLP(MultiscaleTransportLevel<TPrecision>
        *source, MultiscaleTransportLevel<TPrecision> *target,
        TransportPlanSolutions<TPrecision> *prevSol, double p, bool lastScale){



      int nPoints =  source->getNodes().size() + target->getNodes().size();
      
      TransportPlanSolutions<TPrecision> *sols;

      if( nPoints > maxNeighborhoodSize ){
        sols = propagation2->propagate(solver, source, target, prevSol, p,
            lastScale);
      }
      else{
        sols = propagation1->propagate(solver, source, target, prevSol, p,
            lastScale);

        if(prevSol != NULL ){
          for(NSIterator it = neighborhood.begin(); it != neighborhood.end(); ++it){
            TransportPlan<TPrecision> *res = (*it)->solveNeighborhoodLP(source, target, sols->getPrimarySolution(),
               sols->getCombinedPaths(), solver, p); 
	    sols->setPrimarySolution(res); 
          }
          if(lastScale && lastScaleNeighborhood!=NULL){
	    TransportPlan<TPrecision> *res =
		    lastScaleNeighborhood->solveNeighborhoodLP(source, target, sols->getPrimarySolution(),
				    sols->getCombinedPaths(), solver, p);
	    sols->setPrimarySolution(res); 
          }
        }
      } 

      return sols;
    };







};


#endif


