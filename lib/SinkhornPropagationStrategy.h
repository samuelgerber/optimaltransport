#ifndef SINKHORNPROPAGATIONSTRATEGY_H
#define SINKHORNPROPAGATIONSTRATEGY_H

#include "PropagationStrategy.h"

#include "SparseSinkhornTransport.h"
#include "Eigen/Sparse"
#include "Eigen/Dense"

#include <ctime>

#include <vector>


template <typename TPrecision>
class SinkhornPropagationStrategy : public PropagationStrategy<TPrecision>{

  private:
    
    double lambda;
    double tolerance;
    double threshold;
    int iterations;

  public:


    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;

    typedef typename TransportPlan<TPrecision>::Path Path; 



  public:

    SinkhornPropagationStrategy(double l=50, double t=0.000001, double thres=0,
        int iter=100) : lambda(l), tolerance(t), threshold(thres),
                        iterations(iter){ 
    };

    virtual TransportPlanSolutions<TPrecision> *propagate(LPSolver *solver,
        MultiscaleTransportLevel<TPrecision> *source,
        MultiscaleTransportLevel<TPrecision> *target, TransportPlanSolutions<TPrecision>
        *pSol, double p, bool lastScale){

      TransportPlanSolutions<TPrecision> *sols = new TransportPlanSolutions<TPrecision>(source, target); 
      TransportPlan<TPrecision> *sol = sols->getPrimarySolution(); 
      
      TransportPlan<TPrecision> *prevSol = pSol->getPrimarySolution(); 
      
      std::cout << "---- Sinkhorn propagation ------" << std::endl; 
      sol->cost = 0;  
     
      clock_t t1 = clock(); 
      for( prevSol->pathIteratorBegin(); !prevSol->pathIteratorIsAtEnd();
          prevSol->pathIteratorNext() ){

        const Path &path = prevSol->pathIteratorCurrent();
        double mass = std::min( path.from->getMass(), path.to->getMass() );
        if( prevSol->leftScaling.rows() != 0 && path.w < (mass * threshold) ){
          continue;
        }

        const TransportNodeVector &fKids = path.from->getChildren();
        const TransportNodeVector &tKids = path.to->getChildren(); 

        for(TransportNodeVectorCIterator fkIt = fKids.begin(); fkIt !=
            fKids.end(); ++fkIt){
          TransportNode<TPrecision> *f2 = *fkIt;
          for(TransportNodeVectorCIterator tkIt = tKids.begin(); tkIt !=
              tKids.end(); ++tkIt){
            TransportNode<TPrecision> *t2 = *tkIt;

            Path path(f2, t2);
            path.cost =  f2->getTransportCost(t2, p);
            sol->addPath(path);
          }
        }

      }


      TransportNodeVector &sNodes = source->getNodes();
      TransportNodeVector &tNodes = target->getNodes();

      Eigen::VectorXd mu( sNodes.size() );
      Eigen::MatrixXd nu( tNodes.size(), 1);
      for(int i=0; i< sNodes.size(); ++i){ 
        int id = sNodes[i]->getID();
        mu(id) = sNodes[i]->getMass();
      }
      for(int i=0; i< tNodes.size(); ++i){ 
        int id = tNodes[i]->getID();
        nu(id, 0) = tNodes[i]->getMass();
      }


      std::cout << "Filling sparse matrix" << std::endl; 

      typedef Eigen::Triplet<double> T;
      std::vector<T> KList;
      KList.reserve( sol->getNumberOfPaths() );     
      std::vector<T> UList;
      UList.reserve( sol->getNumberOfPaths() );     
      for(sol->pathIteratorBegin(); !sol->pathIteratorIsAtEnd();
          sol->pathIteratorNext() ){
        const Path &path = sol->pathIteratorCurrent();
        double d = path.cost;
        double e = exp(-lambda * d);
        int i = path.from->getID();  
        int j = path.to->getID();
        KList.push_back( T( i, j, e ) );
        UList.push_back( T( i, j, d*e ) );
      }

      Eigen::SparseMatrix<double> K(sNodes.size(), tNodes.size());
      K.setFromTriplets( KList.begin(), KList.end() );
      KList.clear();

      Eigen::SparseMatrix<double> U(sNodes.size(), tNodes.size());
      U.setFromTriplets( UList.begin(), UList.end() );
      UList.clear();


      std::cout << sol->getNumberOfPaths() << std::endl; 

      SparseSinkhornTransport sinkhorn;
      //Initalize from previous solution
      if(prevSol->leftScaling.rows() != 0 ){
        Eigen::MatrixXd leftScaling(sNodes.size(), 1);
        TransportNodeVector &psNodes = prevSol->source->getNodes();
        for(int i = 0; i< psNodes.size(); ++i){

          const TransportNodeVector &kids = psNodes[i]->getChildren();
          int index = psNodes[i]->getID();
          double s = prevSol->leftScaling(index, 0) / kids.size();

          for(TransportNodeVectorCIterator kIt = kids.begin(); kIt !=
              kids.end(); ++kIt){
            TransportNode<TPrecision> *node = *kIt;
            leftScaling(node->getID(), 0) = s;
          }
          sinkhorn.initalizeLeftScaling( leftScaling );
        }

      }
      clock_t t2 = clock(); 
      sinkhorn.transport(mu, nu, K, U, tolerance, iterations); 
      clock_t t3 = clock(); 

      std::cout << " Sinkhorn transport plan " << std::endl; 
      Eigen::SparseMatrix<double> map = sinkhorn.getTransportPlan();
      for(sol->pathIteratorBegin(); !sol->pathIteratorIsAtEnd();
          sol->pathIteratorNext()){
        Path &path = sol->pathIteratorCurrent();
        int i = path.from->getID();  
        int j = path.to->getID();
        path.w = map.coeff(i, j); 
      }

      sol->leftScaling = sinkhorn.getLeftScaling();
      sol->cost = sinkhorn.getDistances()(0); 


      sol->timeSolve += t3 - t2;
      sol->timePropagate += t2 - t1;
      return sols;
    };


};


#endif


