#ifndef MULTISCALESINKHORNTRANSPORT_H
#define MULTISCALESINKHORNTRANSPORT_H

#include "SparseSinkhornTransport.h"
#include "MultiscaleTransport.h"


#include <vector>
#include <queue>
#include <ctime>



template <typename TPrecision>
class SinkhornParameters{


  public:


    SinkhornParameters(){
      lambda=50;
      iterations = 100;
      tolerance = 0.00001;
      threshold = 0;
      maxPathsPerNode = 10;
    };


    double lambda;
    
    double tolerance;
    
    double threshold;

    int iterations;

    int maxPathsPerNode;


};


template <typename TPrecision>
class MultiscaleSinkhornTransport : public MultiscaleTransport<TPrecision> {
  public:


    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;


    typedef typename TransportPlan<TPrecision>::Path Path; 

    typedef SinkhornParameters<TPrecision>  SinkhornParams;


  private:
    SinkhornParams config;


    typedef typename std::pair< TPrecision, Path> MassPath;
    
    typedef std::priority_queue<MassPath> MassPathQueue;

    typedef typename std::map< TransportNode<TPrecision>*, MassPathQueue>
      MassPathMap; 
    typedef typename MassPathMap::iterator MassPathMapIterator;


  public:

    MultiscaleSinkhornTransport(SinkhornParams &params) : config(params) {
    };


  protected:

    TransportPlanSolutions<TPrecision> *solveLP(MultiscaleTransportLevel<TPrecision>
        *source, MultiscaleTransportLevel<TPrecision> *target,
        TransportPlanSolutions<TPrecision> *prevSols, double p, bool lastScale){

      TransportPlanSolutions<TPrecision> *sols = new
        TransportPlanSolutions<TPrecision>(source, target);

      TransportPlan<TPrecision> *prevSol = prevSols->getPrimarySolution();
      TransportPlan<TPrecision> *sol = sols->getPrimarySolution();
      sol->cost = 0;  
      
      //Store solution
      if(prevSol == NULL){
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
      }
      else{

        MassPathMap mpFrom;
        MassPathMap mpTo;
        for(prevSol->pathIteratorBegin(); !prevSol->pathIteratorIsAtEnd();
            prevSol->pathIteratorNext() ){
          Path &path = prevSol->pathIteratorCurrent();

          TPrecision minw  = std::min( path.from->getMass(),
              path.to->getMass());
          if(path.w > minw * config.threshold ){
            MassPath mp(-path.w, path);
            
            mpFrom[path.from].push(mp);
            if(mpFrom[path.from].size () > config.maxPathsPerNode ){
              mpFrom[path.from].pop();
            }
            
            mpTo[path.to].push(mp);
            if(mpTo[path.to].size () > config.maxPathsPerNode ){
              mpTo[path.to].pop();
            }

          }

        }
        /*
        for(MassPathMapIterator it = mpFrom.begin(); it != mpFrom.end(); ++it){
          MassPathVector &v = it->second;
          if(v.size() > config.maxPathsPerNode ){
            std::sort( v.begin(), v.end() );
          }
        }
        for(MassPathMapIterator it = mpTo.begin(); it != mpTo.end(); ++it){
          MassPathVector &v = it->second;
          if(v.size() > config.maxPathsPerNode ){
            std::sort( v.begin(), v.end() );
          }
        }
        */



        for(MassPathMapIterator mit = mpFrom.begin(); mit != mpFrom.end(); ++mit){

          MassPathQueue &q = mit->second;
          while( !q.empty() ) {

            Path path = q.top().second;
            q.pop();
            const TransportNodeVector &fKids = path.from->getChildren();
            const TransportNodeVector &tKids = path.to->getChildren(); 

            for(TransportNodeVectorCIterator fkIt = fKids.begin(); fkIt !=
                fKids.end(); ++fkIt){
              TransportNode<TPrecision> *f2 = *fkIt;
              for(TransportNodeVectorCIterator tkIt = tKids.begin(); tkIt !=
                  tKids.end(); ++tkIt){
                TransportNode<TPrecision> *t2 = *tkIt;

                Path path(f2, t2);
                path.cost = f2->getTransportCost(t2, p) ;
                sol->addPath(path);
              }
            }
          }
        }

        for(MassPathMapIterator mit = mpTo.begin(); mit != mpTo.end(); ++mit){

          MassPathQueue &q = mit->second;
          while( !q.empty() ) {

            Path path = q.top().second;
            q.pop();
            const TransportNodeVector &fKids = path.from->getChildren();
            const TransportNodeVector &tKids = path.to->getChildren(); 

            for(TransportNodeVectorCIterator fkIt = fKids.begin(); fkIt !=
                fKids.end(); ++fkIt){
              TransportNode<TPrecision> *f2 = *fkIt;
              for(TransportNodeVectorCIterator tkIt = tKids.begin(); tkIt !=
                  tKids.end(); ++tkIt){
                TransportNode<TPrecision> *t2 = *tkIt;

                Path path(f2, t2);
                path.cost = f2->getTransportCost(t2, p) ;
                sol->addPath(path);
              }
            }
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


#ifdef VERBOSE
      std::cout << "Filling sparse matrix" << std::endl; 
      std::cout << sol->getNumberOfPaths() << std::endl; 
#endif

      typedef Eigen::Triplet<double> T;
      std::vector<T> KList;
      KList.reserve( sol->getNumberOfPaths() );     
      std::vector<T> UList;
      UList.reserve( sol->getNumberOfPaths() );     
      for(sol->pathIteratorBegin(); !sol->pathIteratorIsAtEnd();
          sol->pathIteratorNext() ){
        Path &path = sol->pathIteratorCurrent();
        double d = path.cost;
        double e = exp(-config.lambda * d);
        int i = path.from->getID();  
        int j = path.to->getID();
        KList.push_back( T(i,j,e) );
        UList.push_back( T(i,j,d*e) );
      }

      Eigen::SparseMatrix<double> K(sNodes.size(), tNodes.size());
      K.setFromTriplets( KList.begin(), KList.end() );
      K.makeCompressed();
      KList.clear();


      Eigen::SparseMatrix<double> U(sNodes.size(), tNodes.size());
      U.setFromTriplets( UList.begin(), UList.end() );
      U.makeCompressed();
      UList.clear();



      SparseSinkhornTransport sinkhorn;

      //Initalize from previous solution
      if(prevSol != NULL ){
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
      sinkhorn.transport(mu, nu, K, U, config.tolerance, config.iterations); 

#ifdef VERBOSE
      std::cout << " Storing sinkhorn transport plan " << std::endl; 
#endif

      Eigen::SparseMatrix<double> map = sinkhorn.getTransportPlan();
      for(sol->pathIteratorBegin(); !sol->pathIteratorIsAtEnd();
          sol->pathIteratorNext() ){
        Path &path = sol->pathIteratorCurrent();
        int i = path.from->getID();  
        int j = path.to->getID();
        path.w = map.coeff(i, j); 
      }

      sol->leftScaling = sinkhorn.getLeftScaling();
      sol->cost = sinkhorn.getDistances()(0); 


      return sols;
    };


};

#endif
