#ifndef NULL
#define NULL 0
#endif

#define R_NO_REMAP


#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdio.h>


#include "IPCATree.h"
#include "IPCAGWT.h"
#include "GMRATree.h"
#include "GMRANeighborhood.h"
#include "RelativeGMRANeighborhood.h"

#include "LPSolver.h"
#include "MultiscaleTransportLP.h"
#include "ExpandNeighborhoodStrategy.h"
#include "RefineNeighborhoodStrategy.h"
#include "PotentialNeighborhoodStrategy.h"

#include "NeighborhoodPropagationStrategy.h"
#include "CapacityPropagationStrategy.h"
#include "IteratedCapacityPropagationStrategy.h"
#include "RandomizedNeighborhoodPropagationStrategy.h"
#include "MaxEntropyPropagationStrategy.h"
#include "SinkhornPropagationStrategy.h"
#include "GMRAMultiscaleTransport.h"


#include "TransportLP.h"

#include "mop_config.h"

#ifdef MOP_USE_GLPK
#include "GLPKSolver.h"
#endif

#ifdef MOP_USE_MOSEK
#include "MOSEKSolver.h"
#endif

#ifdef MOP_USE_CPLEX
#include "CPLEXSolver.h"
#include "CPLEXNetworkSolver.h"
#endif

#include "NodeDistance.h"
#include "WassersteinNodeDistance.h"
#include "EigenL1Metric.h"
#include "EigenMetric.h"
#include "EigenEuclideanMetric.h"
#include "EigenSquaredEuclideanMetric.h"





extern "C" {
  int trpCounter;
  int propCounter;
  
  /* called from .onLoad */
  SEXP mop_load(void){
    trpCounter = 0;
    propCounter = 0;
    
    return R_NilValue;
  };





  enum  Optimizer {
    GLPK_PRIMAL_SIMPLEX = 0, 
    GLPK_DUAL_SIMPLEX   = 1, 
    GLPK_INTPOINT       = 2,
    GLPK                = 10,
    MSK_SIMPLEX_FREE    = 11,
    MSK_NETWORK_SIMPLEX = 12,
    MSK_FREE            = 13,
    MSK_INTPOINT        = 14,
    MSK_CONCURRENT      = 15,
    MSK                 = 20,
    CPLEX_AUTO          = 21,      
    CPLEX_BARRIER       = 22,      
    CPLEX_DUAL          = 23,      
    CPLEX_PRIMAL        = 24,      
    CPLEX_SIFTING       = 25,      
    CPLEX_NETWORK       = 26,
    CPLEX_NETWORK_2     = 27,
    CPLEX_CONCURRENT    = 28,
    CPLEX               = 30
  };



  LPSolver *createSolver(Optimizer optimType, double lambda){

    LPSolver *solver;
    if(optimType < GLPK){
      if(lambda > 0){
        Rprintf("Fuzzy match not supported in this optimizer\n");
      }
#ifdef MOP_USE_GLPK
      solver = new GLPKSolver();
#else
      Rprintf("GLPK not supported in this installation\n");
      return NULL;
#endif
    }
    else if( optimType < MSK){
      if(lambda > 0){
        Rprintf("Fuzzy match not supported in this optimizer\n");
      }
#ifdef MOP_USE_MOSEK
      MSKint32t optimizer = 0;
      switch(optimType){
        case MSK_FREE:
          optimizer = MSK_OPTIMIZER_FREE;
          break;
        case MSK_NETWORK_SIMPLEX:
          optimizer = MSK_OPTIMIZER_NETWORK_PRIMAL_SIMPLEX;
          break;
        case MSK_SIMPLEX_FREE:
          optimizer = MSK_OPTIMIZER_FREE_SIMPLEX;
          break;
        case MSK_INTPOINT:
          optimizer = MSK_OPTIMIZER_INTPNT;
          break;
        case MSK_CONCURRENT:
          optimizer = MSK_OPTIMIZER_CONCURRENT;
          break;
        default:
          break;
      }
      solver = new MOSEKSolver(optimizer);
#else
      Rprintf("MOSEK not supported in this installation\n");
      return NULL;
#endif
    }
    else if(optimType == CPLEX_NETWORK){
#ifdef MOP_USE_CPLEX
      solver = new CPLEXNetworkSolver(lambda);
#else
      Rprintf("CPLEX not supported in this installation\n");
      return NULL;
#endif
    }
    else{
#ifdef MOP_USE_CPLEX
      //IloCplex::Algorithm optimizer = IloCplex::NoAlg;
      int optimizer = CPLEX_AUTO;
      switch(optimType){
        case CPLEX_AUTO:
          //optimizer = IloCplex::AutoAlg;
          optimizer = CPX_ALG_AUTOMATIC;
          break;
        case CPLEX_BARRIER:
          //optimizer = IloCplex::Barrier;
          optimizer = CPX_ALG_BARRIER;
          break;
        case CPLEX_DUAL:
          //optimizer = IloCplex::Dual;
          optimizer = CPX_ALG_DUAL;
          break;
        case CPLEX_PRIMAL:
          //optimizer = IloCplex::Primal;
          optimizer = CPX_ALG_PRIMAL;
          break;
        case CPLEX_SIFTING:
          //optimizer = IloCplex::Sifting;
          optimizer = CPX_ALG_SIFTING;
          break;
        case CPLEX_NETWORK_2:
          //optimizer = IloCplex::Network;
          optimizer = CPX_ALG_NET;
          break;
        case CPLEX_CONCURRENT:
          //optimizer = IloCplex::Concurrent;
          optimizer = CPX_ALG_CONCURRENT;
          break;
        default:
          break;

      }
      solver = new CPLEXSolver(optimizer, lambda);
#else
      Rprintf("CPLEX not supported in this installation\n");
      return NULL;
#endif
    }


    return solver;
  }

  /**/  
  /* Multiscale Transport Calls */
  SEXP multiscaleTransport(GMRANeighborhood<double> *nh1,
      GMRANeighborhood<double> *nh2, int d1, int d2, double p, int nScales1, int
      nScales2, std::vector<double> &w1, std::vector<double> &w2,
      MultiscaleTransport<double> &transport, bool matchScale,
      bool multiscaleCost, bool multiscaleSolution){ 
    
    using namespace Eigen;
    
    typedef TransportPlan<double>::Path Path;


    std::vector< MultiscaleTransportLevel<double> * > t1Levels =
      GMRAMultiscaleTransportLevel<double>::buildTransportLevels(*nh1, w1,
          multiscaleCost);

    std::vector< MultiscaleTransportLevel<double> * > t2Levels =
      GMRAMultiscaleTransportLevel<double>::buildTransportLevels(*nh2, w2,
          multiscaleCost);

  

    std::vector< TransportPlan<double> * > sols = transport.solve( t1Levels,
        t2Levels, p, nScales1, nScales2, matchScale);


    VectorXd cost(sols.size());
    VectorXi vars(sols.size());

    for(unsigned int i=0; i < sols.size(); i++){
      TransportPlan<double> *s = sols[i];
      cost(i) = s->cost;
      vars(i) = s->getNumberOfPaths();
    }

    int nVars = 17;
    SEXP Rres;
    if(multiscaleSolution){
      PROTECT( Rres = Rf_allocVector( VECSXP, 5 + nVars*sols.size() ));
    }
    else{
      PROTECT( Rres = Rf_allocVector( VECSXP, 2 + nVars*sols.size() ));
    }
    
    SEXP Rcost;
    PROTECT(Rcost = Rf_allocVector(REALSXP, cost.size()));
    memcpy( REAL(Rcost), cost.data(), cost.size()*sizeof(double) );
    SET_VECTOR_ELT(Rres, 0, Rcost);

    SEXP Rvars;
    PROTECT(Rvars = Rf_allocVector(INTSXP, vars.size()));
    memcpy( INTEGER(Rvars), vars.data(), vars.size()*sizeof(int) );
    SET_VECTOR_ELT(Rres, 1, Rvars);
   





    
    for(unsigned int i=0; i < sols.size(); i++){
      TransportPlan<double> *s = sols[i];
    
      int n1 = s->source->getNodes().size();
      int n2 = s->target->getNodes().size();
      MatrixXd sPoints(d1, n1);
      MatrixXd tPoints(d2, n2);
      VectorXd fromMass(n1);
      VectorXd toMass(n2);

      typedef  TransportNode<double>::TransportNodeVector TransportNodeVector;
      typedef TransportNodeVector::const_iterator TransportNodeVectorCIterator;
      VectorXd pFrom = VectorXd::Zero(n1);
      VectorXd pTo = VectorXd::Zero(n2);
      
      for(TransportNodeVectorCIterator it = s->source->getNodes().begin(); it !=
          s->source->getNodes().end(); ++it){
        
        GMRATransportNode<double> *node = (GMRATransportNode<double> *)  *it;
        int fromIndex = node->getID();
        
        double pi = node->getPotential();
        pFrom(fromIndex) = pi;        
        
        sPoints.col(fromIndex)= node->getGMRANode()->getCenter();
        fromMass(fromIndex) = node->getMass();
        
      }

      for(TransportNodeVectorCIterator it = s->target->getNodes().begin(); it !=
          s->target->getNodes().end(); ++it){
        
        GMRATransportNode<double> *node = (GMRATransportNode<double> *) *it;
        int toIndex = node->getID();
        
        double pi = node->getPotential();
        pTo(toIndex) = pi;        
        
        tPoints.col(toIndex) = node->getGMRANode()->getCenter();
        
        toMass(toIndex) = node->getMass(); 
      }

      VectorXd costsFrom =  VectorXd::Zero( n1 );
      VectorXd costsTo = VectorXd::Zero( n2 );

      int nPaths = 0;
      for(s->pathIteratorBegin(); ! s->pathIteratorIsAtEnd();
          s->pathIteratorNext() ){
        Path &path = s->pathIteratorCurrent();
        GMRATransportNode<double> *from = (GMRATransportNode<double> *) path.from;
        GMRATransportNode<double> *to = (GMRATransportNode<double> *) path.to;
       
        if(path.w > 0){
          nPaths++;
          double c = path.cost * path.w;
          costsFrom( from->getID() ) += c;        
          costsTo( to->getID() ) += c;
        }
      }


      SEXP RcFrom; 
      PROTECT(RcFrom = Rf_allocVector(REALSXP, costsFrom.size()));
      memcpy( REAL(RcFrom), costsFrom.data(), costsFrom.size() * sizeof(double) );
      SET_VECTOR_ELT( Rres, 2 + nVars*i, RcFrom);

      SEXP RcTo; 
      PROTECT(RcTo = Rf_allocVector(REALSXP, costsTo.size()));
      memcpy( REAL(RcTo), costsTo.data(), costsTo.size() * sizeof(double) );
      SET_VECTOR_ELT( Rres, 2 + nVars*i + 1, RcTo);

      SEXP RcFromPot; 
      PROTECT( RcFromPot = Rf_allocVector(REALSXP, pFrom.size()) );
      memcpy( REAL(RcFromPot), pFrom.data(), pFrom.size() * sizeof(double) );
      SET_VECTOR_ELT( Rres, 2 + nVars*i + 2, RcFromPot);

      SEXP RcToPot; 
      PROTECT( RcToPot = Rf_allocVector(REALSXP, pTo.size()) );
      memcpy( REAL(RcToPot), pTo.data(), pTo.size() * sizeof(double) );
      SET_VECTOR_ELT( Rres, 2 + nVars*i + 3, RcToPot);



      MatrixXd map(nPaths, 4);
      int mapIndex = 0;
      for(s->pathIteratorBegin(); ! s->pathIteratorIsAtEnd();
          s->pathIteratorNext() ){
        Path &path = s->pathIteratorCurrent();
        if(path.w > 0 ){
          GMRATransportNode<double> *from = (GMRATransportNode<double> *)  path.from;
          GMRATransportNode<double> *to = (GMRATransportNode<double> *) path.to;

          map( mapIndex, 0 ) = from->getID() + 1;
          map( mapIndex, 1 ) = to->getID() + 1;
          map( mapIndex, 2 ) = path.w;
          map( mapIndex, 3 ) = path.cost;
          mapIndex++;
        }
      }


      SEXP Rmap;
      PROTECT( Rmap = Rf_allocMatrix(REALSXP, map.rows(), map.cols()));
      memcpy( REAL(Rmap), map.data(), map.rows()*map.cols()*sizeof(double) );
      SET_VECTOR_ELT( Rres, 2+nVars*i+4, Rmap); 

      SEXP RnTotal;
      int totalPaths = s->getNumberOfPaths();
      PROTECT( RnTotal = Rf_allocVector(INTSXP, 1) );
      memcpy( INTEGER(RnTotal), &totalPaths, sizeof(int) );
      SET_VECTOR_ELT(Rres, 2+nVars*i+5, RnTotal); 

      SEXP Rfrom;
      PROTECT(Rfrom = Rf_allocMatrix(REALSXP, sPoints.rows(), sPoints.cols()));
      memcpy( REAL(Rfrom), sPoints.data(), sPoints.rows()*sPoints.cols()*sizeof(double) );
      SET_VECTOR_ELT(Rres, 2+nVars*i+6, Rfrom); 

      SEXP Rto;
      PROTECT(Rto = Rf_allocMatrix(REALSXP, tPoints.rows(), tPoints.cols()));
      memcpy( REAL(Rto), tPoints.data(), tPoints.rows()*tPoints.cols()*sizeof(double) );
      SET_VECTOR_ELT(Rres, 2+nVars*i+7, Rto); 

      SEXP RfromMass;
      PROTECT(RfromMass = Rf_allocVector(REALSXP, fromMass.size()));
      memcpy( REAL(RfromMass), fromMass.data(), fromMass.size()*sizeof(double) );
      SET_VECTOR_ELT(Rres, 2+nVars*i+8, RfromMass); 

      SEXP RtoMass;
      PROTECT(RtoMass = Rf_allocVector(REALSXP, toMass.size()));
      memcpy( REAL(RtoMass), toMass.data(), toMass.size()*sizeof(double) );
      SET_VECTOR_ELT(Rres, 2+nVars*i+9, RtoMass);

      double time[3];
      time[0] = (double)s->timeSolve / CLOCKS_PER_SEC;
      time[1] = (double)s->timeRefine / CLOCKS_PER_SEC;
      time[2] = (double)s->timePropagate / CLOCKS_PER_SEC;
      SEXP Rtime;
      PROTECT( Rtime = Rf_allocVector(REALSXP, 3) );
      memcpy( REAL(Rtime), time, sizeof(double)*3 );
      SET_VECTOR_ELT(Rres, 2+nVars*i+10, Rtime); 



      std::vector<double> radiusFrom;
      std::vector<int> sizeFrom;
      std::vector<int> indexFrom;

      for (TransportNodeVectorCIterator it =
          s->source->getNodes().begin(); it != s->source->getNodes().end(); ++it) {
        GMRATransportNode<double> *node =
          (GMRATransportNode<double> *)  *it;
        std::vector<int> &pts = node->getGMRANode()->getPoints();
        sizeFrom.push_back(pts.size());
        radiusFrom.push_back( node->getNodeRadius() );
        for (unsigned int j = 0 ; j < pts.size() ; j++) {
          indexFrom.push_back(pts[j]+1);
        }
      }


      SEXP RindexFrom;
      PROTECT(RindexFrom = Rf_allocVector(INTSXP, indexFrom.size()));
      memcpy(INTEGER(RindexFrom), indexFrom.data(),
          indexFrom.size()*sizeof(int));
      SET_VECTOR_ELT(Rres, 2+nVars*i+11, RindexFrom);

      SEXP RsizeFrom;
      PROTECT(RsizeFrom = Rf_allocVector(INTSXP, sizeFrom.size()));
      memcpy(INTEGER(RsizeFrom), sizeFrom.data(),
          sizeFrom.size()*sizeof(int));
      SET_VECTOR_ELT(Rres, 2 + nVars*i+12, RsizeFrom);
      
      SEXP RradiusFrom;
      PROTECT(RradiusFrom = Rf_allocVector(REALSXP, radiusFrom.size()));
      memcpy( REAL(RradiusFrom), radiusFrom.data(),
          radiusFrom.size()*sizeof(double));
      SET_VECTOR_ELT(Rres, 2 + nVars*i+13, RradiusFrom);



      std::vector<double> radiusTo;
      std::vector<int> sizeTo;
      std::vector<int> indexTo;

      for (TransportNodeVectorCIterator it =
          s->target->getNodes().begin(); it != s->target->getNodes().end(); ++it) {
        GMRATransportNode<double> *node =
          (GMRATransportNode<double> *)  *it;
        std::vector<int> &pts = node->getGMRANode()->getPoints();
        sizeTo.push_back(pts.size());
        radiusTo.push_back( node->getNodeRadius() );
        for (unsigned int j = 0 ; j < pts.size() ; j++) {
          indexTo.push_back(pts[j]+1);
        }
      }


      SEXP RindexTo;
      PROTECT(RindexTo = Rf_allocVector(INTSXP, indexTo.size()));
      memcpy(INTEGER(RindexTo), indexTo.data(),
          indexTo.size()*sizeof(int));
      SET_VECTOR_ELT(Rres, 2+nVars*i+14, RindexTo);

      SEXP RsizeTo;
      PROTECT(RsizeTo = Rf_allocVector(INTSXP, sizeTo.size()));
      memcpy(INTEGER(RsizeTo), sizeTo.data(),
          sizeTo.size()*sizeof(int));
      SET_VECTOR_ELT(Rres, 2 + nVars*i+15, RsizeTo);

      
      SEXP RradiusTo;
      PROTECT( RradiusTo = Rf_allocVector(REALSXP, radiusTo.size()));
      memcpy( REAL(RradiusTo), radiusTo.data(),
              radiusTo.size()*sizeof(double));
      SET_VECTOR_ELT( Rres, 2 + nVars*i+16, RradiusTo );

    }



    //multiscale solution
    if(multiscaleSolution){

      TransportPlan<double> *finest = sols[sols.size()-1];
      int nPaths = 0;
      for(finest->pathIteratorBegin(); ! finest->pathIteratorIsAtEnd();
          finest->pathIteratorNext() ){
        Path &path = finest->pathIteratorCurrent();
        if(path.w > 0){
          nPaths++;
        }
      }

      MatrixXi pFrom(t1Levels.size(), nPaths );
      MatrixXi pTo(t2Levels.size(), nPaths );

      int pId = 0;
      for(finest->pathIteratorBegin(); ! finest->pathIteratorIsAtEnd();
          finest->pathIteratorNext() ){

        
        Path &path = finest->pathIteratorCurrent();
        if(path.w > 0){       
          TransportNode<double> *from = (GMRATransportNode<double> *)  path.from;
          int fIndex= pFrom.rows()-1;
          while(from != NULL){
            pFrom(fIndex, pId) = from->getID(); 
            from = from->getParent();
            fIndex--;
          }

          TransportNode<double> *to = (GMRATransportNode<double> *) path.to;
          int tIndex= pTo.rows()-1;
          while(to != NULL){
            pTo(tIndex, pId) = to->getID(); 
            to = to->getParent();
            tIndex--;
          }
          ++pId;

        }
        
      }

      int index = nVars*sols.size() + 2;

      SEXP RpFrom;
      PROTECT(RpFrom = Rf_allocMatrix(INTSXP, pFrom.rows(), pFrom.cols()));
      memcpy( INTEGER(RpFrom), pFrom.data(), pFrom.rows()*pFrom.cols()*sizeof(int) );
      SET_VECTOR_ELT(Rres, index, RpFrom);
   
      SEXP RpTo;
      PROTECT(RpTo = Rf_allocMatrix(INTSXP, pTo.rows(), pTo.cols()));
      memcpy( INTEGER(RpTo), pTo.data(), pTo.rows()*pTo.cols()*sizeof(int) );
      SET_VECTOR_ELT(Rres, index+1, RpTo);
      
     
      std::vector<double> mCost = finest->getMultiscaleTransportCost(p);
      SEXP RmCost;
      PROTECT(RmCost = Rf_allocVector(REALSXP, mCost.size() ) );
      memcpy( REAL(RmCost), mCost.data(), mCost.size()*sizeof(double) );
      SET_VECTOR_ELT(Rres, index+2, RmCost);



      UNPROTECT( 6 + nVars*sols.size() );

    }
    else{
      UNPROTECT( 3 + nVars*sols.size() );
    }

    for(unsigned int i=0; i < sols.size(); i++){
      delete sols[i];
    }
    
    for(unsigned int i=0; i < t1Levels.size(); i++){
      delete t1Levels[i];
    }
    for(unsigned int i=0; i < t2Levels.size(); i++){
      delete t2Levels[i];
    }



    return Rres;

  };









  NodeDistance<double> *getNodeDistance(SEXP RdType){
    int distType = *INTEGER(RdType);
    NodeDistance<double> *dist = NULL;
    if(distType == 1){
      dist = new CenterNodeDistance<double>( new EuclideanMetric<double>() );
    }
    else if(distType == 2){
      dist = new CenterNodeDistance<double>( new L1Metric<double>() );
    }
    else if(distType == 3){
      dist = new CenterNodeDistance<double>( new SquaredEuclideanMetric<double>() );
    }
    else if(distType == 4){
      dist = new WassersteinNodeDistance<double>();
    }
    else{
      dist = new CenterNodeDistance<double>( new EuclideanMetric<double>() );
    }
    return dist;
  };











  GMRANeighborhood<double> *getGMRANeighborhood(SEXP Rgmra, SEXP RnType, SEXP RdType){

    SEXP Rgmra_pointer = Rf_getAttrib(Rgmra, Rf_install("gmra_ptr") );
    GMRATree<double> *gmra = static_cast<GMRATree<double> *>( R_ExternalPtrAddr(Rgmra_pointer) );
    

    NodeDistance<double> *dist = getNodeDistance(RdType);
    gmra->computeRadii(dist);
    gmra->computeLocalRadii(dist);

    int nType = *INTEGER(RnType);
    if(nType == 2){
      return new RelativeGMRANeighborhood<double>(gmra, dist);
    }
    else{
      return new GenericGMRANeighborhood<double>(gmra, dist);
    }
  };




  //
  MultiscaleTransportLP<double> *getTransportLP(SEXP Rtrp){    
    
    SEXP Rtrp_pointer = Rf_getAttrib(Rtrp, Rf_install("mop_trp_lp_ptr") ); 
    
    MultiscaleTransportLP<double> *trp =
      static_cast<MultiscaleTransportLP<double> *>(
          R_ExternalPtrAddr(Rtrp_pointer) );
    
    return trp;
  
  };



  //
  SEXP createTransportLP(SEXP RoType, SEXP Rlambda){
 
    double lambda = *REAL(Rlambda);
    Optimizer optimizer = (Optimizer)  * INTEGER(RoType) ;

    LPSolver *solver = createSolver(optimizer, lambda);
    MultiscaleTransportLP<double> *transport = new MultiscaleTransportLP<double>(solver);
      
    SEXP Rid;
    PROTECT(Rid = Rf_allocVector(INTSXP, 1));
    INTEGER(Rid)[0] = trpCounter;
    SEXP ptr = R_MakeExternalPtr(transport, Rf_install("mop"), R_NilValue);
    PROTECT(ptr);
    Rf_setAttrib(Rid, Rf_install("mop_trp_lp_ptr"), ptr);
    
    ++trpCounter;

    UNPROTECT(2);
    return Rid; 
  };




  //Neighborhood strategies
  SEXP addExpandNeighborhoodStrategy(SEXP Rtrp, SEXP Rfactor, SEXP Rtol, SEXP RnIter, SEXP RmaxAdd){
    double factor = *REAL(Rfactor);
    double tolerance = *REAL(Rtol);
    int nIter = *INTEGER(RnIter);
    int maxAdd = *INTEGER(RmaxAdd);

    ExpandNeighborhoodStrategy<double> *expand = new
      ExpandNeighborhoodStrategy<double>(factor, tolerance, nIter, maxAdd);

    MultiscaleTransportLP<double> *trp = getTransportLP(Rtrp);
    trp->addNeighborhodStrategy(expand);

    return R_NilValue;
  };


  SEXP addRefineNeighborhoodStrategy(SEXP Rtrp, SEXP Rfactor, SEXP Rtol, SEXP RnIter, SEXP RmaxAdd){
    double factor = *REAL(Rfactor);
    double tolerance = *REAL(Rtol);
    int nIter = *INTEGER(RnIter);
    int maxAdd = *INTEGER(RmaxAdd);

    RefineNeighborhoodStrategy<double> *refine = new
      RefineNeighborhoodStrategy<double>(factor, tolerance, nIter, maxAdd);

    MultiscaleTransportLP<double> *trp = getTransportLP(Rtrp);
    trp->addNeighborhodStrategy(refine);

    return R_NilValue;
  };


  
  SEXP addPotentialNeighborhoodStrategy(SEXP Rtrp, SEXP Rthres,  SEXP Rtol, SEXP RnIter, SEXP Rsort, SEXP Rexpand, SEXP RmaxAdd){
    double threshold = *REAL(Rthres);
    double tolerance = *REAL(Rtol);
    int nIter = *INTEGER(RnIter);
    bool sort = *INTEGER(Rsort);
    bool expand = *INTEGER(Rexpand);
    int maxAdd = *INTEGER(RmaxAdd);

    PotentialNeighborhoodStrategy<double> *refine = new
      PotentialNeighborhoodStrategy<double>(threshold, tolerance, sort, expand,
          nIter, maxAdd);

    MultiscaleTransportLP<double> *trp = getTransportLP(Rtrp);
    trp->addNeighborhodStrategy(refine);

    return R_NilValue;
  };


  //Propagation strategies
  PropagationStrategy<double> *getPropagationStrategy(SEXP Rprop){    
    
    SEXP Rprop_pointer = Rf_getAttrib(Rprop, Rf_install("mop_trp_prop_ptr") ); 
    
    PropagationStrategy<double> *prop = static_cast<
      PropagationStrategy<double>* > ( R_ExternalPtrAddr(Rprop_pointer) );
    
    return prop;
  };



  //
  SEXP createNeighborhoodPropagationStrategy(SEXP Rfactor){
 
    double factor = *REAL(Rfactor);
     
    NeighborhoodPropagationStrategy<double> *prop = new NeighborhoodPropagationStrategy<double>(factor); 
    
    SEXP Rid;
    PROTECT(Rid = Rf_allocVector(INTSXP, 1));
    INTEGER(Rid)[0] = propCounter;
    SEXP ptr = R_MakeExternalPtr(prop, Rf_install("mop"), R_NilValue);
    PROTECT(ptr);
    Rf_setAttrib(Rid, Rf_install("mop_trp_prop_ptr"), ptr);
    
    ++propCounter;

    UNPROTECT(2);
    return Rid; 
  };

  SEXP createCapacityPropagationStrategy(SEXP Rk, SEXP Rfactor){
 
    double factor = *REAL(Rfactor);
    int k = *INTEGER(Rk);
     
    CapacityPropagationStrategy<double> *prop = new CapacityPropagationStrategy<double>(k, factor); 
    
    SEXP Rid;
    PROTECT(Rid = Rf_allocVector(INTSXP, 1));
    INTEGER(Rid)[0] = propCounter;
    SEXP ptr = R_MakeExternalPtr(prop, Rf_install("mop"), R_NilValue);
    PROTECT(ptr);
    Rf_setAttrib(Rid, Rf_install("mop_trp_prop_ptr"), ptr);
    
    ++propCounter;

    UNPROTECT(2);
    return Rid; 
  };

    SEXP createIteratedCapacityPropagationStrategy(SEXP Rk, SEXP Rfactor){
 
    double factor = *REAL(Rfactor);
    int k = *INTEGER(Rk);
     
    IteratedCapacityPropagationStrategy<double> *prop = new
	    IteratedCapacityPropagationStrategy<double>(k, factor); 
    
    SEXP Rid;
    PROTECT(Rid = Rf_allocVector(INTSXP, 1));
    INTEGER(Rid)[0] = propCounter;
    SEXP ptr = R_MakeExternalPtr(prop, Rf_install("mop"), R_NilValue);
    PROTECT(ptr);
    Rf_setAttrib(Rid, Rf_install("mop_trp_prop_ptr"), ptr);
    
    ++propCounter;

    UNPROTECT(2);
    return Rid; 
  };

  SEXP createRandomizedNeighborhoodPropagationStrategy(SEXP RnRandom, SEXP Rfactor){
 
    int nRadnom = *INTEGER( RnRandom );
    double factor = *REAL( Rfactor );
     
    RandomizedNeighborhoodPropagationStrategy<double> *prop = new
      RandomizedNeighborhoodPropagationStrategy<double>(nRadnom, factor); 
    
    SEXP Rid;
    PROTECT(Rid = Rf_allocVector(INTSXP, 1));
    INTEGER(Rid)[0] = propCounter;
    SEXP ptr = R_MakeExternalPtr(prop, Rf_install("mop"), R_NilValue);
    PROTECT(ptr);
    Rf_setAttrib(Rid, Rf_install("mop_trp_prop_ptr"), ptr);
    
    ++propCounter;

    UNPROTECT(2);
    return Rid; 
  };


  SEXP createMaxEntropyPropagationStrategy(){
     
    MaxEntropyPropagationStrategy<double> *prop = new MaxEntropyPropagationStrategy<double>(); 
    
    SEXP Rid;
    PROTECT(Rid = Rf_allocVector(INTSXP, 1));
    INTEGER(Rid)[0] = propCounter;
    SEXP ptr = R_MakeExternalPtr(prop, Rf_install("mop"), R_NilValue);
    PROTECT(ptr);
    Rf_setAttrib(Rid, Rf_install("mop_trp_prop_ptr"), ptr);
    
    ++propCounter;

    UNPROTECT(2);
    return Rid; 
  };



  SEXP createSinkhornPropagationStrategy(SEXP Rlambda, SEXP Rtol, SEXP Rthres, SEXP Riter){
    double lambda = *REAL(Rlambda);
    double threshold = *REAL(Rthres);
    double tolerance = *REAL(Rtol);
    int nIter = *INTEGER(Riter);


     
    SinkhornPropagationStrategy<double> *prop = new
      SinkhornPropagationStrategy<double>(lambda, tolerance, threshold, nIter); 
    
    SEXP Rid;
    PROTECT(Rid = Rf_allocVector(INTSXP, 1));
    INTEGER(Rid)[0] = propCounter;
    SEXP ptr = R_MakeExternalPtr(prop, Rf_install("mop"), R_NilValue);
    PROTECT(ptr);
    Rf_setAttrib(Rid, Rf_install("mop_trp_prop_ptr"), ptr);
    
    ++propCounter;

    UNPROTECT(2);
    return Rid; 
  };



  SEXP setPropagationStrategy1(SEXP Rtrp, SEXP Rprop){
    PropagationStrategy<double> *prop = getPropagationStrategy(Rprop);
    MultiscaleTransportLP<double> *trp = getTransportLP(Rtrp);
    trp->setPropagationStrategy1(prop);

    return R_NilValue;
  };


  SEXP setPropagationStrategy2(SEXP Rtrp, SEXP Rprop){
    PropagationStrategy<double> *prop = getPropagationStrategy(Rprop);
    MultiscaleTransportLP<double> *trp = getTransportLP(Rtrp);
    trp->setPropagationStrategy2(prop);

    return R_NilValue;
  };


  SEXP setMaxNeighborhoodSize(SEXP Rtrp, SEXP Rn){
    int n = *INTEGER(Rn);
    MultiscaleTransportLP<double> *trp = getTransportLP(Rtrp);
    trp->setMaxNeighborhoodSize(n);

    return R_NilValue;
  };


  //solve the multiscale transport problem
  SEXP multiscaleTransportSolve(SEXP Rtrp, SEXP Rgmra1, SEXP Rgmra2, SEXP Rs1, SEXP Rs2,
      SEXP Rw1, SEXP Rnw1, SEXP Rw2, SEXP Rnw2, SEXP Rp, SEXP RmatchScale, SEXP
      RmultiscaleCost, SEXP RmultiscaleSolution, SEXP RdType, SEXP RnType){



    GMRANeighborhood<double> *nh1 = getGMRANeighborhood(Rgmra1, RnType, RdType);
    GMRANeighborhood<double> *nh2 = getGMRANeighborhood(Rgmra2, RnType, RdType);

    int m1 = nh1->getTree()->getRoot()->getCenter().size();
    int m2 = nh2->getTree()->getRoot()->getCenter().size();

    int n1 = nh1->getTree()->getRoot()->getPoints().size();
    int n2 = nh2->getTree()->getRoot()->getPoints().size();

    int s1 = *INTEGER(Rs1);
    int s2 = *INTEGER(Rs2);


    double *w1 = REAL(Rw1);
    double *w2 = REAL(Rw2);

    int nw1 = *INTEGER(Rnw1);
    int nw2 = *INTEGER(Rnw2);

    bool matchScale = *INTEGER(RmatchScale);
    bool multiscaleCost = *INTEGER(RmultiscaleCost);
    bool multiscaleSolution = *INTEGER(RmultiscaleSolution);

    std::vector<double> weights1;
    if(nw1 == n1){
      weights1.resize(nw1);
      for(int i=0; i<nw1; i++){
        weights1[i] = w1[i];
      }
    }
    std::vector<double> weights2;
    if(nw2 == n2){
      weights2.resize(nw2);
      for(int i=0; i<nw2; i++){
        weights2[i] = w2[i];
      }
    }



    double p = *REAL(Rp);

    MultiscaleTransportLP<double> *transport = getTransportLP(Rtrp);


    SEXP res = multiscaleTransport(nh1, nh2, m1, m2, p, s1, s2, weights1,
        weights2, *transport, matchScale, multiscaleCost, multiscaleSolution); 

    return res;

  };












  //Single scale transport for a given cost matrix c between densities x1 and
  //x2
  SEXP transport(SEXP Rx1, SEXP Rn1, SEXP Rx2, SEXP Rn2, SEXP Rc, SEXP RoType, SEXP Rlambda){
    using namespace Eigen;

    int n1 = *INTEGER(Rn1);
    int n2 = *INTEGER(Rn2);
    double *x1 = REAL(Rx1);
    double *x2 = REAL(Rx2);
    double *c = REAL(Rc);
    double lambda = *REAL(Rlambda);

    Map<MatrixXd> C(c, n1, n2);
    Map<VectorXd> from(x1, n1);
    Map<VectorXd> to(x2, n2);

    Optimizer optimizer = (Optimizer)  * INTEGER(RoType) ;
    LPSolver *solver = createSolver(optimizer, lambda);

    TransportLP<double> tp(solver);


    std::map< std::pair<int, int>, double> plan = tp.solve(C, from, to);
    MatrixXd map(plan.size(), 3);
    int mapIndex = 0;

    double cost = 0;
    for(std::map< std::pair<int, int>, double>::iterator it = plan.begin(); it!=plan.end(); ++it){

      const std::pair<int, int> &p = it->first;
      map( mapIndex, 0 ) = p.first;
      map( mapIndex, 1 ) = p.second;
      map( mapIndex, 2 ) = it->second;;

      cost += C(p.first, p.second) * it->second;

      mapIndex++;
    }


    SEXP Rres;
    PROTECT( Rres = Rf_allocVector( VECSXP, 2) );


    SEXP Rmap;
    PROTECT( Rmap = Rf_allocMatrix(REALSXP, map.rows(), map.cols()));
    memcpy( REAL(Rmap), map.data(), map.rows()*map.cols()*sizeof(double) );
    SET_VECTOR_ELT( Rres, 0, Rmap );

    SEXP Rcost;
    PROTECT( Rcost = Rf_allocVector(REALSXP, 1) );
    memcpy( REAL(Rcost), &cost, sizeof(double) );
    SET_VECTOR_ELT( Rres, 1, Rcost );


    UNPROTECT(3);


    delete solver;

    return Rres;

  };





}//end extern C
