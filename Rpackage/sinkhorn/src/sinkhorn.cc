#ifndef NULL
#define NULL 0
#endif

#define R_NO_REMAP


#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdio.h>


#include "IPCATree.h"
#include "GMRATree.h"
#include "GMRANeighborhood.h"
#include "RelativeGMRANeighborhood.h"

#include "GMRAMultiscaleTransport.h"
#include "MultiscaleSinkhornTransport.h"
#include "NodeDistance.h"
#include "EigenL1Metric.h"
#include "EigenSquaredEuclideanMetric.h"
#include "EigenEuclideanMetric.h"

#include "SinkhornTransport.h"



extern "C" {


  /**/  
  /* Multiscale Transport Calls */
  SEXP multiscaleSinkhornTransport(GMRATree<double> *t1, GMRATree<double> *t2, int d1,
      int d2, double p, int nScales1, int nScales2, std::vector<double> &w1,
      std::vector<double> &w2, SinkhornParameters<double> &params, int distType,
      int nType, bool multiscale){ 
    
    using namespace Eigen;
    
    typedef  TransportPlan<double>::Path Path;


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

    t1->computeRadii(dist);
    t1->computeLocalRadii(dist);

    t2->computeRadii(dist);
    t2->computeLocalRadii(dist);


    GMRANeighborhood<double> *nh1;
    GMRANeighborhood<double> *nh2;
    if(nType == 0){
      nh1 = new GenericGMRANeighborhood<double>(t1, dist);
      nh2 = new GenericGMRANeighborhood<double>(t2, dist);
    }
    else{
      nh1 = new RelativeGMRANeighborhood<double>(t1, dist);
      nh2 = new RelativeGMRANeighborhood<double>(t2, dist);
    }



    std::cout << "building transport levels" << std::endl;
    std::vector< MultiscaleTransportLevel<double> * > t1Levels =
      GMRAMultiscaleTransportLevel<double>::buildTransportLevels(*nh1, w1, multiscale);

    std::cout << "building transport levels" << std::endl;
    std::vector< MultiscaleTransportLevel<double> * > t2Levels =
      GMRAMultiscaleTransportLevel<double>::buildTransportLevels(*nh2, w2, multiscale);


    MultiscaleSinkhornTransport<double> *transport = new
      MultiscaleSinkhornTransport<double>(params);

    std::vector< TransportPlan<double> * > sols = transport->solve( t1Levels,
        t2Levels, p, nScales1, nScales2, false);


    VectorXd cost(sols.size());
    VectorXi vars(sols.size());

    for(unsigned int i=0; i < sols.size(); i++){
      TransportPlan<double> *s = sols[i];
      cost(i) = s->cost;
      vars(i) = s->getNumberOfPaths();
    }

    SEXP Rres;
    PROTECT( Rres = Rf_allocVector( VECSXP, 2 + 10*sols.size() ));

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
      VectorXd pTo= VectorXd::Zero(n2);
      
      for(TransportNodeVectorCIterator it = s->source->getNodes().begin(); it !=
          s->source->getNodes().end(); ++it){
        
        GMRATransportNode<double> *node = (GMRATransportNode<double> *)  *it;
        int fromIndex = node->getID();
        
        double pi = node->getPotential();
        pFrom(fromIndex) = pi;        
        
        sPoints.col(fromIndex) = node->getGMRANode()->getCenter();
        
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

      VectorXd costsFrom = VectorXd::Zero( n1 );
      VectorXd costsTo = VectorXd::Zero( n2 );

      int nPaths = 0;
      for(s->pathIteratorBegin(); ! s->pathIteratorIsAtEnd();
          s->pathIteratorNext() ){
        Path &path = s->pathIteratorCurrent();
        GMRATransportNode<double> *from = (GMRATransportNode<double> *)  path.from;
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
      SET_VECTOR_ELT( Rres, 2 + 10*i, RcFrom);

      SEXP RcTo; 
      PROTECT(RcTo = Rf_allocVector(REALSXP, costsTo.size()));
      memcpy( REAL(RcTo), costsTo.data(), costsTo.size() * sizeof(double) );
      SET_VECTOR_ELT( Rres, 2 + 10*i + 1, RcTo);

      SEXP RcFromPot; 
      PROTECT( RcFromPot = Rf_allocVector(REALSXP, pFrom.size()) );
      memcpy( REAL(RcFromPot), pFrom.data(), pFrom.size() * sizeof(double) );
      SET_VECTOR_ELT( Rres, 2 + 10*i + 2, RcFromPot);

      SEXP RcToPot; 
      PROTECT( RcToPot = Rf_allocVector(REALSXP, pTo.size()) );
      memcpy( REAL(RcToPot), pTo.data(), pTo.size() * sizeof(double) );
      SET_VECTOR_ELT( Rres, 2 + 10*i + 3, RcToPot);


      MatrixXd map(nPaths, 3);
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
          mapIndex++;
        }
      }


      SEXP Rmap;
      PROTECT( Rmap = Rf_allocMatrix(REALSXP, map.rows(), map.cols()));
      memcpy( REAL(Rmap), map.data(), map.rows()*map.cols()*sizeof(double) );
      SET_VECTOR_ELT( Rres, 2+10*i+4, Rmap); 


      int nTotal = s->getNumberOfPaths();
      SEXP RnTotal;
      PROTECT( RnTotal = Rf_allocVector(INTSXP, 1) );
      memcpy( INTEGER(RnTotal), &nTotal, sizeof(int) );
      SET_VECTOR_ELT(Rres, 2+10*i+5, RnTotal); 

      SEXP Rfrom;
      PROTECT(Rfrom = Rf_allocMatrix(REALSXP, sPoints.rows(), sPoints.cols()));
      memcpy( REAL(Rfrom), sPoints.data(), sPoints.rows()*sPoints.cols()*sizeof(double) );
      SET_VECTOR_ELT(Rres, 2+10*i+6, Rfrom); 

      SEXP Rto;
      PROTECT(Rto = Rf_allocMatrix(REALSXP, tPoints.rows(), tPoints.cols()));
      memcpy( REAL(Rto), tPoints.data(), tPoints.rows()*tPoints.cols()*sizeof(double) );
      SET_VECTOR_ELT(Rres, 2+10*i+7, Rto); 


      SEXP RfromMass;
      PROTECT(RfromMass = Rf_allocVector(REALSXP, fromMass.size()));
      memcpy( REAL(RfromMass), fromMass.data(), fromMass.size()*sizeof(double) );
      SET_VECTOR_ELT(Rres, 2+10*i+8, RfromMass); 

      SEXP RtoMass;
      PROTECT(RtoMass = Rf_allocVector(REALSXP, toMass.size()));
      memcpy( REAL(RtoMass), toMass.data(), toMass.size()*sizeof(double) );
      SET_VECTOR_ELT(Rres, 2+10*i+9, RtoMass);





    }


    UNPROTECT( 3 + 10*sols.size() );




    for(unsigned int i=0; i < sols.size(); i++){
      delete sols[i];
    }
    
    for(unsigned int i=0; i < t1Levels.size(); i++){
      delete t1Levels[i];
    }
    for(unsigned int i=0; i < t2Levels.size(); i++){
      delete t2Levels[i];
    }

    delete transport;
    delete nh1;
    delete nh2;
    delete dist;


    return Rres;

  };






  SEXP multiscale_sinkhorn_transport(SEXP Rgmra1, SEXP Rgmra2, SEXP Rs1, SEXP
      Rs2, SEXP Rw1, SEXP Rnw1, SEXP Rw2, SEXP Rnw2, SEXP Rp,  SEXP RdType,
      SEXP RnType, SEXP Rlambda, SEXP Rtolerance, SEXP RmaxIter, SEXP
      Rthreshold, SEXP RmaxPathsPerNode, SEXP RmsCost){
        
    SEXP Rgmra_pointer1 = Rf_getAttrib(Rgmra1, Rf_install("gmra_ptr") );
    GMRATree<double> *gmra1 = static_cast<GMRATree<double> *>( R_ExternalPtrAddr(Rgmra_pointer1) );\
    
    SEXP Rgmra_pointer2 = Rf_getAttrib(Rgmra2, Rf_install("gmra_ptr") );
    GMRATree<double> *gmra2 = static_cast<GMRATree<double> *>( R_ExternalPtrAddr(Rgmra_pointer2) );\



    int m1 = gmra1->getRoot()->getCenter().size();
    int m2 = gmra2->getRoot()->getCenter().size();

    int n1 = gmra2->getRoot()->getPoints().size();
    int n2 = gmra2->getRoot()->getPoints().size();
    
    double *w1 = REAL(Rw1);
    int nw1 = *INTEGER(Rnw1);

    double *w2 = REAL(Rw2);
    int nw2 = *INTEGER(Rnw2);
    bool msCost = *INTEGER(RmsCost);



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
    



    SinkhornParameters<double> config;
    config.iterations = *INTEGER(RmaxIter);
    config.maxPathsPerNode = *INTEGER(RmaxPathsPerNode);
    config.tolerance = *REAL(Rtolerance);
    config.lambda = *REAL(Rlambda);
    config.threshold = *REAL(Rthreshold);
    
    int dType = *INTEGER(RdType);
    int nType= *INTEGER(RnType);
    double p = *REAL(Rp);
    int s1 = *INTEGER(Rs1);
    int s2 = *INTEGER(Rs2);
    return multiscaleSinkhornTransport(gmra1, gmra2, m1, m2, p, s1, s2,
        weights1, weights2, config, dType, nType, msCost); 

  };
     




  SEXP sinkhorn_transport(SEXP Rmu, SEXP Rn1, SEXP Rnu, SEXP Rn2, SEXP Rm2, SEXP
      RD, SEXP Rlambda, SEXP Rtol, SEXP RmaxIter){
    using namespace Eigen;
    
    int n1 = *INTEGER(Rn1);
    int n2 = *INTEGER(Rn2);
    int m2 = *INTEGER(Rm2);
    double *muData = REAL(Rmu);
    double *nuData = REAL(Rnu);
    double *mData = REAL(RD);
    
    double lambda = *REAL(Rlambda);
    double tol = *REAL(Rtol);
    int iter = *INTEGER(RmaxIter);


    Map<VectorXd> muTmp(muData, n1); 
    VectorXd mu(n1);
    for(int i=0; i<n1; i++){
      mu(i) = muTmp(i);
    }

    Map<MatrixXd> nuTmp(nuData,n2,m2);
    MatrixXd nu(nuTmp.rows(), nuTmp.cols());
    for(int i=0; i<nuTmp.rows(); i++){ 
      for(int j=0; j<nuTmp.cols(); j++){ 
        nu(i, j) = nuTmp(i, j); 
      } 
    }

    Map<MatrixXd> M(mData, n1,n2 );
    MatrixXd K(M.rows(), M.cols());
    MatrixXd U(M.rows(), M.cols());
    for(int i=0; i<M.rows(); i++){ 
      for(int j=0; j<M.cols(); j++){ 
        double e = exp(-lambda * M(i, j)); 
        K(i, j) = e;
        U(i, j) = e*M(i,j);
      } 
    }

    clock_t t0 = clock();
    SinkhornTransport sinkhorn;
    sinkhorn.transport(mu, nu, K, U, tol, iter);
    clock_t t1 = clock();

    
    SEXP Rres;
    PROTECT( Rres = Rf_allocVector( VECSXP, 2 ));


    SEXP Rcosts; 
    PROTECT( Rcosts = Rf_allocVector( REALSXP, nu.cols() ) );
    memcpy( REAL(Rcosts), sinkhorn.getDistances().data(), nu.cols() * sizeof(double) );
    SET_VECTOR_ELT( Rres, 0, Rcosts);

    SEXP Rtime; 
    PROTECT(Rtime = Rf_allocVector(REALSXP, 1) );
    double time = ((double) t1 - t0) / CLOCKS_PER_SEC;
    memcpy( REAL(Rtime), &time,  sizeof(double) );
    SET_VECTOR_ELT( Rres, 1, Rtime);

    UNPROTECT(3);

    return Rres;


  }; 

   

}//end extern C
