#ifndef NULL
#define NULL 0
#endif

#define R_NO_REMAP

//#define VERBOSE

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>


#include <stdio.h>


#include "IPCATree.h"
#include "IKMTree.h"

#include "GMRATree.h"
#include "GMRAPrune.h"

#include "GMRADensityEstimator.h"
#include "MeanGMRAPredictor.h"
#include "GMRANeighborhood.h"
#include "EigenL1Metric.h"
#include "EigenEuclideanMetric.h"
#include "EigenSquaredEuclideanMetric.h"
#include "WassersteinNodeDistance.h"



extern "C" {
  //GMRA tree storage
  //std::vector< GMRATree<double> *> gmras; 
  int idCounter;
  int nodeCounter;
  
  
  /* called from .onLoad */
  SEXP gmra_init(void){
    
    idCounter = 0;  
    nodeCounter = 0;  
    return R_NilValue;

  };



  //Get node distance type
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
    return dist;
  };






  //Create IPCA Tree setup methods
  IPCANodeFactory<double> *getNodeFactory(GMRADataObject<double> *D, int d, double t){
    if(t < 0){
      return new FixedNodeFactory<double>(D, d);
    }
    else if(t < 1){
      return new RelativePrecisionNodeFactory<double>(D, d, t);
    }
    else{
      return new RelativeRatioNodeFactory<double>(D, d, t);
    }
  };



  IPCANodeFactory<double> *createIPCAFactory(GMRADataObject<double> *D, SEXP Rt, SEXP Reps, SEXP RmaxKids,
      SEXP Rsplit, SEXP RsplitDir, SEXP Rstop, SEXP Rd, SEXP RminPoints){ 

    IPCANodeFactory<double> *factory = getNodeFactory(D, *INTEGER(Rd), *REAL(Rt) );
    
    factory->epsilon = *REAL(Reps);
    factory->maxKidDim = *INTEGER(RmaxKids);
    factory->setSplitStrategy( *INTEGER(Rsplit) );
    factory->setSplitDirectionStrategy( *INTEGER(RsplitDir) );
    factory->setStoppingCriterium(*INTEGER(Rstop) );
    factory->minPoints = *INTEGER(RminPoints);

    return factory;
  };




  //Create IPCA trees, if splitStartegy/stopStrategy is random it can be used to
  //build multiple trees on the same data without increasing storage beyond the
  //tree structures
  SEXP gmra_create_ipca(SEXP Rx, SEXP Rn, SEXP Rm, SEXP Rd, SEXP Rt, SEXP Reps,
      SEXP RmaxKids, SEXP Rsplit, SEXP RsplitDir, SEXP Rstop, SEXP RnRuns, SEXP RminPoints){


    int m = *INTEGER(Rm);
    int n = *INTEGER(Rn);
    double *x = REAL(Rx);

    Eigen::MatrixXd X = Eigen::Map<Eigen::MatrixXd>(x,m,n);

    MatrixGMRADataObject<double> *data = new MatrixGMRADataObject<double>(X);

    int nRuns = *INTEGER(RnRuns);
    
    SEXP list;
    PROTECT( list = Rf_allocVector(VECSXP, nRuns) );
   
    std::vector<int> pts(X.cols() );
    for(unsigned int i=0; i<pts.size(); i++){
      pts[i] = i;
    };
    for(int i=0; i<nRuns; i++){
      IPCANodeFactory<double> *factory = 
        createIPCAFactory(data, Rt, Reps, RmaxKids, Rsplit, RsplitDir, Rstop, Rd, RminPoints); 
      
      IPCATree<double> *gmra = new IPCATree<double>(data, factory);

      gmra->addPoints(pts);


      SEXP Rid;
      PROTECT(Rid = Rf_allocVector(INTSXP, 1));
      INTEGER(Rid)[0] = idCounter;
      SEXP ptr = R_MakeExternalPtr(gmra, Rf_install("GMRA"), R_NilValue);
      PROTECT(ptr);
      Rf_setAttrib(Rid, Rf_install("gmra_ptr"), ptr); 
      SET_VECTOR_ELT(list, i, Rid);

      ++idCounter;
    }

    UNPROTECT(2*nRuns+1);

    return list;

  };




  IKMKmeansDataFactory<double> *getIKMDataFactory(int similarity, int nCor, int nSpatial){
    std::cout << similarity << std::endl;
    
    switch(similarity){
      case 1:
       return new L2GMRAKmeansDataFactory<double>();
      case 2:
       return new CorrelationGMRAKmeansDataFactory<double>();
      case 3:
       return new JointCorrelationAndSpatialGMRAKmeansDataFactory<double>(nCor, nSpatial);
    }
       
    return new L2GMRAKmeansDataFactory<double>();

  };


  //Create IKmeans trees, since Kmans typically finds a local minimum multiple
  //trees can be cerated on the same data without increasing storage beyond the
  //tree structures
  SEXP gmra_create_ikm(SEXP Rx, SEXP Rn, SEXP Rm, SEXP Rt, SEXP Reps, SEXP
      RnKids, SEXP RmaxIter, SEXP Rsplit, SEXP Rstop, SEXP RnRuns, SEXP
      Rsimilarity, SEXP RnSpatial, SEXP RnCor, SEXP RminPoints){


    int m = *INTEGER(Rm);
    int n = *INTEGER(Rn);
    double *x = REAL(Rx);

    
    int nSpatial = *INTEGER(RnSpatial);
    int nCor = *INTEGER(RnCor);

    Eigen::MatrixXd X = Eigen::Map<Eigen::MatrixXd>(x,m,n);

    MatrixGMRADataObject<double> *data = new MatrixGMRADataObject<double>(X);

    int nRuns = *INTEGER(RnRuns);
    
    SEXP list;
    PROTECT( list = Rf_allocVector(VECSXP, nRuns) );
   
    std::vector<int> pts(X.cols() );
    for(unsigned int i=0; i<pts.size(); i++){
      pts[i] = i;
    };
    for(int i=0; i<nRuns; i++){
      
      IKMTree<double> *gmra = new IKMTree<double>(data);
      gmra->setStoppingCriterium( *INTEGER(Rstop) );
      gmra->setSplitCriterium( *INTEGER(Rsplit) );
      gmra->dataFactory = getIKMDataFactory( *INTEGER(Rsimilarity), nCor, nSpatial );
      gmra->epsilon = *REAL(Reps);
      gmra->nKids = *INTEGER(RnKids);
      gmra->threshold = *REAL(Rt);
      gmra->maxIter = *INTEGER(RmaxIter);
      gmra->minPoints = *INTEGER(RminPoints);

      gmra->addPoints(pts);


      SEXP Rid;
      PROTECT(Rid = Rf_allocVector(INTSXP, 1));
      INTEGER(Rid)[0] = idCounter;
      SEXP ptr = R_MakeExternalPtr(gmra, Rf_install("GMRA"), R_NilValue);
      PROTECT(ptr);
      Rf_setAttrib(Rid, Rf_install("gmra_ptr"), ptr); 
      SET_VECTOR_ELT(list, i, Rid);

      UNPROTECT(2);
      ++idCounter;
    }

    UNPROTECT(1);

    return list;

  };





  //Get pointer to gmra
  GMRATree<double> *getGMRA(SEXP Rgmra){     
    SEXP Rgmra_pointer = Rf_getAttrib(Rgmra, Rf_install("gmra_ptr") );
    GMRATree<double> *gmra = static_cast<GMRATree<double> *>( R_ExternalPtrAddr(Rgmra_pointer) );
    return gmra;
  };


  SEXP createNodePointer( GMRANode<double> *node){
      if(node == NULL){
        return R_NilValue;
      }
      SEXP Rid;
      PROTECT(Rid = Rf_allocVector(INTSXP, 1));
      INTEGER(Rid)[0] = nodeCounter;
      ++nodeCounter;
      SEXP ptr = R_MakeExternalPtr(node, Rf_install("GMRANode"), R_NilValue);
      PROTECT(ptr);
      Rf_setAttrib(Rid, Rf_install("gmra_node_ptr"), ptr); 
      UNPROTECT(2);
      return Rid; 
  };


  GMRANode<double> *getGMRANode(SEXP Rnode){     
    SEXP Rnode_pointer = Rf_getAttrib(Rnode, Rf_install("gmra_node_ptr") );
    GMRANode<double> *node = static_cast<GMRANode<double> *>( R_ExternalPtrAddr(Rnode_pointer) );
    return node;
  };


  void getGMRANodes(SEXP Rnodes, std::vector< GMRANode<double> *> &nodes){
    for(int i=0; i < Rf_length(Rnodes); i++){
      SEXP Rnode = VECTOR_ELT(Rnodes, i);
      nodes.push_back( getGMRANode(Rnode) );
    }
  };




  //gmra based density estimation
  SEXP gmra_density_estimate(SEXP RX, SEXP Rm, SEXP Rn, SEXP Rgmras, SEXP
      RminPoints, SEXP RmaxPoints, SEXP Rsigma, SEXP RdType){ 

    NodeDistance<double> *dist = getNodeDistance(RdType);

    double *data = REAL(RX);
    int n= *INTEGER(Rn);
    int m= *INTEGER(Rm);

    MeanGMRADensityEstimator<double> density;
    int nPoints = 0;
    for(int i=0; i < Rf_length(Rgmras); i++){
      SEXP Rgmra = VECTOR_ELT(Rgmras, i);
      GMRATree<double> *gmra = getGMRA(Rgmra);

      gmra->computeRadii(dist);
      gmra->computeLocalRadii(dist);

      if(gmra != NULL){
        density.addTree(gmra);
        if(i == 0){
          nPoints = gmra->getRoot()->getPoints().size();
        }
        else{
          int tmp = gmra->getRoot()->getPoints().size();
          if(tmp != nPoints){
            return R_NilValue;
          }
        }
      }
    }


    double sigma = *REAL(Rsigma);
    int minP = *INTEGER(RminPoints);
    int maxP = *INTEGER(RmaxPoints);
    KernelNodeDensityEstimator<double> estimator(minP, maxP, sigma);
    //MinMaxNodeDensityEstimator<double> estimator(minP, maxP);
    //ScaleNodeDensityEstimator<double> estimator(minP, maxP);
    //MinMaxDensityEstimatorVisitor<double> estimator(nPoints, minP, maxP);
    //int scale = *INTEGER(Rscale);
    //PointsDensityEstimatorVisitor<double> estimator(nPoints, scale);    
    
    Eigen::MatrixXd X = Eigen::Map<Eigen::MatrixXd>(data,m,n);
    std::vector<double> p = density.estimate(X, estimator);


    SEXP Rp;
    PROTECT(Rp = Rf_allocVector(REALSXP, p.size()));
    memcpy(REAL(Rp), p.data(), p.size()*sizeof(double) );
    UNPROTECT(1);

    delete dist;

    return Rp;


  };






  //--- Mean signal predictor ----//
  SEXP gmra_mean_signal_estimate(SEXP Rsignal, SEXP Rsn, SEXP Rx, SEXP Rm, SEXP
      Rn, SEXP Rgmras, SEXP RminPoints){ 
 

    int m = *INTEGER(Rm);
    int n = *INTEGER(Rn);
    double *x = REAL(Rx);

    Eigen::MatrixXd X = Eigen::Map<Eigen::MatrixXd>(x,m,n);

    double *sig = REAL(Rsignal);
    int sn= *INTEGER(Rsn);

    std::vector< double > signal(sn, 0);
    for(int i=0; i<sn; i++){
      signal[i] = sig[i];
    };
    int minP = *INTEGER(RminPoints);

    SEXP list;
    PROTECT( list = Rf_allocVector(VECSXP, Rf_length(Rgmras) ) );
    
    for(int i=0; i < Rf_length(Rgmras); i++){
      SEXP Rgmra = VECTOR_ELT(Rgmras, i);
      GMRATree<double> *gmra = getGMRA(Rgmra);

      MeanGMRAPredictor<double> meanP(*gmra, signal, minP);
      std::vector<double> p= meanP.getPredictions(X);

      SEXP Rp;
      PROTECT(Rp = Rf_allocVector(REALSXP, p.size()));
      memcpy(REAL(Rp), p.data(), p.size()*sizeof(double) );
      SET_VECTOR_ELT(list, i, Rp);

    }


    UNPROTECT( Rf_length(Rgmras) + 1);
    return list;

  };






  //--- GMRA query methds ---//
  

  void collectNodes(GMRATree<double> *tree, std::vector<GMRANode<double> *>
      &collected, int scale = std::numeric_limits<int>::max()) {

    std::list<GMRANode<double>*> nodes;
    nodes.push_back( tree->getRoot() );

    while( !nodes.empty() ){
      GMRANode<double> *n = nodes.front();
      nodes.pop_front();

      std::vector<GMRANode<double> *> kids = n->getChildren();

      if(n->getScale() == scale || kids.size() == 0){
        collected.push_back(n);
      }
      else{
        //For each kid check if nearest neighbors within epsilon are possible.
        for(std::vector< GMRANode<double>* >::iterator it = kids.begin(); it != kids.end(); ++it){
          nodes.push_back(*it);
        }
      }
    }
  };



  SEXP getCenters(std::vector<GMRANode<double>*> &nodes){
    using namespace Eigen;

    if(nodes.size() == 0){
      return R_NilValue;
    }

    MatrixXd C( nodes.front()->getCenter().size(), nodes.size() );
    int index = 0;
    for(std::vector<GMRANode<double> *>::iterator it = nodes.begin(); it !=
        nodes.end(); ++it, ++index){
      C.col(index) = (*it)->getCenter();
    }

    SEXP RC;
    PROTECT( RC = Rf_allocMatrix(REALSXP, C.rows(), C.cols()) );
    memcpy( REAL(RC), C.data(), C.rows()*C.cols()*sizeof(double) );
    UNPROTECT(1);

    return RC;
  
  };
  


  SEXP getNumberOfPoints(std::vector<GMRANode<double>*> &nodes){
    using namespace Eigen;

    if(nodes.size() == 0){
      return R_NilValue;
    }

    VectorXi C( nodes.size() );
    int index = 0;
    for(std::vector<GMRANode<double> *>::iterator it = nodes.begin(); it !=
        nodes.end(); ++it, ++index){
      C(index) = (*it)->getPoints().size();
    }

    SEXP RC;
    PROTECT( RC = Rf_allocVector(INTSXP, C.size()) );
    memcpy( INTEGER(RC), C.data(), C.size()*sizeof(int) );
    UNPROTECT(1);

    return RC;
  
  };



  SEXP getRadii(std::vector<GMRANode<double>*> &nodes, NodeDistance<double> *dist){
    using namespace Eigen;

    if(nodes.size() == 0){
      return R_NilValue;
    }

    VectorXd C( nodes.size() );
    int index = 0;
    for(std::vector<GMRANode<double> *>::iterator it = nodes.begin(); it !=
        nodes.end(); ++it, ++index){
      (*it)->computeRadius(dist);
      C(index) = (*it)->getRadius( );
    }

    SEXP RC;
    PROTECT( RC = Rf_allocVector(REALSXP, C.size()) );
    memcpy( REAL(RC), C.data(), C.size()*sizeof(double) );
    UNPROTECT(1);

    return RC;
  };




  SEXP getPoints(std::vector<GMRANode<double>*> &nodes){
    using namespace Eigen;

    if(nodes.size() == 0){
      return R_NilValue;
    }

    SEXP list;
    PROTECT( list = Rf_allocVector( VECSXP, nodes.size() ) );

    int index = 0;
    for(std::vector<GMRANode<double> *>::iterator it = nodes.begin(); it !=
        nodes.end(); ++it, ++index){
      std::vector<int> &tmp = (*it)->getPoints();
      SEXP RC;
      PROTECT( RC = Rf_allocVector(INTSXP, tmp.size()) );
      memcpy( INTEGER(RC), tmp.data(), tmp.size()*sizeof(int) );
      SET_VECTOR_ELT(list, index, RC);
      UNPROTECT(1);
    }

    UNPROTECT(1);

    return list;
  
  };




  SEXP getNeighborhoodCenters(GMRANeighborhood<double>::NeighborList &nodes){
    using namespace Eigen;

    if(nodes.size() == 0){
      return R_NilValue;
    }

    MatrixXd C( nodes.front().second->getCenter().size(), nodes.size() );
    int index = 0;
    for(GMRANeighborhood<double>::NeighborListIterator it = nodes.begin(); it !=
        nodes.end(); ++it, ++index){
      C.col(index) = it->second->getCenter();
    }
    SEXP RC;
    PROTECT( RC = Rf_allocMatrix(REALSXP, C.rows(), C.cols()) );
    memcpy( REAL(RC), C.data(), C.rows()*C.cols()*sizeof(double) );
    UNPROTECT(1);

    return RC;
  };









  SEXP getNeighborhood(const Eigen::VectorXd &x, GMRATree<double> *tree, double eps, int scale, SEXP RdType){

    NodeDistance<double> *dist = getNodeDistance(RdType);
    QueryGMRANode<double> query(x);
    GenericGMRANeighborhood<double> nh(tree, dist);
    

    GMRANeighborhood<double>::NeighborList nodes;
    nh.neighbors(&query, eps, nodes, scale);
    
    delete dist;
    return getNeighborhoodCenters(nodes);

  };








  

  SEXP gmra_neighborhood(SEXP Rx, SEXP Rn, SEXP Rgmra, SEXP Reps, SEXP Rscale, SEXP RdType){
    using namespace Eigen;
    double eps =  *REAL(Reps);
    int scale = *INTEGER(Rscale);
    int n = *INTEGER(Rn);
    Map<VectorXd> x(REAL(Rx), n);

    GMRATree<double> *gmra = getGMRA(Rgmra);
    return getNeighborhood(x, gmra, eps, scale, RdType);
  };



  SEXP gmra_get_root(SEXP Rgmra){
    GMRATree<double> *gmra = getGMRA(Rgmra);
    return createNodePointer( gmra->getRoot() );
  };


  
  SEXP gmra_node_parent(SEXP Rnode){
    GMRANode<double> *node = getGMRANode(Rnode);
    GMRANode<double> *parent = node->getParent();
    return createNodePointer(parent);
  };



  SEXP gmra_node_children(SEXP Rnode){
    GMRANode<double> *node = getGMRANode(Rnode);
    std::vector< GMRANode<double>  *> &nodes = node->getChildren();
        
    SEXP list;
    PROTECT( list = Rf_allocVector(VECSXP, nodes.size() ) );
    for(int i=0; i < nodes.size(); i++){

      SET_VECTOR_ELT(list, i, createNodePointer(nodes[i]) );
    }
    UNPROTECT( 1);
    return list;
  };



  SEXP gmra_centers(SEXP Rgmra, SEXP Rscale){
    int scale = *INTEGER(Rscale);

    GMRATree<double> *gmra = getGMRA(Rgmra);
    std::vector<GMRANode<double> *> nodes;
    collectNodes(gmra, nodes, scale);
    return getCenters(nodes);
  };
  


  SEXP gmra_nodes(SEXP Rgmra, SEXP Rscale){
    int scale = *INTEGER(Rscale);

    GMRATree<double> *gmra = getGMRA(Rgmra);
    std::vector<GMRANode<double> *> nodes;
    collectNodes(gmra, nodes, scale);
        
    SEXP list;
    PROTECT( list = Rf_allocVector( VECSXP, nodes.size() )  );
    for(int i=0; i < nodes.size(); i++){
      SET_VECTOR_ELT(list, i, createNodePointer(nodes[i]) );
    }
    UNPROTECT( 1);
    return list;
  };



  SEXP gmra_node_centers(SEXP Rnodes){
    std::vector< GMRANode<double> *> nodes;
    getGMRANodes(Rnodes, nodes);
    return getCenters(nodes);    
  };



  SEXP gmra_radii(SEXP Rgmra, SEXP Rscale, SEXP RdType){
    int scale = *INTEGER(Rscale);
    NodeDistance<double> *dist = getNodeDistance(RdType);
    GMRATree<double> *gmra = getGMRA(Rgmra);
    std::vector< GMRANode<double> *> nodes;
    collectNodes(gmra, nodes, scale);
    return getRadii(nodes, dist);
  };



  SEXP gmra_node_radii(SEXP Rnodes, SEXP RdType){
    NodeDistance<double> *dist = getNodeDistance(RdType);
    std::vector< GMRANode<double> *> nodes;
    getGMRANodes(Rnodes, nodes);
    return getRadii(nodes, dist);    
  };



  SEXP gmra_npoints( SEXP Rgmra, SEXP Rscale){
    int scale = *INTEGER(Rscale);

    GMRATree<double> *gmra = getGMRA(Rgmra);
    std::vector<GMRANode<double> *> nodes;
    collectNodes(gmra, nodes, scale);
    return getNumberOfPoints(nodes);
  };



  SEXP gmra_node_npoints(SEXP Rnodes){
    std::vector< GMRANode<double> *> nodes;
    getGMRANodes(Rnodes, nodes);
    return getNumberOfPoints(nodes);    
  };



  SEXP gmra_partition(SEXP Rgmra, SEXP Rscale){
    int scale = *INTEGER(Rscale);

    GMRATree<double> *gmra = getGMRA(Rgmra);
    std::vector<GMRANode<double> *> nodes;
    collectNodes(gmra, nodes, scale);
    return getPoints(nodes);
  };



  SEXP gmra_node_points(SEXP Rnodes){
    std::vector< GMRANode<double> *> nodes;
    getGMRANodes(Rnodes, nodes);
    return getPoints(nodes);    
  };







  SEXP gmra_prune_min_points_scale(SEXP Rgmra, SEXP Rscale, SEXP Rn){
    GMRATree<double> *gmra = getGMRA(Rgmra);
    MinPointsPruneVisitor<double> cutter(*INTEGER(Rscale), *INTEGER(Rn) );
    gmra->depthFirstVisitor(&cutter);
    return R_NilValue;
  };
    

/*
std::vector< GWT<double>* > active;

SEXP gwtIteratedPCA(SEXP Rdata, SEXP Rm, SEXP Rn, SEXP Rd, SEXP Reps) {
  int m = *INTEGER(Rm);
  int n = *INTEGER(Rn);
  double *data = REAL(Rdata);

  double eps = *REAL(Reps);
  int d = *INTEGER(Rd);


  DenseMatrix<double> X(m, n, data);
 
  GWTTree<double> *gwtTree =  new IteratedPCATree<double>(d, eps);
  gwtTree->construct(X);

  GWT<double> *gwt = new ThinGWT<double>();
  gwt->setTree(gwtTree); 

  active.push_back(gwt);
  
  SEXP RgwtID;
  PROTECT(RgwtID = Rf_allocVector(INTSXP, 1));
  INTEGER(RgwtID)[0] = active.size()-1;
  UNPROTECT(1);
  
  return RgwtID;
};





SEXP gwtEncode(SEXP RgwtID, SEXP Rdata, SEXP Rm, SEXP Rn){
  int gwtID = *INTEGER(RgwtID);
  int m = *INTEGER(Rm);
  int n = *INTEGER(Rn);
  double *data = REAL(Rdata);

  GWT<double> *gwt = active[gwtID];

  DenseMatrix<double> X(m, n, data);

  DenseVector<double> x = Linalg<double>::ExtractColumn(X, 0);
  GWTCoefficents<double> c = gwt->encode(x);
  int cm = c.coeff.M(); 
  int cn = c.coeff.N();
 

  DenseMatrix<double>  C(cm, cn*n);
  Linalg<double>::SetColumns(C, 0, cn, c.coeff, 0);
  c.coeff.deallocate();

  for(int i=1; i<n; i++){
    Linalg<double>::ExtractColumn(X, 0, x);
    GWTCoefficents<double> c = gwt->encode(x);
    Linalg<double>::SetColumns(C, i*cn, (i+1)*cn, c.coeff, 0);
    c.coeff.deallocate();
  }

  x.deallocate();
  

  SEXP list;
  PROTECT( list = Rf_allocVector(VECSXP, 2) );

  SEXP RC;
  PROTECT(RC = Rf_allocMatrix(REALSXP, C.M(), C.N()));
  memcpy( REAL(RC), C.data(), C.M()*C.N()*sizeof(double) );
  SET_VECTOR_ELT(list, 0, RC);
  
  SEXP RCdim;
  PROTECT(RCdim = Rf_allocVector(INTSXP, 2));
  int *Cdim =INTEGER(RCdim);
  Cdim[0] = cm;
  Cdim[1] = cn;
  SET_VECTOR_ELT(list, 2, RCdim);
   
  UNPROTECT(3);

  C.deallocate();
 
  return list;
};


SEXP gwtProject(SEXP RgwtID, SEXP Rdata, SEXP Rm, SEXP Rn, SEXP Rscale){
  int gwtID = *INTEGER(RgwtID);
  int m = *INTEGER(Rm);
  int n = *INTEGER(Rn);
  double *data = REAL(Rdata);
  int scale= *INTEGER(Rscale);

  GWT<double> *gwt = active[gwtID];

  DenseMatrix<double> X(m, n, data);
  DenseMatrix<double> Xr(m, n);

  DenseVector<double> x(m);
  
  for(int i=0; i<n; i++){
    Linalg<double>::ExtractColumn(X, i, x);
    GWTCoefficents<double> c = gwt->encode(x);
    DenseVector<double> xr = gwt->decode(c, scale);
    c.coeff.deallocate();
    Linalg<double>::SetColumn(Xr, i, xr);
    xr.deallocate();
  }

  x.deallocate();
  
  SEXP RXr;
  PROTECT(RXr = Rf_allocMatrix(REALSXP, Xr.M(), Xr.N()));
  memcpy( REAL(RXr), Xr.data(), Xr.M()*Xr.N()*sizeof(double) );
  UNPROTECT(1); 

  Xr.deallocate();
  return RXr;
};

*/


}//end extern C
