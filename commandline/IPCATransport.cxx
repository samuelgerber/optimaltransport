#include "Precision.h"

#include <stdio.h>

#include "EigenEuclideanMetric.h"
#include "MultiscaleSinkhornTransport.h"

#include "CPLEXNetworkSolver.h"
#include "CPLEXSolver.h"

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
#include "NodeDistance.h"
#include "WassersteinNodeDistance.h"
#include "RelativeGMRANeighborhood.h"
#include "IPCAGWT.h"

#include "EigenLinalgIO.h"

#include <tclap/CmdLine.h>

int main(int argc, char **argv){

  //Command line parsing
  TCLAP::CmdLine cmd("IPCA Transport", ' ', "1");

  TCLAP::ValueArg<std::string> x1Arg("","X1", "Point set 1", true, "",
      "matrix header file");
  cmd.add(x1Arg);

  TCLAP::ValueArg<std::string> x2Arg("","X2", "Point set 2", true, "",
      "matrix header file");
  cmd.add(x2Arg);

  TCLAP::ValueArg<std::string> outArg("o","out", "data file name", true, "", "file to store results in"); 
  cmd.add(outArg);
  
  TCLAP::ValueArg<int> s1Arg("","s1", "number of scales of tree 1", false, -1, "integer"); 
  cmd.add(s1Arg);
  
  TCLAP::ValueArg<int> s2Arg("","s2", "number of scales of tree 2", false, -1, "integer"); 
  cmd.add(s2Arg);

  TCLAP::ValueArg<Precision> rArg("r","rFactor", "Multiplicative factor for size of neighborhoods to consider when moving to finer scale", false, 0, "real"); 
  cmd.add(rArg);

  TCLAP::ValueArg<int> pArg("p","nRandom", "N randomized trials", false, 10, "integer"); 
  cmd.add(pArg);

  TCLAP::ValueArg<Precision> lArg("l","lambda", "Regularization factor", false, 0, "real"); 
  cmd.add(lArg);

  TCLAP::ValueArg<int> dArg("d","Wd", "W_d distance computation", false, 1, "integer"); 
  cmd.add(dArg);

/*
  TCLAP::ValueArg<int> pArg("p","propagation", "Solution proagation across scale", false, 0, "0 or 1"); 
  cmd.add(pArg);
  */


  TCLAP::ValueArg<int> sArg("s","strategy", "Neighborhood expansion strategy 0 =  Potential, 1 = Refine, 2 = Expand", false, 0, "0, 1 or 2"); 
  cmd.add(sArg);

  
  TCLAP::ValueArg<int> nArg("n","neighborhood", "Neighborhood type 0 = Generic, 1 = Relative", false, 0, "0 or 1"); 
  cmd.add(nArg);



 // TCLAP::SwitchArg mcfArg("m","mcf", "Use network simplex", false, 0); 
 // cmd.add(mcfArg);

  try{ 
    cmd.parse( argc, argv );
  } 
  catch (TCLAP::ArgException &e){ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    return -1;
  }



  std::string outFileName = outArg.getValue();
  int nScales1 = s1Arg.getValue();
  int nScales2 = s2Arg.getValue();
  int p = dArg.getValue();
  int sType = sArg.getValue();
  
  typedef GMRANode<Precision>::MatrixXp MatrixXp;
  
  MatrixXp X1 = EigenLinalg::LinalgIO<Precision>::readMatrix(x1Arg.getValue());
  MatrixXp X2 = EigenLinalg::LinalgIO<Precision>::readMatrix(x2Arg.getValue());
  
  MatrixGMRADataObject<Precision> data1(X1);
  MatrixGMRADataObject<Precision> data2(X2);
  
  IPCANodeFactory<Precision> *factory1 = new FixedNodeFactory<Precision>(&data1, 2);
  IPCANodeFactory<Precision> *factory2 = new FixedNodeFactory<Precision>(&data2, 2);

  IPCATree<Precision> *t1 = new IPCATree<Precision>(&data1, factory1);  
  std::vector<int> pts1(X1.cols() );
  for(int i=0; i<pts1.size(); i++){
    pts1[i] = i;
  };
  t1->addPoints(pts1);

  IPCATree<Precision> *t2 = new IPCATree<Precision>(&data2, factory2);
  std::vector<int> pts2(X2.cols() );
  for(int i=0; i<pts2.size(); i++){
    pts2[i] = i;
  };
  t2->addPoints(pts2);



  EuclideanMetric<Precision> *metric = new EuclideanMetric<Precision>();
  NodeDistance<Precision> *d = new CenterNodeDistance<Precision>( metric);
  NodeDistance<Precision> *d2 =  new WassersteinNodeDistance<Precision>();
  
  t1->computeRadii(d);
  t1->computeLocalRadii(d);
  
  t2->computeRadii(d);
  t2->computeLocalRadii(d);

  int nType = nArg.getValue();
  GMRANeighborhood<Precision> *nh1;
  GMRANeighborhood<Precision> *nh2;
  if(nType == 0){
   nh1 = new GenericGMRANeighborhood<Precision>(t1, d);
   nh2 = new GenericGMRANeighborhood<Precision>(t2, d);
  }
  else{
    nh1 = new RelativeGMRANeighborhood<Precision>(t1, d);
    nh2 = new RelativeGMRANeighborhood<Precision>(t2, d);
  }

  std::vector<double> weights;

  std::vector< MultiscaleTransportLevel<Precision> * > t1Levels =
    GMRAMultiscaleTransportLevel<Precision>::buildTransportLevels(*nh1, weights,
        true);

  std::vector< MultiscaleTransportLevel<Precision> * > t2Levels =
    GMRAMultiscaleTransportLevel<Precision>::buildTransportLevels(*nh2, weights,
        true);
  
  LPSolver *lpSolver = new CPLEXNetworkSolver( lArg.getValue() );
  //LPSolver *lpSolver = new CPLEXSolver(CPX_ALG_AUTOMATIC, lArg.getValue());
  MultiscaleTransportLP<Precision> transport(lpSolver);
 
  ExpandNeighborhoodStrategy<double> expand(2, 0, 1);
  RefineNeighborhoodStrategy<double> refine(2, 0, 1);
  PotentialNeighborhoodStrategy<double> potential(0, 0, false, false);
  //transport.addNeighborhodStrategy(&expand);

  NeighborhoodPropagationStrategy<double> prop1(0);
  transport.setPropagationStrategy1( &prop1);
  
  CapacityPropagationStrategy<double> cprop(3, 0);
  IteratedCapacityPropagationStrategy<double> icprop(3, 0);
  RandomizedNeighborhoodPropagationStrategy<double> rprop1(10);

  SinkhornPropagationStrategy<double> sprop;
  transport.setPropagationStrategy2( &sprop);
  MaxEntropyPropagationStrategy<double> meprop;
  transport.setPropagationStrategy2( &meprop);
     
  std::vector< TransportPlan<Precision> * > sols = transport.solve( t1Levels,
      t2Levels, p, nScales1, nScales2);

  

  for(unsigned int i=0; i < t1Levels.size(); ++i){
    delete t1Levels[i];
  }
  
  for(unsigned int i=0; i < t2Levels.size(); ++i){
    delete t2Levels[i];
  }


  t1Levels.clear();
  t2Levels.clear();

  for(int i=0; i<sols.size(); i++){
    delete sols[i];
  }
  sols.clear();
  sols = transport.solve( t1Levels,
      t2Levels, p, nScales1, nScales2);

 for(unsigned int i=0; i < t1Levels.size(); ++i){
    delete t1Levels[i];
  }
  
  for(unsigned int i=0; i < t2Levels.size(); ++i){
    delete t2Levels[i];
  }


  t1Levels.clear();
  t2Levels.clear();

  for(int i=0; i<sols.size(); i++){
    delete sols[i];
  }
  sols.clear();


  delete d;
  delete metric;
  delete nh1;
  delete nh2;
  delete t1;
  delete t2;
  delete factory1;
  delete factory2;
  //delete lpSolver;

  return 0;

}
