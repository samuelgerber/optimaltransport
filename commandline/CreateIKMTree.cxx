#define VERBOSE


#include "Precision.h"

#include <stdio.h>

#include "CmdLine.h"

#include "EigenEuclideanMetric.h"
#include "EigenLinalgIO.h"
#include "IKMTree.h"
#include "GMRAPrune.h"
#include "NodeDistance.h"

#include "MeanGMRAPredictor.h"





int main(int argc, char **argv){
  
  //Command line parsing
  TCLAP::CmdLine cmd("Create Iterated PCA Tree. Create a tree by iteratively splitting a point cloud according to local principal components.", ' ', "1");

  TCLAP::ValueArg<int> dArg("d","dim","Maximal number of principal comonents to compute per node", true, 3, "integer");
  cmd.add(dArg);

  TCLAP::ValueArg<int> cArg("c","childDim","Maximal number of children per node. I.e. each node can have at most 2^c children", true, 3, "integer");
  cmd.add(cArg);

  TCLAP::ValueArg<Precision> tArg("t", "threshold", "Threshold to determine dimensionality of a node. For t is between 0 and 1 min(d, k) principal components ar retained where k is the smallest number of principal components to capture 100*t percent of the variance of the data points in the node. For negative t each node has 2^d children if possible. For t > 1 it is min(d, k) the ratio between the variance of the k'th over the k'th+1 principal components has to be less than t.", true, -1, "double");
  cmd.add(tArg);

  TCLAP::ValueArg<Precision> eArg("e","epsilon","Parameter for stopping criterion, stop iterative splitting procedure when the specified stopping criterion is less than epsilon", true, 0.1, "double");
  cmd.add(eArg);
  
  TCLAP::ValueArg<int> stopArg("","stop","Stopping criterion. 1: Local R^2, stop if the R^2 of a node, i.e. R^2 with respect to the variance in the node is less than epsilon. 2: Total R^2, stop if the total R^2, i.e. R^2 with respect to the global varianc of the data is less than epsilon. 3: Node MSE, stop if the mean squared error of a node is less than epsilon 4: Node Radius, stop if the radius of a node is less than epsilon, 5: Relative node radius, stop if the radius of the node over the global radius is less than epsilon", true, 2, "integer");
  cmd.add(stopArg);
  

  TCLAP::ValueArg<std::string> dataArg("x","data", "File with matrix for data point cloud, with each column a data point",  true, "", "matrix header file");
  cmd.add(dataArg);

  TCLAP::ValueArg<std::string> outArg("o","out", "File to store IPCATreeConfiguration in binary format in",  true, "", "file name");
  cmd.add(outArg);

  try{
	  cmd.parse( argc, argv );
	} 
  catch (TCLAP::ArgException &e){ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    return -1;
  }

try{

  typedef GMRANode<Precision>::MatrixXp MatrixXp;
  typedef GMRANode<Precision>::VectorXp VectorXp;

  std::string fileName = outArg.getValue();
  MatrixXp X = EigenLinalg::LinalgIO<Precision>::readMatrix(dataArg.getValue());

  MatrixGMRADataObject<Precision> data(X);
  IKMTree<Precision> ikmTree(&data);
  ikmTree.dataFactory = new L2GMRAKmeansDataFactory<Precision>(); 
  ikmTree.dataFactory = new JointCorrelationAndSpatialGMRAKmeansDataFactory<Precision>(X.rows()-3, 3); 
  ikmTree.nKids=16; 

  std::vector<int> pts(X.cols() );
  for(int i=0; i<pts.size(); i++){
    pts[i] = i;
  };

  ikmTree.addPoints( pts );


  std::vector< double > signal(pts.size(), 0);
  MeanGMRAPredictor<Precision> pred(ikmTree, signal, 10);  
  //MinPointsPruneVisitor<Precision> cutter(4, 10);
  //ikmTree.breadthFirstVisitor(&cutter);

  return 0;
}
catch(...){
  std::cout << "Not a matrix header file" << std::endl;
  return -1;
}

}
