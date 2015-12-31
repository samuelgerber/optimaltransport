#define VERBOSE


#include "Precision.h"

#include <stdio.h>

#include "EigenLinalgIO.h"
#include "GMRADensityEstimator.h"
#include "IPCATree.h"

#include <tclap/CmdLine.h>
#include "NodeDistance.h"

IPCANodeFactory<Precision> *getNodeFactory(GMRADataObject<Precision> *D, int d, double t){
  if(t < 0){
    return new FixedNodeFactory<Precision>(D, d);
  }
  else if(t < 1){
    return new RelativePrecisionNodeFactory<Precision>(D, d, t);
  }
  else{
    return new RelativeRatioNodeFactory<Precision>(D, d, t);
  }
};








int main(int argc, char **argv){
  
  //Command line parsing
  TCLAP::CmdLine cmd("Create Iterated PCA Tree denisty estimator. Create a tree by iteratively splitting a point cloud according to local principal components.", ' ', "1");

  TCLAP::ValueArg<int> dArg("d","dim","Maximal number of principal comonents to compute per node", false, 3, "integer");
  cmd.add(dArg);

  TCLAP::ValueArg<int> cArg("c","childDim","Maximal number of children per node. I.e. each node can have at most 2^c children", false, 3, "integer");
  cmd.add(cArg);

  TCLAP::ValueArg<Precision> tArg("t", "threshold", "Threshold to determine dimensionality of a node. For t is between 0 and 1 min(d, k) principal components ar retained where k is the smallest number of principal components to capture 100*t percent of the variance of the data points in the node. For negative t each node has 2^d children if possible. For t > 1 it is min(d, k) the ratio between the variance of the k'th over the k'th+1 principal components has to be less than t.", false, -1, "double");
  cmd.add(tArg);

  TCLAP::ValueArg<Precision> eArg("e","epsilon","Parameter for stopping criterion, stop iterative splitting procedure when the specified stopping criterion is less than epsilon", false, 0, "double");
  cmd.add(eArg);
  
  TCLAP::ValueArg<int> stopArg("","stop","Stopping criterion. 1: Local R^2, stop if the R^2 of a node, i.e. R^2 with respect to the variance in the node is less than epsilon. 2: Total R^2, stop if the total R^2, i.e. R^2 with respect to the global varianc of the data is less than epsilon. 3: Node MSE, stop if the mean squared error of a node is less than epsilon 4: Node Radius, stop if the radius of a node is less than epsilon, 5: Relative node radius, stop if the radius of the node over the global radius is less than epsilon", false, 2, "integer");
  cmd.add(stopArg);
  
  TCLAP::ValueArg<int> splitArg("","split","Splitting strategy determines how to split a node into its children along the local split direction. 1: Mean, split nodes according to the positing of points with respect to the mean of the node. 2: Midpoint, split across the midpoint of each principal component. 3: Random Mean, add random shift to the mean. 4: Random Midpoint, add random shift to the midpoint location.", false, 3, "integer");
  cmd.add(splitArg);

    TCLAP::ValueArg<int> splitDirArg("","splitDir","Splitting direction selects local split direction. 1: Principal components 2: Randomized principal component.", false, 2, "integer");
  cmd.add(splitDirArg);

  TCLAP::ValueArg<int> nArg("n","nRuns","Number of trees to build", true, 2, "integer");
  cmd.add(nArg);

  TCLAP::ValueArg<int> sArg("s","scale","Scale at which to estimate density", false, std::numeric_limits<int>::max(), "integer");
  cmd.add(sArg);

  TCLAP::ValueArg<std::string> dataArg("x","data", "File with matrix for data point cloud, with each column a data point",  true, "", "matrix header file");
  cmd.add(dataArg);


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

  MatrixXp X = EigenLinalg::LinalgIO<Precision>::readMatrix(dataArg.getValue());
 
  MatrixGMRADataObject<Precision> data(X);


  
  IPCANodeFactory<Precision> *factory = getNodeFactory( &data, dArg.getValue(), tArg.getValue() );
  factory->epsilon = eArg.getValue();


  factory->stop = IPCANodeFactory<Precision>::TOTAL_R2;
  switch( stopArg.getValue() ){
    case 1:
      factory->stop = IPCANodeFactory<Precision>::NODE_R2;
      break;
    case 2:
      factory->stop = IPCANodeFactory<Precision>::TOTAL_R2;
      break;
    case 3:
      factory->stop = IPCANodeFactory<Precision>::NODE_MSE;
      break;
    case 4:
      factory->stop = IPCANodeFactory<Precision>::NODE_RADIUS;
      break;
    case 5:
      factory->stop = IPCANodeFactory<Precision>::RELATIVE_NODE_RADIUS;
      break;    
    case 6:
      factory->stop = IPCANodeFactory<Precision>::MASS_RADIUS;
      break;
  }



  factory->splitStrategy = IPCANodeFactory<Precision>::MEAN;
  switch( splitArg.getValue() ){
    case 1:
       factory->splitStrategy = IPCANodeFactory<Precision>::MEAN;
      break;
    case 2:
       factory->splitStrategy = IPCANodeFactory<Precision>::MIDPOINT;
      break;
    case 3:
       factory->splitStrategy = IPCANodeFactory<Precision>::RANDOM_MEAN;
      break;
    case 4:
       factory->splitStrategy = IPCANodeFactory<Precision>::RANDOM_MIDPOINT;
      break;
  }

  factory->maxKidDim = cArg.getValue();
    
  int scale = sArg.getValue();

  std::vector<int> pts(X.cols() );
  for(int i=0; i<pts.size(); i++){
    pts[i] = i;
  };


  MaxGMRADensityEstimator<Precision> density;
  int nRuns = nArg.getValue();
  for(int i=0; i<nRuns; i++){
    IPCATree<Precision> *tree = new IPCATree<Precision>(&data, factory);
    tree->addPoints(pts); 
    density.addTree(tree);
  }
      

  //NodeDensityEstimator<Precision> *estimator = new MinMaxNodeDensityEstimator<Precision>(5, 20);
  //NodeDensityEstimator<Precision> *estimator = new ScaleNodeDensityEstimator<Precision>(5, 20);
  NodeDensityEstimator<Precision> *estimator = new KernelNodeDensityEstimator<Precision>(20, 100, 0.5);
  std::vector<Precision> p = density.estimate(X, *estimator);

  std::cout << p[0] << std::endl;



  return 0;
}
catch(...){
  std::cout << "Not a matrix header file" << std::endl;
  return -1;
}

}
