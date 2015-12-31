#define VERBOSE


#include "Precision.h"

#include <stdio.h>

#include "EuclideanMetric.h"
#include "LinalgIO.h"
#include "AxisAlignedTree.h"
#include "DenseMatrix.h"
#include <tclap/CmdLine.h>
#include "NodeDistance.h"
  







int main(int argc, char **argv){
  
  //Command line parsing
  TCLAP::CmdLine cmd("Create Axis Aligned Tree. Create a tree by iteratively splitting a point cloud according to local principal components.", ' ', "1");


  TCLAP::ValueArg<Precision> eArg("e","epsilon","Parameter for stopping criterion, stop iterative splitting procedure when the specified stopping criterion is less than epsilon", true, 0.1, "double");
  cmd.add(eArg);
  
  TCLAP::ValueArg<int> stopArg("","stop"," 1: Node Radius, stop if the radius of a node is less than epsilon, 2: Relative node radius, stop if the radius of the node over the global radius is less than epsilon, 3: Mass * Raiuds of node less than epsilon", true, 2, "integer");
  cmd.add(stopArg);
  
  TCLAP::ValueArg<int> splitArg("","split","Splitting strategy determines how to split a node into its children along the local principal components. 1: Mean, split nodes according to the positing of points with respect to the mean of the node. 2: Midpoint, split across the midpoint of each principal component. 3: Random Mean, add random shift to the mean. 4: Random Midpoint, add random shift to the midpoint location.", true, 2, "integer");
  cmd.add(splitArg);

  TCLAP::ValueArg<std::string> dataArg("x","data", "File with matrix for data point cloud, with each column a data point",  true, "", "matrix header file");
  cmd.add(dataArg);

  TCLAP::ValueArg<std::string> outArg("o","out", "File to store IPCATree in binary format in",  true, "", "file name");
  cmd.add(outArg);

  try{
	  cmd.parse( argc, argv );
	} 
  catch (TCLAP::ArgException &e){ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    return -1;
  }

try{

  FortranLinalg::DenseMatrix<Precision> X = FortranLinalg::LinalgIO<Precision>::readMatrix(dataArg.getValue());
  Precision eps = eArg.getValue();
  std::string fileName = outArg.getValue();




  AxisAlignedTree<Precision>::StoppingCriterium stop = AxisAlignedTree<Precision>::RELATIVE_NODE_RADIUS;
  switch( stopArg.getValue() ){
    case 1:
      stop = AxisAlignedTree<Precision>::NODE_RADIUS;
      break;
    case 2:
      stop = AxisAlignedTree<Precision>::RELATIVE_NODE_RADIUS;
      break;    
    case 3:
      stop = AxisAlignedTree<Precision>::MASS_RADIUS;
      break;
  }



  AANode<Precision>::SplitStrategy split = AANode<Precision>::MEAN;
  switch( splitArg.getValue() ){
    case 1:
      split = AANode<Precision>::MEAN;
      break;
    case 2:
      split = AANode<Precision>::MIDPOINT;
      break;
    case 3:
      split = AANode<Precision>::RANDOM_MEAN;
      break;
    case 4:
      split = AANode<Precision>::RANDOM_MIDPOINT;
      break;
  }

  CenterNodeDistance<Precision> dist(new EuclideanMetric<Precision>());

  AxisAlignedTree<Precision> ipcaTree(eps, stop, split, X.M() );
  ipcaTree.construct(X.data(), X.M(), X.N() );

  std::ofstream file(fileName.c_str(), std::ios_base::out | std::ios_base::binary);
  ipcaTree.flatten(file);

  file.close();

  X.deallocate();


  return 0;
}
catch(...){
  std::cout << "Not a matrix header file" << std::endl;
  return -1;
}

}
