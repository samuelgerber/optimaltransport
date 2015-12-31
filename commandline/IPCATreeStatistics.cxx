#include "Precision.h"

#include <stdio.h>

#include "LinalgIO.h"
#include "IPCATree.h"
#include "DenseMatrix.h"
#include "CmdLine.h"

int main(int argc, char **argv){
  using namespace FortranLinalg;

  //Command line parsing
  TCLAP::CmdLine cmd("GWT", ' ', "1");

  TCLAP::ValueArg<std::string> tArg("t","tree", "Tree data file",  true, "",
      "tree data file");
  cmd.add(tArg);

  TCLAP::ValueArg<std::string> outArg("o","out", "Data file",  true, "", "file to store results in");
  cmd.add(outArg);

  try{
	  cmd.parse( argc, argv );
	} 
  catch (TCLAP::ArgException &e){ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    return -1;
  }

  
  std::string outFileName = outArg.getValue();
  std::string treeFileName = tArg.getValue();

  IPCATree<Precision> ipcaTree;

  std::ifstream treeFile(treeFileName.c_str(), std::ios_base::in | std::ios_base::binary);
  ipcaTree.unflatten(treeFile);
  treeFile.close();

  std::cout << "read tree " << std::endl;
  
  std::list< GMRANode<Precision> * > nodes;
  nodes.push_back( ipcaTree.getRoot() );

  std::list< int > depths;
  depths.push_back( 0 );
 
  std::list< DenseVector<Precision> > mse;
  std::list< int > scale;
  std::list< int > pID;
  std::list< Precision > rads;
  pID.push_back(-1);


  int ID = 0;
  while( !nodes.empty() ){
    IPCANode<Precision> *node = (IPCANode<Precision> *) nodes.front();
    nodes.pop_front();

    int depth = depths.front();   
    depths.pop_front();

    scale.push_back(depth);
    mse.push_back(node->mse);
    rads.push_back(node->getRadius() );

    std::vector< GMRANode<Precision>* > children = node->getChildren();
    for(std::vector< GMRANode<Precision>* >::iterator it = children.begin(); it !=
        children.end(); ++it){
      pID.push_back(ID);
      depths.push_back(depth+1);
      nodes.push_back(*it); 
    }
    ++ID;
  }


  DenseMatrix<Precision> A = Linalg<Precision>::ToMatrix(mse);
  std::stringstream ss1;
  ss1 << outFileName << "-mse.data";
  LinalgIO<Precision>::writeMatrix(ss1.str(), A);

  DenseVector<Precision> rs = Linalg<Precision>::ToVector(rads);
  std::stringstream ss2;
  ss2 << outFileName << "-radii.data";
  LinalgIO<Precision>::writeVector(ss2.str(), rs);

  DenseVector<int> depth = Linalg<int>::ToVector(scale);
  std::stringstream ss3;
  ss3 << outFileName << "-scale.data";
  LinalgIO<int>::writeVector(ss3.str(), depth);

  DenseVector<int> parent = Linalg<int>::ToVector(pID);
  std::stringstream ss4;
  ss4 << outFileName << "-pids.data";
  LinalgIO<int>::writeVector(ss4.str(), parent);

  return 0;

}
