#include "Precision.h"

#include <stdio.h>

#include "Wasserstein.h"
#include "LinalgIO.h"
#include "IPCATree.h"
#include "DenseMatrix.h"
#include <tclap/CmdLine.h>

int main(int argc, char **argv){
  
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
  
  std::list< GMRANode * > nodes;
  nodes.push_back( ipcaTree.getRoot() );

  std::vector< DenseMatrix<Precision> > phi;
  std::vector< DenseVector<Precision> > sigma;

  while( !nodes.empty() ){
    IPCANode<Precision> *node = (IPCANode<Precision> *) nodes.front();
    nodes.pop_front();


    phi.push_back(node->phi);
    DenseVector<Precision> sigma2 = Linalg<Precision>::Copy(node->sigma);
    for(int i=0; i<sigma2.N(); i++){
      sigma2(i) *= sigma2(i);
    }
    sigma.push_back(sigma2);
    std::list< GMRANode* > children = node->getChildren();
    for(std::list< GMRANode* >::iterator it = children.begin(); it !=
        children.end(); ++it){
      nodes.push_back(*it); 
    }
  }



  DenseMatrix<Precision> D(phi.size(), phi.size());
  for(int i=0; i < (phi.size()-1); i++){
    D(i, i) = 0;
    for(int j=i+1; j<phi.size(); j++){
      D(i, j) = Wasserstein<Precision>::distance2(phi[i], sigma[i], phi[j], sigma[j]); 
      D(j, i) = D(i, j);
    }
  }
  D(D.N()-1, D.N()-1) = 0;


  LinalgIO<Precision>::writeMatrix(outFileName, D);

  return 0;

}
