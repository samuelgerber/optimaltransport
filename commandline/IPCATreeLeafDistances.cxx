#include "Precision.h"

#include <stdio.h>

#include "Wasserstein.h"
#include "LinalgIO.h"
#include "IPCATree.h"
#include "EuclideanMetric.h"
#include "DenseMatrix.h"
#include "CmdLine.h"



int findParent(int id, int np, std::vector<int> &parent){
    for(int k=0; k < np && id != -1; k++){
      id = parent[id];
    } 
    return id;
};


int main(int argc, char **argv){
  
  //Command line parsing
  TCLAP::CmdLine cmd("GWT", ' ', "1");

  TCLAP::ValueArg<std::string> tArg("t","tree", "Tree data file",  true, "",
      "tree data file");
  cmd.add(tArg);

  TCLAP::ValueArg<int> pArg("p","parents", "How many levels up to see if points are connected",  false, 0,
      "integer");
  cmd.add(pArg);


  TCLAP::ValueArg<std::string> outArg("o","out", "Data file prefirx",  true, "", "file to store results in");
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
  int np = pArg.getValue();

  IPCATree<Precision> ipcaTree;

  std::ifstream treeFile(treeFileName.c_str(), std::ios_base::in | std::ios_base::binary);
  ipcaTree.unflatten(treeFile);
  
  treeFile.close();

  std::cout << "read tree " << std::endl;
  
  std::list< GMRANode * > nodes;
  nodes.push_back( ipcaTree.getRoot() );
  
  std::list< int > scales;
  scales.push_back( 0 );

  std::list< int > pids;
  pids.push_back( -1 );
  
  std::vector< DenseMatrix<Precision> > phi;
  std::vector< DenseVector<Precision> > sigma;
  std::vector< DenseVector<Precision> > centers;
  std::vector< Precision > radii;
  std::vector< int > scale;
  std::vector< int > ids;
  std::vector< int > parent;

  int nodeID = 0;
  while( !nodes.empty() ){
    IPCANode<Precision> *node = (IPCANode<Precision> *) nodes.front();
    nodes.pop_front();

    int s = scales.front();
    scales.pop_front();

    int pid = pids.front();
    pids.pop_front();
  
    parent.push_back(pid);
    
    std::list< GMRANode* > children = node->getChildren(); 
    if(children.empty()){
      scale.push_back(s);

      radii.push_back( node->r );

      ids.push_back(nodeID);

      phi.push_back(node->phi);
      sigma.push_back(node->sigma);
      centers.push_back(node->center);
    }
      
    
    for(std::list< GMRANode* >::iterator it = children.begin(); it !=
        children.end(); ++it){
      pids.push_back(nodeID);
      nodes.push_back(*it);
      scales.push_back(s+1); 
    }
    nodeID++;
  }



  DenseMatrix< Precision>  D( centers.size(), centers.size() );
  EuclideanMetric< Precision > l2;
  for(int i=0; i < (phi.size()-1); i++){
    D(i, i) = 0;

    int pid1 = findParent(ids[i], np, parent); 

    for(int j=i+1; j<phi.size(); j++){
      int pid2 = findParent(ids[j], np, parent); 
      if(pid2 == pid1){
        D(i, j) =  l2.distance(centers[i], centers[j]);
      }
      else{
        D(i, j) = -1;
      }
      D(j, i) = D(i, j);
    }
  }
  D( D.N()-1, D.N()-1 ) = 0;

  LinalgIO<Precision>::writeMatrix(outFileName, D);

  return 0;

}


