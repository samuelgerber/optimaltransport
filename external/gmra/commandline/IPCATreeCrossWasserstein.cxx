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
  
  TCLAP::ValueArg<double> pArg("p","percent", "Maximum percentage difference of raddi in order to compute Wasserstein distance between nodes",  true, 0.1,"file to store results in"); 
  cmd.add(pArg);
 
  try{
	  cmd.parse( argc, argv );
	} 
  catch (TCLAP::ArgException &e){ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    return -1;
  }

  
  std::string outFileName = outArg.getValue();
  std::string treeFileName = tArg.getValue();
  Precision p = pArg.getValue();

  IPCATree<Precision> ipcaTree;

  std::ifstream treeFile(treeFileName.c_str(), std::ios_base::in | std::ios_base::binary);
  ipcaTree.unflatten(treeFile);
  
  treeFile.close();

  std::cout << "read tree " << std::endl;
  
  std::list< GMRANode * > nodes;
  nodes.push_back( ipcaTree.getRoot() );
  
  std::list< int > scales;
  scales.push_back( 0 );

  std::vector< DenseMatrix<Precision> > phi;
  std::vector< DenseVector<Precision> > sigma;
  std::vector< Precision > radii;
  std::vector< int > scale;
  std::vector< int > ids;

  int nodeID = 0;
  while( !nodes.empty() ){
    IPCANode<Precision> *node = (IPCANode<Precision> *) nodes.front();
    nodes.pop_front();

    int s = scales.front();
    scales.pop_front();
    
    scale.push_back(s);

    radii.push_back( node->r );

    ids.push_back(nodeID);
    nodeID++;

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
      scales.push_back(s+1); 
    }
  }



  std::list<DenseVector<Precision> > dists;
  for(int i=0; i < (phi.size()-1); i++){
    for(int j=i+1; j<phi.size(); j++){
      Precision r1 = radii[i];
      Precision r2 = radii[j];
      Precision rd = p*std::max(r1, r2);
      
      if( fabs( r1 - r2 ) < rd ){
        
        DenseVector<Precision> v(7);
        v(0) = ids[i];
        v(1) = ids[j];
        v(2) = radii[i];
        v(3) = radii[j];

        v(4) = scale[i];
        v(5) = scale[j];
      
        v(6) =Wasserstein<Precision>::distance2(phi[i], sigma[i], phi[j], sigma[j]); 
        dists.push_back(v);
      }
    }
  }

  DenseMatrix<Precision> D = Linalg<Precision>::ToMatrix(dists);

  LinalgIO<Precision>::writeMatrix(outFileName, D);

  return 0;

}
