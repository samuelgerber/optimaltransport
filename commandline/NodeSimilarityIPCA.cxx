#include "Precision.h"

#include <stdio.h>


#include "LinalgIO.h"
#include "IPCAGWT.h"
#include "RelativeSubspaceSimilarity.h"
#include "WassersteinSimilarity.h"
#include "RigidPatchSimilarity.h"
#include "GrassmannianNodeSimilarity.h"
#include "DenseMatrix.h"
#include "CmdLine.h"


int main(int argc, char **argv){
  using namespace FortranLinalg;

  //Command line parsing
  TCLAP::CmdLine cmd("IPCA GWT", ' ', "1");

  TCLAP::ValueArg<std::string> tArg("t","tree", "Tree data file", true, "",
      "tree data file");
  cmd.add(tArg);

  TCLAP::ValueArg<int> mArg("m","metric", "Metric type, 1=Wasserstein, 2=Relative subspace", true, 1,
      "integer");
  cmd.add(mArg);

  TCLAP::ValueArg<std::string> outArg("o","out", "data file name", true, "", "file to store results in"); 
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
  int metric = mArg.getValue();

  SimilarityVisitor<Precision> *dist;
  if(metric == 1){
    dist = new WassersteinSimilarity<Precision>();
  }
  else if(metric == 2){
    dist = new RelativeSubspaceSimilarity<Precision>();
  }
  else if(metric == 3){
    dist = new RigidPatchSimilarity<Precision>();
  }
  else if(metric == 4){
    dist = new GrassmannianNodeSimilarity<Precision>();
  }
  else{
    std::cerr << "No valid metric specified" << std::endl;
    exit(1);
  }

  IPCATree<Precision> ipcaTree;

  std::ifstream treeFile(treeFileName.c_str(), std::ios_base::in | std::ios_base::binary);
  ipcaTree.unflatten(treeFile);

  treeFile.close();

  std::cout << "read tree " << std::endl;


  IPCAGWT<Precision> gwt;
  gwt.setTree(&ipcaTree);

   

  ipcaTree.breadthFirstVisitor(dist);

  DenseMatrix<Precision> D = dist->distances();
   
  LinalgIO<Precision>::writeMatrix(outFileName, D);

  
  return 0;

}
