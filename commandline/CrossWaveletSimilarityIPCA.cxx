#include "Precision.h"

#include <stdio.h>


#include "LinalgIO.h"
#include "IPCAGWT.h"
#include "WaveletSimilarity.h"
#include "DenseMatrix.h"
#include <tclap/CmdLine.h>


int main(int argc, char **argv){

  //Command line parsing
  TCLAP::CmdLine cmd("IPCA GWT", ' ', "1");

  TCLAP::ValueArg<std::string> t1Arg("","t1", "Tree data file", true, "",
      "tree data file");
  cmd.add(t1Arg);

  TCLAP::ValueArg<std::string> t2Arg("","t2", "Tree data file", true, "",
      "tree data file");
  cmd.add(t2Arg);
  
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
  std::string tree1FileName = t1Arg.getValue();
  std::string tree2FileName = t2Arg.getValue();

  IPCATree<Precision> ipcaTree1;
  std::ifstream treeFile1(tree1FileName.c_str(), std::ios_base::in | std::ios_base::binary);
  ipcaTree1.unflatten(treeFile1);
  treeFile1.close();

  IPCATree<Precision> ipcaTree2;
  std::ifstream treeFile2(tree2FileName.c_str(), std::ios_base::in | std::ios_base::binary);
  ipcaTree2.unflatten(treeFile2);
  treeFile2.close();
  std::cout << "read tree " << std::endl;


  IPCAGWT<Precision> gwt1;
  gwt1.setTree(&ipcaTree1);

  IPCAGWT<Precision> gwt2;
  gwt2.setTree(&ipcaTree2);
   
  WaveletSimilarity<Precision> ws1;
  ipcaTree1.breadthFirstVisitor(&ws1);

  WaveletSimilarity<Precision> ws2;
  ipcaTree2.breadthFirstVisitor(&ws2);
  
  std::map<int, DenseVector<Precision> > &s_psi1 = ws1.getSigmaPsi();
  std::map<int, DenseVector<Precision> > &s_psi2 = ws2.getSigmaPsi();

  EuclideanMetric<Precision> l2;
  DenseMatrix<Precision> D( s_psi1.size(), s_psi2.size() );
  for(int i=0; i<s_psi1.size(); i++){
    DenseVector<Precision> v1 = s_psi1[i];
    for(int j=0; j<s_psi2.size(); j++){
      DenseVector<Precision> v2 = s_psi2[j];
      D(i, j) = l2.distance(v1, v2);
    }
  }
   
  LinalgIO<Precision>::writeMatrix(outFileName, D);

  
  return 0;

}
