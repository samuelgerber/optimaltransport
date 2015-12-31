#include "Precision.h"

#include <stdio.h>


#include "LinalgIO.h"
#include "IPCAGWT.h"
#include "WaveletSimilarity2.h"
#include "DenseMatrix.h"
#include <tclap/CmdLine.h>

int main(int argc, char **argv){

  //Command line parsing
  TCLAP::CmdLine cmd("IPCA GWT", ' ', "1");

  TCLAP::ValueArg<std::string> tArg("t","tree", "Tree data file", true, "",
      "tree data file");
  cmd.add(tArg);

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

  IPCATree<Precision> ipcaTree;

  std::ifstream treeFile(treeFileName.c_str(), std::ios_base::in | std::ios_base::binary);
  ipcaTree.unflatten(treeFile);

  treeFile.close();

  std::cout << "read tree " << std::endl;


  IPCAGWT<Precision> gwt;
  gwt.setTree(&ipcaTree);

   
  WaveletSimilarity2<Precision> ws;

  ipcaTree.breadthFirstVisitor(&ws);

  DenseMatrix<Precision> D = ws.distances();
   
  LinalgIO<Precision>::writeMatrix(outFileName, D);

  
  return 0;

}
