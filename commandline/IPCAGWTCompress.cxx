#include "Precision.h"

#include <stdio.h>


#include "LinalgIO.h"
#include "IPCAGWT.h"
#include "CompressGWT.h"
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

  TCLAP::ValueArg<std::string> xArg("x","data", "matrix header file", true, "", "input data matrix"); 
  cmd.add(xArg);

  TCLAP::ValueArg<double> sArg("","thres", "threshold for combinig planes", 0, true, "double"); 
  cmd.add(sArg);

  try{ 
    cmd.parse( argc, argv );
  } 
  catch (TCLAP::ArgException &e){ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    return -1;
  }


  std::string outFileName = outArg.getValue();
  std::string treeFileName = tArg.getValue();
  std::string xFile = xArg.getValue();
  double thres = sArg.getValue();

  DenseMatrix<Precision> X = LinalgIO<Precision>::readMatrix(xFile);

  IPCATree<Precision> ipcaTree;

  std::ifstream treeFile(treeFileName.c_str(), std::ios_base::in | std::ios_base::binary);
  ipcaTree.unflatten(treeFile);

  treeFile.close();

  std::cout << "Read tree " << std::endl;


  IPCAGWT<Precision> gwt;
  gwt.setTree(&ipcaTree);

  CompressGWT<Precision> *compress = new CompressGWT<Precision>(X, thres);
  ipcaTree.breadthFirstVisitor(compress);

  std::cout << "Compress size: " << compress->dict.size() << std::endl;
  std::cout << "Original size: " << gwt.getNumberOfNodes() << std::endl;
   
  return 0;

}
