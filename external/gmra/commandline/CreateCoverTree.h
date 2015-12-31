#define VERBOSE


#include "Precision.h"

#include <stdio.h>

#include "LinalgIO.h"
#include "IPCATree.h"
#include "DenseMatrix.h"
#include "CmdLine.h"


int main(int argc, char **argv){
  
  //Command line parsing
  TCLAP::CmdLine cmd("GWT", ' ', "1");

  TCLAP::ValueArg<int> dArg("d","dim","Dimension of GWT", true, 3, "integer");
  cmd.add(dArg);

  TCLAP::ValueArg<Precision> tArg("t", "threshold", "Dimensionality thresholds (MSE), negative uses fixed dimensionality", true, -1, "double");
  cmd.add(tArg);

  
  TCLAP::ValueArg<Precision> eArg("e","epsilon","Precision of GWT", true, 0.1, "double");
  cmd.add(eArg);
  
  TCLAP::ValueArg<int> sArg("s","stop","1: Local R^2, 2: Total R^2, 3: Node MSE, 4: Node Radius, 5: Relative node radius", true, 2, "integer");
  cmd.add(sArg);

  TCLAP::ValueArg<std::string> dataArg("x","data", "Data file",  true, "", "matrix header file");
  cmd.add(dataArg);

  TCLAP::ValueArg<std::string> outArg("o","out", "Data file",  true, "", "file to store tree in");
  cmd.add(outArg);

  try{
	  cmd.parse( argc, argv );
	} 
  catch (TCLAP::ArgException &e){ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    return -1;
  }


  FortranLinalg::DenseMatrix<Precision> X = FortranLinalg::LinalgIO<Precision>::readMatrix(dataArg.getValue());
  int d = dArg.getValue();
  Precision eps = eArg.getValue();
  std::string fileName = outArg.getValue();

  NodeFactory<Precision> *nf;
  Precision t = tArg.getValue();
  if( t > 0 ){
      nf = new RelativePrecisionNodeFactory<Precision>(d,t);
  }
  else{
    nf = new FixedNodeFactory<Precision>(d); 
  }

  IPCATree<Precision>::StoppingCriterium stop = IPCATree<Precision>::TOTAL_R2;
  if(sArg.getValue() == 1){
    stop = IPCATree<Precision>::NODE_R2;
  }
  else if(sArg.getValue() == 2){
    stop = IPCATree<Precision>::TOTAL_R2;
  }
  else if(sArg.getValue() == 3){
    stop = IPCATree<Precision>::NODE_MSE;
  }
  else if(sArg.getValue() == 4){
    stop = IPCATree<Precision>::NODE_RADIUS;
  }
  else if(sArg.getValue() == 5){
    stop = IPCATree<Precision>::RELATIVE_NODE_RADIUS;
  }

  IPCATree<Precision> ipcaTree(d, eps, stop, nf);
  ipcaTree.construct(X.data(), X.M(), X.N() );

  std::ofstream file(fileName.c_str(), std::ios_base::out | std::ios_base::binary);
  ipcaTree.flatten(file);

  file.close();

  X.deallocate();

  return 0;

}
