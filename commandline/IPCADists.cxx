#include "Precision.h"

#include <stdio.h>


#include "LinalgIO.h"
#include "IPCAGWT.h"
#include "DenseMatrix.h"
#include "CmdLine.h"

#include <set>

int main(int argc, char **argv){
  using namespace FortranLinalg;

  //Command line parsing
  TCLAP::CmdLine cmd("GMRA invariance distances", ' ', "1");

  TCLAP::ValueArg<std::string> tArg("t","tree", "Tree data file", true, "",
      "tree data file");
  cmd.add(tArg);



  TCLAP::ValueArg<std::string> outArg("o","out", "data file name", true, "", "file prefix to store results in"); 
  cmd.add(outArg);

  TCLAP::ValueArg<std::string> xArg("x","data", "matrix header file", false, "", "input data matrix"); 
  cmd.add(xArg);

  
  TCLAP::ValueArg<std::string> dArg("d","distance", "node distances", false, "", "input data matrix"); 
  cmd.add(dArg);



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
  std::string dFile = dArg.getValue();


  IPCATree<Precision> ipcaTree;
  std::ifstream treeFile(treeFileName.c_str(), std::ios_base::in | std::ios_base::binary);
  ipcaTree.unflatten(treeFile);
  treeFile.close();


  IPCAGWT<Precision> gwt;
  gwt.setTree(&ipcaTree);

  std::cout << "Constructed GWT information" << std::endl;


  DenseMatrix<Precision> X = LinalgIO<Precision>::readMatrix(xFile);
  DenseMatrix<Precision> D = LinalgIO<Precision>::readMatrix(dFile);
  DenseVector<Precision> x(X.M());


  std::vector< GWTCoefficients<Precision> > C;


  for(int i=0; i<X.N(); i++){
    Linalg<Precision>::ExtractColumn(X, i, x);
    GWTCoefficients<Precision> c = gwt.encode(x);
    C.push_back(c);
  }



  DenseMatrix<Precision> D2(X.N(), X.N() );
  for(int i=0; i<X.N(); i++){
    D2(i, i) = 0;
    GWTCoefficients<Precision> c1 = C[i];

    for(int j=i+1; j<X.N(); j++){
      GWTCoefficients<Precision> c2 = C[j];

      Precision d = 0;
      int n=0;
      for(; n< std::min(c1.ids.size(), c2.ids.size()); n++){
        Precision tmp = D(c1.ids[n], c2.ids[n]);
        d += tmp*tmp;
      }
      int nStop = n-1;
      for(; n< c1.ids.size(); n++){
        Precision tmp = D(c1.ids[n], c2.ids[nStop]);
        d += tmp*tmp;
      }
      for(; n< c2.ids.size(); n++){
        Precision tmp = D(c1.ids[nStop], c2.ids[n]);
        d += tmp*tmp;
      }

      d = sqrt(d/n);
      D2(i, j) = d;
      D2(j, i) = d;
    }

  }
  
  LinalgIO<Precision>::writeMatrix(outFileName, D2);







  return 0;

}
