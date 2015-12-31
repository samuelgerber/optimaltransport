#include "Precision.h"

#include <stdio.h>


#include "LinalgIO.h"
#include "IPCAGWT.h"
#include "DenseMatrix.h"
#include <tclap/CmdLine.h>
#include "SquaredEuclideanMetric.h"

#include <set>

int main(int argc, char **argv){
  using namespace FortranLinalg;

  //Command line parsing
  TCLAP::CmdLine cmd("GMRA Squared distances", ' ', "1");

  TCLAP::ValueArg<std::string> t1Arg("t","tree1", "Tree data file", true, "",
      "tree data file");
  cmd.add(t1Arg);

  TCLAP::ValueArg<std::string> t2Arg("u","tree2", "Tree data file", true, "",
      "tree data file");
  cmd.add(t2Arg);


  TCLAP::ValueArg<std::string> outArg("o","out", "data file name", true, "", "file prefix to store results in"); 
  cmd.add(outArg);

  TCLAP::ValueArg<std::string> xArg("x","xdata", "matrix header file", false, "", "input data matrix"); 
  cmd.add(xArg);

  
  TCLAP::ValueArg<std::string> yArg("y","ydata", "matrix header file", false, "", "input data matrix"); 
  cmd.add(yArg);


  try{ 
    cmd.parse( argc, argv );
  } 
  catch (TCLAP::ArgException &e){ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    return -1;
  }


  std::string outFileName = outArg.getValue();
  std::string treeFileName1 = t1Arg.getValue();
  std::string treeFileName2 = t2Arg.getValue();
  std::string xFile = xArg.getValue();
  std::string yFile = yArg.getValue();


  IPCATree<Precision> ipcaTree1;
  std::ifstream treeFile1(treeFileName1.c_str(), std::ios_base::in | std::ios_base::binary);
  ipcaTree1.unflatten(treeFile1);
  treeFile1.close();


  IPCATree<Precision> ipcaTree2;
  std::ifstream treeFile2(treeFileName2.c_str(), std::ios_base::in | std::ios_base::binary);
  ipcaTree2.unflatten(treeFile2);
  treeFile2.close();


  IPCAGWT<Precision> gwt1;
  gwt1.setTree(&ipcaTree1);

  IPCAGWT<Precision> gwt2;
  gwt2.setTree(&ipcaTree2);

  std::cout << "Constructed GWT information" << std::endl;


  DenseMatrix<Precision> X1 = LinalgIO<Precision>::readMatrix(xFile);
  DenseMatrix<Precision> X2 = LinalgIO<Precision>::readMatrix(yFile);
  DenseVector<Precision> x1(X1.M());

  DenseVector<Precision> x2(X2.M());

  std::vector< GWTCoefficients<Precision> > C1;
  std::vector< GWTCoefficients<Precision> > C2;


  for(int i=0; i<X1.N(); i++){
    Linalg<Precision>::ExtractColumn(X1, i, x1);
    GWTCoefficients<Precision> c = gwt1.encode(x1);
    C1.push_back(c);
  }

  for(int i=0; i<X2.N(); i++){
    Linalg<Precision>::ExtractColumn(X2, i, x2);
    GWTCoefficients<Precision> c = gwt2.encode(x2);
    C2.push_back(c);
  }


  DenseMatrix<Precision> D(X1.N(), X1.N() );
  for(int i=0; i<X1.N(); i++){
    D(i, i) = 0;
    GWTCoefficients<Precision> &c1 = C1[i];
    std::map<int, DenseVector<Precision> > z1;
    for(int k =0; k<c1.ids.size(); k++){
      int id = c1.ids[k];

      GWTCoefficients<Precision> d1 = C2[id];
      for(int l=0; l<d1.ids.size(); l++){
        std::map<int, DenseVector<Precision> >::iterator it =
          z1.find(d1.ids[l]);
        DenseVector<Precision> v;
        if(it == z1.end() ){
          v  = DenseVector<Precision>(d1.coeff[l].N() );
          Linalg<Precision>::Zero(v);
          z1[d1.ids[l]] = v;
        } 
        else{
          v = it->second;
        }
        Linalg<Precision>::AddScale(v,
            Linalg<Precision>::Sum(c1.coeff[k]), d1.coeff[l], v);
      }
    }

    for(int j=i+1; j<X1.N(); j++){

      GWTCoefficients<Precision> c2 = C1[j];
      std::map<int, DenseVector<Precision> > z2;
      for(int k =0; k<c2.ids.size(); k++){
        int id = c2.ids[k];

        GWTCoefficients<Precision> d2 = C2[id];
        for(int l=0; l<d2.ids.size(); l++){
          std::map<int, DenseVector<Precision> >::iterator it = z2.find(d2.ids[l]);
          DenseVector<Precision> v;
          if(it == z2.end() ){
            v  = DenseVector<Precision>(d2.coeff[l].N() );
            Linalg<Precision>::Zero(v);
            z2[d2.ids[l]] = v;
          } 
          else{
            v = it->second;
          }
          Linalg<Precision>::AddScale(v, Linalg<Precision>::Sum(c2.coeff[k]),
              d2.coeff[l], v);
        }
      }

      Precision tmp = 0;
      std::map<int, DenseVector<Precision> >::iterator it1 = z1.begin();
      std::map<int, DenseVector<Precision> >::iterator it2 = z2.begin();
      while(it1!=z1.end() && it2 != z2.end()){
        if(it1->first < it2->first){
          tmp += Linalg<Precision>::SquaredLength(it1->second);
          ++it1;
        }
        else if(it1->first > it2->first){
          tmp += Linalg<Precision>::SquaredLength(it2->second);
          ++it2;
        }
        else{
          static SquaredEuclideanMetric<Precision> sl2;
          tmp += sl2.distance(it1->second, it2->second);
          ++it1;
          ++it2;
        }
      }
      while(it1!=z1.end()){
        tmp += Linalg<Precision>::SquaredLength(it1->second);
        ++it1;
      }
      while(it2!=z2.end()){
        tmp += Linalg<Precision>::SquaredLength(it2->second);
        ++it2;
      }

      tmp = sqrt(tmp);
      D(i, j) = tmp;
      D(j, i) = tmp;
    }

  }
  
  
  LinalgIO<Precision>::writeMatrix(outFileName, D);




  return 0;

}
