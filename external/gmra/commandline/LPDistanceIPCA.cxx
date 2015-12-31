#include "Precision.h"

#include <stdio.h>


#include "LinalgIO.h"
#include "IPCAGWT.h"
#include "WaveletSimilarity.h"
#include "DenseMatrix.h"
#include <tclap/CmdLine.h>
#include "LPDistance.h"
#include "HierarchicalLPDistance.h"

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

  TCLAP::ValueArg<int> cArg("c", "constrained", "Maximum number of variables for starting to use hierarchical constraints. -1 = no constraints (default)", false, -1, "integer");
  cmd.add(cArg); 

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
   
  WassersteinNodeDistance<Precision> baseMetric;
  std::list< LPSolution > sols;

  int maxVar = cArg.getValue();
  if(maxVar == -1 ){ 
   sols = LPDistance<Precision>::distance( (GMRATree<Precision> *) &ipcaTree1, 
                                           (GMRATree<Precision> *)  &ipcaTree2, 
                                           (NodeDistance<Precision> *) &baseMetric  );
  }
  else{
   sols = HierarchicalLPDistance<Precision>::distance( (GMRATree<Precision> *) &ipcaTree1, 
                                                       (GMRATree<Precision> *)  &ipcaTree2, 
                                                       (NodeDistance<Precision> *) &baseMetric, maxVar  );
  }


  std::list< LPSolution >::iterator it = sols.begin();
  DenseVector<double> dist(sols.size());
  DenseVector<double> times(sols.size());
  DenseVector<double> errors(sols.size());
  
  for(int i=0; i < sols.size(); i++, ++it){
    LPSolution &s = *it;
    dist(i) = s.objective;
    times(i) = s.time / 1000.0;
    errors(i) = s.errorBound;

    std::stringstream ss1;
    ss1 << outFileName << "-weights-" << i << ".data";
    LinalgIO<Precision>::writeMatrix(ss1.str(), s.W);
    
    std::stringstream ss2;
    ss2 << outFileName << "-cost-" << i << ".data";
    LinalgIO<Precision>::writeMatrix(ss2.str(), s.C);
    
    std::stringstream ss3;
    ss3 << outFileName << "-pts1-" << i << ".data";
    LinalgIO<int>::writeVector(ss3.str(), s.pts1);
    
    std::stringstream ss4;
    ss4 << outFileName << "-pts2-" << i << ".data";
    LinalgIO<int>::writeVector(ss4.str(), s.pts2);


  }
  std::stringstream ss;
  ss << outFileName << ".data";
  LinalgIO<Precision>::writeVector(ss.str(), dist);
  
  std::stringstream ss2;
  ss2 << outFileName << "-times.data";
  LinalgIO<Precision>::writeVector(ss2.str(), times);

    
  std::stringstream ss3;
  ss3 << outFileName << "-errors.data";
  LinalgIO<Precision>::writeVector(ss3.str(), errors);
  return 0;

}
