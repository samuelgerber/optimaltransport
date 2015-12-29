
#include "EigenRandomSVD.h"
#include "EigenLinalgIO.h"
#include "IO.h"

#include <tclap/CmdLine.h>

int main(int argc, char **argv){

  //Command line parsing
  TCLAP::CmdLine cmd("SinkhornTransport", ' ', "1");

  TCLAP::ValueArg<std::string> xArg("x","X", "list of vector header files for streaming", true, "",
      "text file");
  cmd.add(xArg);
  
  TCLAP::ValueArg<int> pArg("p","power", "Number of power iterations", false, 3, "integer"); 
  cmd.add(pArg);
    
  
  TCLAP::ValueArg<int> dArg("d","dimension", "Number of dimensions for SVD. To get k accurate dimension set d = k+5", false, 8, "integer");
  cmd.add(dArg);
  
   TCLAP::SwitchArg mArg("m","mean", "Subtract mean", false);
  cmd.add(mArg);
  
  
  TCLAP::ValueArg<std::string> fArg("f","F", "file prefix for storing svd files and random range files", false, "row-svd", "string");
  cmd.add(fArg);

  try{ 
    cmd.parse( argc, argv );
  } 
  catch (TCLAP::ArgException &e){ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    return -1;
  }

  int p = pArg.getValue();
  int d = dArg.getValue();
  bool subtractMean = mArg.getValue();
  std::string prefix = fArg.getValue();

  std::vector<std::string> files = IO<double>::readStringList(xArg.getValue()); 

  Eigen::VectorXd mean =
      EigenLinalg::LinalgIO<double>::readVector( files[0] );
  int m = files.size();
  int n = mean.size();

  Eigen::MatrixXd X(m, n);
  for(int i=0; i< files.size(); i++){
     X.row(i) = EigenLinalg::LinalgIO<double>::readVector( files[i] );
  }


  if(subtractMean){
    Eigen::VectorXd mean = X.array().rowwise().sum() / X.cols();
    X.array().colwise() -= mean.array();
    std::stringstream ss6;
    ss6 << prefix << "-mean.data";
    EigenLinalg::LinalgIO<double>::writeVector( ss6.str(), mean );
  }
  std::cout << "Randomized svd" << std::endl;
  EigenLinalg::RandomSVD svd( X, d, p);

  std::stringstream ss1;
  ss1 << prefix << "-U.data";
  std::stringstream ss2;
  ss2 << prefix << "-S.data";
  std::stringstream ss3;
  ss3 << prefix << "-B.data";

  EigenLinalg::LinalgIO<double>::writeMatrix( ss1.str(), svd.GetU() );
  EigenLinalg::LinalgIO<double>::writeVector( ss2.str(), svd.GetS() );
  EigenLinalg::LinalgIO<double>::writeMatrix( ss3.str(), svd.GetProjected() );

  return 0;
}
