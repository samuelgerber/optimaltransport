
#include "EigenStreamingRandomSVD.h"
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
  
  
  TCLAP::ValueArg<std::string> fArg("f","F", "file prefix for storing svd files and random range files", false, "row-streaming-svd", "string");
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

  std::vector<std::string> files = IO<double>::readStringList(xArg.getValue()); 

  Eigen::VectorXd tmp =
      EigenLinalg::LinalgIO<double>::readVector( files[0] );
  int m = files.size();
  int n = tmp.size();

  bool subtractMean =mArg.getValue();
  Eigen::VectorXd mean(m);


  std::cout << "Computing randomized range"<< std::endl;
  EigenLinalg::RowStreamingRandomRange<double> range(m, n, d);
  for(int i=0; i< files.size(); i++){
    Eigen::VectorXd x =
      EigenLinalg::LinalgIO<double>::readVector( files[i] );
    mean(i) = x.sum() / x.size();
    if(subtractMean){
      x.array() -= mean(i);
    }
    range.Add(x);
    std::cout << i << "," << std::flush;
  }
  std::cout << std::endl;
  std::cout << "QR decomposition of range" << std::endl;
  range.Compute();

  std::cout << "Randomized svd" << std::endl;
  EigenLinalg::RowStreamingRandomSVD<double> svd( range, n );
  for(int i=0; i< files.size(); i++){
    Eigen::VectorXd x =
      EigenLinalg::LinalgIO<double>::readVector( files[i] );
    if(subtractMean){
      x.array() -= mean(i);
    }
    svd.Add(x);
    std::cout << i << "," << std::flush;
  }
  std::cout << std::endl;
  std::cout << "Doing svd" << std::endl;
  svd.Compute();

  std::string prefix = fArg.getValue();
  std::stringstream ss1;
  ss1 << prefix << "-U.data";
  std::stringstream ss2;
  ss2 << prefix << "-S.data";
  std::stringstream ss3;
  ss3 << prefix << "-B.data";
  std::stringstream ss4;
  ss4 << prefix << "-N.data";
  std::stringstream ss5;
  ss5 << prefix << "-Y.data";

  EigenLinalg::LinalgIO<double>::writeMatrix( ss1.str(), svd.GetU() );
  EigenLinalg::LinalgIO<double>::writeVector( ss2.str(), svd.GetS() );
  EigenLinalg::LinalgIO<double>::writeMatrix( ss3.str(), svd.GetProjected() );
  EigenLinalg::LinalgIO<double>::writeMatrix( ss4.str(), range.GetRandomProjector() );
  EigenLinalg::LinalgIO<double>::writeMatrix( ss5.str(), range.GetProjected() );

  if(subtractMean){
    std::stringstream ss6;
    ss6 << prefix << "-mean.data";
    EigenLinalg::LinalgIO<double>::writeVector( ss6.str(), mean );
  }

  return 0;
}
