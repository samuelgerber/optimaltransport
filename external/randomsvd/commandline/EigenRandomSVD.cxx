
#include "EigenRandomSVD.h"
#include "EigenLinalgIO.h"

#include <tclap/CmdLine.h>

int main(int argc, char **argv){

  //Command line parsing
  TCLAP::CmdLine cmd("Eigen SVD", ' ', "1");

  TCLAP::ValueArg<std::string> xArg("x","X", "Matrix to decompose", true, "",
      "matrix header file");
  cmd.add(xArg);
  
  TCLAP::ValueArg<int> pArg("p","power", "Number of power iterations", false, 3, "integer"); 
  cmd.add(pArg);
  
  TCLAP::ValueArg<int> dArg("d","dimension", "Number of dimensions for SVD. To get k accurate dimension set d = k+5", false, 3, "integer");
  cmd.add(dArg);

  try{ 
    cmd.parse( argc, argv );
  } 
  catch (TCLAP::ArgException &e){ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    return -1;
  }

  int p = pArg.getValue();
  int d = dArg.getValue();

  Eigen::MatrixXd X = EigenLinalg::LinalgIO<double>::readMatrix( xArg.getValue() );
  
  
  clock_t t1 = clock();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
  clock_t t2 = clock();
  std::cout << "SVD" << std::endl;
  std::cout << (t2-t1)/(double)CLOCKS_PER_SEC << std::endl;
  Eigen::MatrixXd U = svd.matrixU();
  Eigen::VectorXd S = svd.singularValues();
  EigenLinalg::LinalgIO<double>::writeVector("ES.data", S );
  EigenLinalg::LinalgIO<double>::writeMatrix("EU.data", U );
  
  t1 = clock();
  EigenLinalg::RandomSVD<double> rsvd(X, d, p);
  t2 = clock();
  std::cout << "Random SVD" << std::endl;
  std::cout << (t2-t1)/(double)CLOCKS_PER_SEC << std::endl;
  EigenLinalg::LinalgIO<double>::writeVector("ErS.data", rsvd.GetS());
  EigenLinalg::LinalgIO<double>::writeMatrix("ErU.data", rsvd.GetU());
  return 0;
}
