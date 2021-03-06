#include "SparseSinkhornTransport.h"

#include "EigenEuclideanMetric.h"
#include "EigenLinalgIO.h"

#include <tclap/CmdLine.h>

int main(int argc, char **argv){

  //Command line parsing
  TCLAP::CmdLine cmd("SinkhornTransport", ' ', "1");

  TCLAP::ValueArg<std::string> muArg("","mu", "Source point set", true, "",
      "matrix header file");
  cmd.add(muArg);

  TCLAP::ValueArg<std::string> nuArg("","nu", "Target point set", true, "",
      "matrix header file");
  cmd.add(nuArg);


  TCLAP::ValueArg<std::string> outArg("o","out", "data file name", true, "", "file to store results in"); 
  cmd.add(outArg);
  
  TCLAP::ValueArg<int> iArg("i","iter", "Maximum number of iterations", false, 1000, "integer"); 
  cmd.add(iArg);

  TCLAP::ValueArg<double> lArg("l","lambda", "Entropy penalty actor", true, 1, "real"); 
  cmd.add(lArg);

  TCLAP::ValueArg<double> tArg("t","tolerance", "Tolerance for stopping", false,
      0.0001, "real"); 
  cmd.add(tArg);
  
  TCLAP::ValueArg<double> zArg("z","threhsold", "Threshold for zero entreis", false,
      0.00001, "real"); 
  cmd.add(zArg);

  try{ 
    cmd.parse( argc, argv );
  } 
  catch (TCLAP::ArgException &e){ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    return -1;
  }

  double lambda = lArg.getValue();
  double tol = tArg.getValue();
  double thres = zArg.getValue();
  int iter = iArg.getValue();

  std::string muFilename = muArg.getValue();
  std::string nuFilename = nuArg.getValue();

  Eigen::MatrixXd Xmu =
    EigenLinalg::LinalgIO<double>::readMatrix( muFilename );
  Eigen::MatrixXd Xnu =
    EigenLinalg::LinalgIO<double>::readMatrix( nuFilename );



  EuclideanMetric<double> l2metric;

  Eigen::SparseMatrix<double> K( Xmu.cols(), Xnu.cols() );
  Eigen::SparseMatrix<double> U( Xmu.cols(), Xnu.cols() );
  K.reserve( ( Xmu.cols()*Xnu.cols()) );
  U.reserve( ( Xmu.cols()*Xnu.cols()) );
  for(int i=0; i<Xmu.cols(); i++){
    for(int j=0; j<Xnu.cols(); j++){
      double d =l2metric.distance(Xmu.col(i), Xnu.col(j)); 
      double e = exp(-lambda*d);
      if( e > thres){
        K.insert(i,j) = e;
        U.insert(i,j) = e * d;
      }
    }
  }

  std::cout << "size: " << K.cols()*K.rows() << std::endl;
  std::cout << "nnz: " << K.nonZeros() << std::endl;

  Eigen::VectorXd mu = Eigen::VectorXd::Constant( Xmu.cols(),    1.0 / Xmu.cols() );
  Eigen::MatrixXd nu = Eigen::MatrixXd::Constant( Xnu.cols(), 2, 1.0 / Xnu.cols() );
  SparseSinkhornTransport sinkhorn;
  sinkhorn.transport(mu, nu, K, U, tol, iter); 

  std::string outFileName = outArg.getValue();

  std::cout << sinkhorn.getDistances() << std::endl;
  return 0;
}
