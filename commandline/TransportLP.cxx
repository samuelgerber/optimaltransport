#include "Precision.h"

#include <stdio.h>

#include "CPLEXNetworkSolver.h"
#include "TransportLP.h"
#include <Eigen/Dense>
#include "EigenLinalgIO.h"

#include <tclap/CmdLine.h>

int main(int argc, char **argv){

  //Command line parsing
  TCLAP::CmdLine cmd("LPTransport", ' ', "1");

  TCLAP::ValueArg<std::string> dArg("d","distances", "cross distance from source to target",  true, "", "matrix header file");
  cmd.add(dArg);

  try{
	  cmd.parse( argc, argv );
	} 
  catch (TCLAP::ArgException &e){ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    return -1;
  }


  Eigen::MatrixXd C = EigenLinalg::LinalgIO<Precision>::readMatrix(dArg.getValue());
  
  Eigen::VectorXd from= Eigen::VectorXd::Constant( C.rows(), 1.0/C.rows() );

  Eigen::VectorXd to= Eigen::VectorXd::Constant( C.cols(), 1.0/C.cols() );
  
  LPSolver *solver= new CPLEXNetworkSolver();
  TransportLP<Precision> trp(solver);

  std::map< std::pair<int, int>, Precision> plan = trp.solve(C, from, to);
 
 
  Precision cost = 0;
  for(std::map< std::pair<int, int>, Precision>::iterator it = plan.begin(); it!=plan.end(); ++it){
    const std::pair<int, int> &p = it->first;
    cost += C(p.first, p.second) * it->second;
  }

  std::cout << cost << std::endl;

  
    
  return 0;

}
