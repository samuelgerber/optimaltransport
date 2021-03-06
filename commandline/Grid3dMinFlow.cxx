#include "Precision.h"

#include <stdio.h>

#include "Grid3dMinFlow.h"
#include <Eigen/Dense>
#include "EigenLinalgIO.h"

#include <tclap/CmdLine.h>

int main(int argc, char **argv){

  //Command line parsing
  TCLAP::CmdLine cmd("Matrix Minimum FLow", ' ', "1");

  TCLAP::ValueArg<std::string> sArg("s","source", "2D source distribution" ,  true, "", "matrix header file");
  cmd.add(sArg);

  TCLAP::ValueArg<std::string> tArg("t","target", "2D target distribution" ,  true, "", "matrix header file");
  cmd.add(tArg);

  try{
	  cmd.parse( argc, argv );
	} 
  catch (TCLAP::ArgException &e){ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    return -1;
  }


  //Eigen::MatrixXd S = EigenLinalg::LinalgIO<Precision>::readMatrix(sArg.getValue());
  //Eigen::MatrixXd T = EigenLinalg::LinalgIO<Precision>::readMatrix(tArg.getValue());
  Grid3dMinFlow<Precision>::TensorXp S;
  Grid3dMinFlow<Precision>::TensorXp T;
   
  Grid3dMinFlow<Precision> mmf;

  std::map< std::pair<int, int>, Precision> plan = mmf.solve(S, T);
 
 
  return 0;

}
