#include "Precision.h"

#include <stdio.h>

#include "CPLEXNetworkSolver.h"
#include "TransportLP.h"
#include "DenseMatrix.h"
#include "LinalgIO.h"

#include "CmdLine.h"

int main(int argc, char **argv){
  using namespace FortranLinalg;

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


  DenseMatrix<Precision> C = LinalgIO<Precision>::readMatrix(dArg.getValue());
  
  DenseVector<Precision> from(C.M());
  Linalg<Precision>::Set(from, 1.0/from.N() );

  DenseVector<Precision> to(C.N());
  Linalg<Precision>::Set(to, 1.0/to.N() );
  
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
