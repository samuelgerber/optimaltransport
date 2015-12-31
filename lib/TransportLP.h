#ifndef TRANSPORTLP_H
#define TRANSPORTLP_H

#include "LPSolver.h"
#include <Eigen/Dense>
#include <map>



template <typename TPrecision>
class TransportLP {

  




  private:
    
    LPSolver *solver;
    typedef LPSolver::Status Status;


  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
    
    
    TransportLP(LPSolver *lps) : solver(lps) {
    };

    virtual ~TransportLP(){
    };





    std::map< std::pair<int, int>, double > solve(const MatrixXp &C, const VectorXp &from, const VectorXp &to){

        solver->createLP(C.rows(), C.cols());

        solver->addRows( C.rows() + C.cols() );
        solver->addCols( C.rows()*C.cols() );

        //set constraints bounds
        for(int i=0; i<C.rows(); i++){ 
          solver->setRowBounds(i, from(i) );
        }
        for(int i=0; i<C.cols(); i++){ 
          solver->setRowBounds(C.rows()+i, -to(i) );
        }


        int index = 0;
        for(int i = 0; i < C.rows(); i++){
          for(int j = 0; j < C.cols(); j++){
            solver->setColBounds(index);
            solver->setCoefficent(index, C(i, j) );
            solver->setColConstraints(index, i, C.rows() + j ); 
            ++index;
          }
        }


        //solver->setupStdBasis();
        solver->solveLP();        
        
        std::map<std::pair<int, int>, double> plan;
        index = 0;
        for(int i = 0; i < C.rows(); i++){
          for(int j = 0; j < C.cols(); j++){
            double w = solver->getColPrimal(index);
            if(w > 0){
              plan[ std::pair<int, int>(i, j) ] = w;
            }
            ++index;
          }
        }

        return plan;

      };




};


#endif


