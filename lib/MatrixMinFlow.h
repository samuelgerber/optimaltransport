#ifndef MATRIXMINFLOW_H
#define MATRIXMINFLOW_H

#include <Eigen/Dense>
#include <map>
#include <vector>
#include <iostream>

#include <ilcplex/cplex.h>


template <typename TPrecision>
class MatrixMinFlow {

  


  public:

    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> MatrixXi;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;




  private:

    TPrecision cost;
    MatrixXi id2spatial;




  public:
    
    MatrixMinFlow(){
      cost=-1;
    };


    virtual ~MatrixMinFlow(){
    };





    std::map< std::pair<int, int>, double > solve(const MatrixXp &X1, 
        const MatrixXp &X2, const MatrixXp &weights, double lambda = 0){ 

      int status=-1;
      CPXENVptr env = CPXopenCPLEX( &status );
      CPXsetintparam( env, CPX_PARAM_SCRIND, CPX_ON );
      CPXsetintparam( env, CPX_PARAM_NETDISPLAY, 1 );
      CPXsetdblparam( env, CPX_PARAM_NETEPOPT, 1e-5 );
      CPXsetdblparam( env, CPX_PARAM_NETEPRHS, 1e-7 );
      CPXsetintparam( env, CPX_PARAM_ADVIND, 2 );
      CPXsetdblparam( env, CPX_PARAM_EPPER, 1e-7 );

      CPXNETptr prob = CPXNETcreateprob( env, &status, "mmflow" );
      
      int nNodes = X1.rows() * X1.cols();

      MatrixXi nId( X1.rows(), X1.cols() );
      std::vector<TPrecision> diff(nNodes);
      id2spatial = MatrixXi(nNodes, 2);
      
      MatrixXp X1n = X1.array() * weights.array(); 
      X1n /= X1n.sum();
      MatrixXp X2n = X2.array() * weights.array();
      X2n /= X2n.sum();
      
      int idCounter = 0;
      for(int i=0; i<X1.rows(); i++){ 
        for(int j=0; j<X1.cols(); j++){
          if(weights(i, j) != 0 ){
            nId(i, j) = idCounter;
            diff[idCounter] = X1n(i, j) - X2n(i, j);
            id2spatial(idCounter, 0) = i;
            id2spatial(idCounter, 1) = j;
            ++idCounter;
          }
          else{
            nId(i, j) = -1;
          }
        }
      } 
      
      nNodes = idCounter;
      diff.resize(nNodes);
      id2spatial = id2spatial.topRows(nNodes);

      CPXNETaddnodes( env, prob, diff.size(), diff.data(), NULL);

      std::vector<int> asInd(nNodes*9);
      std::vector<int> atInd(nNodes*9);
      std::vector<double> coeff(nNodes*9);

      int index = 0;
      for(int i=0; i<X1.rows(); i++){
        for(int j=0; j<X1.cols(); j++){
          

          if( nId(i, j) != -1 ){


            for(int ii=-1; ii< 2; ii+=1){
              int iii = i + ii;
              if(iii < 0 ){
                iii = 0; //X1.rows()-1;
              }
              else if(iii == X1.rows()){
                iii=X1.rows()-1; //0;
              }
              for(int jj=-1; jj< 2; jj+=1){
                int jjj = j + jj;
                if(jjj < 0 ){
                  jjj = 0; //X1.cols()-1;
                }
                else if(jjj == X1.cols()){
                  jjj=X1.cols()-1; //0;
                }

                if( nId(iii, jjj) != -1 ){
                  asInd[index] = nId(i, j);
                  atInd[index] = nId(iii, jjj);

                  coeff[index] = sqrt( ii*ii + jj*jj );
                  if(lambda > 0){
                    //TPrecision idelta = std::min( fabs( X1(i, j) - X1(iii, jjj) ), fabs( X2(i, j) - X2(iii, jjj) ) );
                    //idelta = std::min( fabs( X1(i, j) - X2(iii, jjj) ), idelta );
                    TPrecision idelta = fabs( X1(i, j) - X2(iii, jjj) );
                    coeff[index] += lambda * idelta;
                  }

                  ++index;
                }
              }
            }


          }


        }
      }
     
      asInd.resize(index);
      atInd.resize(index);
      coeff.resize(index);

      CPXNETaddarcs( env, prob, asInd.size(), asInd.data(), atInd.data(), 
                     NULL, NULL, coeff.data(), NULL );

 
      CPXNETchgobjsen( env, prob, CPX_MIN );
       
      status = CPXNETprimopt( env, prob );
      std::cout << "CPXNET optimization status: " << status << std::endl;

      std::vector<double> primal( coeff.size(), 0 );
      std::vector<double> dual( diff.size(), 0 );
      int solstat = 0;
      cost = -1;
      status = CPXNETsolution( env, prob, &solstat, &cost, primal.data(),
         dual.data(), NULL, NULL); 
      std::cout << "CPXNET solution status: " << status << std::endl;

      std::map<std::pair<int, int>, double> plan;
      index = 0;
      for(int i = 0; i < asInd.size(); i++){
        double w = primal[ i ];
        if(w > 0){
          plan[ std::pair<int, int>( asInd[i], atInd[i] ) ] = w;
        }
      }

      CPXNETfreeprob(env, &prob);
      CPXcloseCPLEX (&env);
      
  
      return plan;
    };




    TPrecision getCost(){
      return cost;
    };




    MatrixXi &getId2Spatial(){
      return id2spatial;
    };



};


#endif


