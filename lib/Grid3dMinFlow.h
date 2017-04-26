#ifndef GRID3DMINFLOW_H
#define GRID3DMINFLOW_H

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include <map>
#include <vector>
#include <iostream>

#include <ilcplex/cplex.h>


template <typename TPrecision>
class Grid3dMinFlow {

  


  public:
    typedef typename Eigen::Tensor<TPrecision, 3> TensorXp;
    typedef typename Eigen::Tensor<int, 3> TensorXi;
    typedef typename Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> MatrixXi;


  private:
    TPrecision cost;
    MatrixXi id2spatial;

  public:
    
    Grid3dMinFlow(){
      cost=-1;
    };

    virtual ~Grid3dMinFlow(){
    };





    std::map< std::pair<int, int>, double > solve(const TensorXp &X1, 
        const TensorXp &X2, double lambda = 0){ 

      int status=-1;
      CPXENVptr env = CPXopenCPLEX( &status );
      CPXsetintparam( env, CPX_PARAM_SCRIND, CPX_ON );
      CPXsetintparam( env, CPX_PARAM_NETDISPLAY, 1 );
      CPXsetdblparam( env, CPX_PARAM_NETEPOPT, 1e-5 );
      CPXsetdblparam( env, CPX_PARAM_NETEPRHS, 1e-7 );
      CPXsetintparam (env, CPX_PARAM_ADVIND, 2);
      CPXsetdblparam( env, CPX_PARAM_EPPER, 1e-7 );

      CPXNETptr prob = CPXNETcreateprob (env, &status, "mmflow3d");
      
      auto dims = X1.dimensions();
      int nNodes = dims[0]*dims[1]*dims[2];
      TensorXi nId(dims[0], dims[1], dims[2]);


      std::vector<TPrecision> diff(nNodes);
      id2spatial = MatrixXi(nNodes, 3);
      
      //These statements crash
      //TensorXp X1n = X1 / X1.sum();
      //TensorXp X2n = X2 / X2.sum();
     
      Eigen::Tensor<TPrecision, 0> s1 = X1.sum();
      Eigen::Tensor<TPrecision, 0> s2 = X2.sum();
      TPrecision sum1 = s1(0);
      TPrecision sum2 = s2(0);

      int idCounter = 0;
      for(int i=0; i<dims[0]; i++){ 
        for(int j=0; j<dims[1]; j++){
          for(int k=0; k<dims[2]; k++){
            TPrecision v1 = X1(i, j, k);
            TPrecision v2 = X2(i, j, k);

            if(v1 == 0 && v2 == 0){
              nId(i, j, k) = -1;
              //continue;
            }

            nId(i, j, k) = idCounter;
            
            //diff[idCounter] = X1n(i, j, k) - X2n(i, j, k);
            diff[idCounter] = v1 / sum1 - v2 / sum2;
            
            id2spatial(idCounter, 0) = i;
            id2spatial(idCounter, 1) = j;
            id2spatial(idCounter, 2) = k;

            ++idCounter;
          }
        }
      }
     
      diff.resize(idCounter);

      std::cout << "Nodes: " << idCounter << std::endl;

      CPXNETaddnodes( env, prob, diff.size(), diff.data(), NULL);

      std::vector<int> asInd(nNodes*27);
      std::vector<int> atInd(nNodes*27);
      std::vector<double> coeff(nNodes*27);

      int index = 0;
      for(int i=0; i<dims[0]; i++){ 
        for(int j=0; j<dims[1]; j++){
          for(int k=0; k<dims[2]; k++){


            for(int ii=-1; ii < 2; ii+=1){
              int iii = i + ii;
              if(iii < 0 ){
                iii = 0; 
              }
              else if(iii == dims[0]){
                iii=dims[0]-1; 
              }

              for(int jj=-1; jj < 2; jj+=1){
                int jjj = j + jj;
                if(jjj < 0 ){
                  jjj = 0; 
                }
                else if(jjj == dims[1] ){
                  jjj = dims[1]-1;
                }

                for(int kk=-1; kk < 2; kk+=1){
                  int kkk = k + kk;
                  if(kkk < 0 ){
                    kkk = 0; 
                  }
                  else if(kkk == dims[2]){
                    kkk=dims[2]-1; 
                  }


                  if( ii == 0 && jj == 0 && kk == 0){
                    continue;
                  }
                  int id1 = nId(i, j, k);
                  int id2 = nId(iii, jjj, kkk);
                  if( id1 == -1 || id2 == -1){
                    //continue;
                  };
                  asInd[index] = id1;
                  atInd[index] = id2;

                  coeff[index] = sqrt( ii*ii + jj*jj + kk*kk);
                  if(lambda > 0){
                    coeff[index] += lambda * fabs( X1(i, j, k) - X2(iii, jjj, kkk) );
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
          plan[ std::pair<int, int>(asInd[i], atInd[i]) ] = w;
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


