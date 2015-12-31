

#ifndef CPLEXNETWORKSOLVER_H
#define CPLEXNETWORKSOLVER_H

#include "LPSolver.h"

#include <ilcplex/cplex.h>
#include <stdlib.h>
#include <iostream>

#include <vector>


class CPLEXNetworkSolver : public LPSolver{

  private:
    typedef  LPSolver::Status Status;
 

    CPXENVptr env;
    
    std::vector<int> sInd;
    std::vector<int> tInd;

    std::vector<double> coeff;
    std::vector<double> mass;
    std::vector<double> primal;
    std::vector<double> dual;
    std::vector<double> colLB;
    std::vector<double> colUB;

    std::vector<int> colStatus;
    std::vector<int> rowStatus;

    int solstat;
    int iCount;
    double objValue;
    double lambda;
    int ns;
    int nt;


  public:

    CPLEXNetworkSolver(double l=-1) {

        int status=-1;
        env = CPXopenCPLEX( &status );
        CPXsetintparam( env, CPX_PARAM_SCRIND, CPX_ON );
        CPXsetintparam( env, CPX_PARAM_NETDISPLAY, 1 );
        CPXsetdblparam( env, CPX_PARAM_NETEPOPT, 1e-8 );
        CPXsetdblparam( env, CPX_PARAM_NETEPRHS, 1e-8 );
        CPXsetintparam (env, CPX_PARAM_ADVIND, 2);
        lambda = l;
        ns = 0;
        solstat = CPX_STAT_ABORT_USER;
    };


    ~CPLEXNetworkSolver(){
      deleteLP();
      CPXcloseCPLEX (&env);
    };





   virtual void createLP(int nSource, int nTarget){
     deleteLP();
     ns=nSource;
     nt=nTarget;
     
     if(lambda > 0){
       rowStatus.resize( 1, CPX_BASIC );
       dual.resize( 1, 0);
       
       primal.resize( nt, 0);
       colStatus.resize( nt, CPX_AT_LOWER );
     }

   };



   virtual void solveLP(){

     int status;
     CPXNETptr prob = CPXNETcreateprob (env, &status, "netex1");

     if(lambda > 0 ){
       std::vector<double> mtmp = mass;
       //Set target nodes to 0, i.e. flow in =  flow out
       for(unsigned int i=ns; i<mtmp.size(); i++){
         mtmp[i] = 0;
       }
       //Add a terminal node with all the mass
       mtmp.push_back(-1);

       //Add nodes
       CPXNETaddnodes( env, prob, mtmp.size(), mtmp.data(), NULL);
       
       //add arcs from source to target nodes with bounded capacity
       CPXNETaddarcs( env, prob, sInd.size(), sInd.data(), tInd.data(), colLB.data(), colUB.data(),
          coeff.data(), NULL );

       //add arc from target nodes to terminal node with upper and lower bounds accoridng to lambda *
       //target mass.
       //Thus each target node is within delta = m*lambda of its mass
       std::vector<int> ttInd(nt, mass.size());
       std::vector<int> tsInd( nt );
       std::vector<double> lb(nt);
       std::vector<double> ub(nt);

       for(int i=0; i<nt; i++){
         double m = fabs(mass[ns+i]);
         double delta = m*lambda; 
         lb[i] =  m - delta;
         ub[i] =  m + delta;
         tsInd[i] = ns + i;       
       }


       CPXNETaddarcs( env, prob, tsInd.size(), tsInd.data(), ttInd.data(), lb.data(), ub.data(),
          NULL, NULL );
     } 
     else{
       CPXNETaddnodes( env, prob, mass.size(), mass.data(), NULL);
       CPXNETaddarcs( env, prob, sInd.size(), sInd.data(), tInd.data(), colLB.data(), colUB.data(), 
          coeff.data(), NULL );
     }

     CPXNETchgobjsen( env, prob, CPX_MIN );
     CPXNETcopybase( env, prob, colStatus.data(), rowStatus.data() );
    

     status = CPXNETprimopt (env, prob);
     //std::cout << "CPXNET optimization status: " << status << std::endl;

     status = CPXNETsolution( env, prob, &solstat, &objValue, primal.data(),
         dual.data(), NULL, NULL);

     //std::cout << "CPXNET get solution status: " << status << std::endl;
     //std::cout << "CPXNET solution status: " << solstat << " / " <<
     //CPX_STAT_OPTIMAL << ", " <<  CPX_STAT_INFEASIBLE << std::endl;

     if(solstat != CPX_STAT_OPTIMAL ){
       //std::cout << "Arg" << std::endl;
     }

     iCount = CPXNETgetitcnt( env, prob);
     
     CPXNETgetbase( env, prob, colStatus.data(), rowStatus.data());
    
    /* 
     for(int i=0; i<colStatus.size(); i++){
       std::cout << colStatus[i] << ", ";
     }
     std::cout << std::endl;
     std::cout << std::endl;
     */

     CPXNETfreeprob(env, &prob);

     //std::cout << "Obj: " << objValue << std::endl;
     //std::cout << "#iter: " << iCount << std::endl << std::endl;
   };



   virtual void deleteLP(){
     sInd.clear();
     tInd.clear();
     coeff.clear();
     mass.clear();
     dual.clear();
     primal.clear();
     rowStatus.clear();
     colStatus.clear();
     colLB.clear();
     colUB.clear();
     ns = 0;
   };


   virtual bool isOptimal(){
     return solstat == CPX_STAT_OPTIMAL;
   };

   virtual void addCols(int n){
     sInd.resize( sInd.size() + n, -1  );
     tInd.resize( tInd.size() + n, -1 );
     coeff.resize( coeff.size() + n, 0 );
     primal.resize( primal.size() + n, 0 );
     colStatus.resize( colStatus.size() + n, CPX_AT_LOWER );
     colLB.resize( colLB.size() + n, 0 );
     colUB.resize( colUB.size() + n, CPX_INFBOUND );
   };



   virtual void addRows(int n){
     mass.resize( mass.size() + n , 0);
     rowStatus.resize( rowStatus.size() + n, CPX_BASIC );
     dual.resize( dual.size() + n, 0);
   };
 


   virtual double getRowDual(int row){
     return dual[row];
   };



   virtual double getColPrimal(int col){
     return primal[col];
   };
   

   virtual void setRowBounds(int i, double m){
     mass[i] = m;
   };


   virtual double getRowBounds(int i){
     return mass[i]; 
   };

   virtual void setColBounds(int i){
     setColBounds(i, 0, CPX_INFBOUND);
   };

   virtual void setColBounds(int col, double lb, double ub){
     colLB[col] = lb;
     colUB[col] = ub;
   };

   virtual Status getColStatus(int col){
     //std::cout << "get: " << colStatus[col] << std::endl;
     return convertFromCPLEX( colStatus[col] );

   };
   
   virtual Status getRowStatus(int row){
     return convertFromCPLEX( rowStatus[row] );
   };
   
   virtual void setColStatus(int col, Status s){
     //std::cout << s << " -> ";
     //std::cout << convertToCPLEX(s) << std::endl;
     colStatus[col] = convertToCPLEX(s);
   };

   
   virtual void setRowStatus(int row, Status s){
     if(s == LPSolver::BASIC){
       rowStatus[row] = CPX_BASIC;
     }
     else{
       rowStatus[row] = CPX_AT_LOWER;
     }
   };

   
   virtual void setCoefficent(int i, double cost){
     coeff[i] = cost;
   };

   virtual void setColConstraints(int col, int s, int t){
     ns = std::max(ns, s+1);
     tInd[col] = t;
     sInd[col] = s;
   };
   
   virtual int getColConstraints(int col, int *ind, double *val){
     ind[0] = sInd[col];
     ind[1] = tInd[col];
     val[0] = 1;
     val[1] = -1;
     return 2;    
   };
   

   virtual double getObjectiveValue(){
     return objValue;
   };


   virtual int getIterationCount(){
     return iCount;
   };


   virtual int getNumRows(){     
     return dual.size();
   };


   virtual int getNumCols(){
     return primal.size();
   };
   
   virtual void setupStdBasis(){
     for(int i= 0; i< getNumCols(); i++){
       colStatus[i] = CPX_AT_LOWER;
     }
     for(int i= 0; i< getNumRows(); i++){
       rowStatus[i] =  CPX_BASIC;
     }
   };



  private:




   Status convertFromCPLEX(int s){

     switch(s){
       case CPX_BASIC:
         return LPSolver::BASIC;
       case CPX_AT_UPPER:
         return LPSolver::UPPER;
       case CPX_AT_LOWER:
         return LPSolver::LOWER;
       case CPX_FREE_SUPER:
         return LPSolver::FREE;
     }
     return LPSolver::END; 
   };




   int convertToCPLEX(Status s){

     switch(s){
       case LPSolver::BASIC:
         return CPX_BASIC;
       case LPSolver::UPPER:
         return CPX_AT_UPPER;
       case LPSolver::LOWER:
         return CPX_AT_LOWER;
       case LPSolver::SUPERBASIC:
       case LPSolver::FREE:
         return CPX_FREE_SUPER;
       case LPSolver::FIXED:
       case LPSolver::INF:
       case LPSolver::UNKNOWN:
       case LPSolver::END:
         return -1;
     }
     return -1;

   };



};


#endif
