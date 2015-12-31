

#ifndef CPLEXSOLVER_H
#define CPLEXSOLVER_H

#include "LPSolver.h"

#include <ilcplex/cplex.h>
#include <stdlib.h>



class CPLEXSolver : public LPSolver{

  private:
    typedef LPSolver::Status Status;
 

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

    CPLEXSolver(int optimizer = CPX_ALG_AUTOMATIC, double l = -1 ){

      int status=-1;
      env = CPXopenCPLEX (&status);
      CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
      CPXsetdblparam (env, CPX_PARAM_NETEPOPT, 1e-10);
      CPXsetintparam (env, CPX_PARAM_LPMETHOD, optimizer);

      lambda = l;
      solstat = CPX_STAT_ABORT_USER;
    };

    ~CPLEXSolver(){
      deleteLP();
      CPXcloseCPLEX (&env);
    };


   void setLambda(double l){
     lambda = l;  
   };


   virtual void createLP(int nSource, int nTarget){
     deleteLP();
     ns=nSource;
     nt=nTarget;
     
     if(lambda > 0){
       rowStatus.resize( 1, CPX_BASIC );
       dual.resize( 1, 0);
     }

   };


   virtual void solveLP(){

     int status;
     CPXLPptr prob = CPXcreateprob (env, &status, "ot");
   
    /* 
     if(lambda > 0 ){
       std::vector<double> ctmp = coeff;
       for(unsigned int i=0; i<ctmp.size(); i++){
         ctmp[i] += lambda * mass[ tInd[i] ];
       }
       CPXnewcols(env, prob, ctmp.size(), ctmp.data(), NULL, NULL, NULL, NULL);
     }
     else{
     */
       CPXnewcols(env, prob, coeff.size(), coeff.data(), colLB.data(), colUB.data(), NULL, NULL);
     //}

     if(lambda > 0){
       std::vector<double> mtmp = mass;
       std::vector<double> range(mtmp.size(), 0);
       std::vector<char> sense(mtmp.size(), 'E');
       int n= ns;
       for(unsigned int i=n; i<mtmp.size(); i++){
         double b = mtmp[i];
         double r = lambda*b;
         mtmp[i] = std::min(b-r, 0.0);
         range[i] = 2*r;
         sense[i] = 'R';
       }
       //probability function constraint
       mtmp.push_back(1);
       sense.push_back('E');
       range.push_back(0);
       CPXnewrows(env, prob, mtmp.size(), mtmp.data(), sense.data(), range.data(), NULL);
     }else{
       CPXnewrows(env, prob, mass.size(), mass.data(), NULL, NULL, NULL);
     }

     std::vector<int> colIndex(tInd.size());
     for(unsigned int i=0; i< colIndex.size(); i++){
       colIndex[i] = i;
     }
     std::vector<double> ones(sInd.size(), 1);
     std::vector<double> nones(tInd.size(), -1);

     CPXchgcoeflist(env, prob, colIndex.size(), sInd.data(), colIndex.data(),
         ones.data());
     CPXchgcoeflist(env, prob, colIndex.size(), tInd.data(), colIndex.data(),
           nones.data());

     if(lambda>0){
       std::vector<int> row( colIndex.size(), mass.size() );
       CPXchgcoeflist(env, prob, colIndex.size(), row.data(), colIndex.data(),
           ones.data());
     }

     //if(lambda <= 0){
       CPXcopybase (env, prob, colStatus.data(), rowStatus.data());
     //}

     /*
     if(lambda > 0){
       std::map<int, std::vector<int> > marginal;
       for(unsigned int i=0; i<tInd.size(); i++){
         int row = tInd[i];
         marginal[row].push_back(i);
       }

       std::map<int, std::vector<int> >::iterator it =marginal.begin();
       for(;it != marginal.end(); ++it){
         std::vector<int> &m = it->second;
         for(unsigned int i=0; i<m.size(); i++){
           CPXchgqpcoef(env, prob, m[i], m[i], lambda);
           for(unsigned int j=(i+1); j<m.size(); j++){
             CPXchgqpcoef(env, prob, m[i], m[j], lambda/2.0);
           }
         }
       } 
      
     }
*/

  //   if(lambda > 0 ){
  //     status = CPXqpopt (env, prob);
  //   }
  //   else{
       status = CPXlpopt (env, prob);
  //   }


     status = CPXsolution (env, prob, &solstat, &objValue, primal.data(),
         dual.data(), NULL, NULL);

     std::cout << "CPX status: " << status << std::endl;
     std::cout << "CPX solution status: " << solstat << std::endl;

     iCount = CPXgetitcnt(env, prob);
    
    // if(lambda  <= 0){ 
       CPXgetbase (env, prob, colStatus.data(), rowStatus.data());
    // }
    
    /* 
     for(int i=0; i<colStatus.size(); i++){
       std::cout << colStatus[i] << ", ";
     }
     std::cout << std::endl;
     std::cout << std::endl;
     */

     CPXfreeprob(env, &prob);

     std::cout << "Obj: " << objValue << std::endl;
     std::cout << "#iter: " <<iCount << std::endl << std::endl;
   };
   
   virtual bool isOptimal(){
     return solstat == CPX_STAT_OPTIMAL;
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
   };
   
   virtual void setCoefficent(int i, double cost){
     coeff[i] = cost;
   };
   
  

   virtual void setColConstraints(int col, int s, int t){
     ns = std::max(ns, s+1);
     tInd[col] = t;
     sInd[col] = s;
   };
   

   virtual void setColBounds(int col, double lb, double ub){
     colLB[col] = lb;
     colUB[col] = ub;
   };


   virtual int getColConstraints(int col, int *ind, double *val){
     ind[0] = sInd[col];
     ind[1] = tInd[col];
     val[0] = colLB[col];
     val[1] = colUB[col];
     return 2;    
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
       case LPSolver::END:
       case LPSolver::UNKNOWN:
         return -1;
     }
     return -1;

   };



};


#endif
