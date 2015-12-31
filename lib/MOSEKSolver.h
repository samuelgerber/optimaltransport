

#ifndef MOSEKSOLVER_H
#define MOSEKSOLVER_H


#include "mosek.h"



static void MSKAPI printstr(void *handle, 
                            MSKCONST char str[]) 
{ 
  printf("%s",str); 
} /* printstr */ 

class MOSEKSolver : public LPSolver{

  private:
    typedef LPSolver::Status Status;
  
    MSKenv_t     env; 
    MSKtask_t    task; 
    double *primal;
    double *dual;
    MSKint32t optimType;
    std::vector< MSKstakeye > colStatus;
    std::vector< MSKstakeye > rowStatus;
    


  public:

    MOSEKSolver(MSKint32t optimizer = MSK_OPTIMIZER_NETWORK_PRIMAL_SIMPLEX ) {
        MSK_makeenv(&env,NULL); 
        task = NULL;
        primal = NULL;
        dual=NULL;
        optimType=optimizer;
    };

    ~MOSEKSolver(){
      deleteLP();
      MSK_deleteenv(&env);
    };




  

   virtual void createLP(int nSource, int nTarget){
     deleteLP();

     MSK_maketask(env, 0, 0,&task);
     MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
     MSK_putintparam(task, MSK_IPAR_OPTIMIZER,
     //    MSK_OPTIMIZER_FREE_SIMPLEX);
        optimType); 
     MSK_putintparam(task, MSK_IPAR_SIM_HOTSTART, MSK_SIM_HOTSTART_STATUS_KEYS);
     MSK_putintparam(task, MSK_IPAR_PRESOLVE_USE, MSK_PRESOLVE_MODE_OFF);
     //MSK_putintparam(task, MSK_IPAR_PRESOLVE_ELIMINATOR_USE, MSK_PRESOLVE_MODE_OFF);
     MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL, printstr);  
   };



   virtual void solveLP(){

     if(primal != NULL){
       delete[] primal;
     }
     if(dual != NULL){
       delete[] dual;
     }
     primal = new double[ getNumCols() ];
     dual = new double[ getNumRows() ];

     
     MSK_putskx (
         task,
         MSK_SOL_BAS,
         colStatus.data() );

     MSK_putskc (
         task,
         MSK_SOL_BAS,
         rowStatus.data() );
         

     MSKrescodee res = MSK_optimize(task);
     std::cout << "MSK result code: " << res << std::endl;

     MSK_getxx(task, 
         MSK_SOL_BAS,    /* Request the basic solution. */ 
         primal);

     MSK_gety(task, 
         MSK_SOL_BAS,    /* Request the basic solution. */ 
         dual);

     MSK_getskx(task, MSK_SOL_BAS, colStatus.data() );  
     MSK_getskc(task, MSK_SOL_BAS, rowStatus.data() );  
   };



   virtual void deleteLP(){
     if(task != NULL){
       MSK_deletetask(&task);
       task = NULL;
       delete[] primal;
       delete[] dual;
       colStatus.clear();
       rowStatus.clear();
       dual = NULL;
       primal = NULL;
     }
   };



   virtual void addCols(int n){
     MSK_appendvars(task,n);
     for(int i=0; i<n; i++){
       colStatus.push_back( MSK_SK_LOW );
     }
   };



   virtual void addRows(int n){
     MSK_appendcons(task, n); 
     for(int i=0; i<n; i++){
       rowStatus.push_back( MSK_SK_LOW );
     }
   };
 


   virtual double getRowDual(int row){
     return dual[row];
   };



   virtual double getColPrimal(int col){
     return primal[col];
   };
   


   virtual void setRowBounds(int i, double mass){

     MSK_putconbound(task, 
                          i,           /* Index of constraint.*/ 
                          MSK_BK_FX,      /* Bound key.*/ 
                          mass,      /* Numerical value of lower bound.*/ 
                          mass);     /* Numerical value of upper bound.*/ 
   };

   virtual double getRowBounds(int i){
     MSKboundkeye bk;
     MSKrealt bl;
     MSKrealt bu;
     MSK_getbound (
         task,
         MSK_ACC_CON,
         i,
         &bk,
         &bl,
         &bu);
     return bu;
     
   };

   virtual void setColBounds(int i){
     MSK_putvarbound(task, 
         i,           /* Index of variable.*/ 
         MSK_BK_LO,      /* Bound key.*/ 
         0,      /* Numerical value of lower bound.*/ 
         0);     /* Numerical value of upper bound.*/ 
   };
   
   virtual void setCoefficent(int i, double cost){
     MSK_putcj(task, i, cost);
   };
   
   virtual void setColConstraints(int col, int s, int t){
     int ind[2] = {s,t};
     double val[2] = {1,-1};
     MSK_putacol(task, 
                        col,                 /* Variable (column) index.*/ 
                        2, /* Number of non-zeros in column j.*/ 
                        ind,     /* Pointer to row indexes of column j.*/ 
                        val);    /* Pointer to Values of column j.*/ 
   };

   virtual int getColConstraints(int col, int *ind, double *val){
     int n = 0;
     MSK_getacol(task, col, &n, ind, val);
     return n;    
   };
   
   virtual Status getColStatus(int col){
     return convertFromMOSEK( colStatus[col] );
   };
   
   virtual Status getRowStatus(int row){
     return convertFromMOSEK( rowStatus[row] );
   };
   
   virtual void setColStatus(int col, Status s){
     colStatus[col] = convertToMOSEK(s);
   };
   
   virtual void setRowStatus(int row, Status s){
     rowStatus[row] = convertToMOSEK(s);
   };

   virtual double getObjectiveValue(){
     double obj;
     MSK_getprimalobj(task, MSK_SOL_BAS, &obj);
     return obj;
   };


   virtual int getIterationCount(){
     static int n=16;
     static MSKiinfiteme iters[]= {
       MSK_IINF_INTPNT_ITER,
       MSK_IINF_SIM_DUAL_DEG_ITER,
       MSK_IINF_SIM_DUAL_INF_ITER,
       MSK_IINF_SIM_DUAL_ITER,
       MSK_IINF_SIM_NETWORK_DUAL_DEG_ITER,
       MSK_IINF_SIM_NETWORK_DUAL_INF_ITER,
       MSK_IINF_SIM_NETWORK_DUAL_ITER,
       MSK_IINF_SIM_NETWORK_PRIMAL_DEG_ITER,
       MSK_IINF_SIM_NETWORK_PRIMAL_INF_ITER,
       MSK_IINF_SIM_NETWORK_PRIMAL_ITER,
       MSK_IINF_SIM_PRIMAL_DEG_ITER,
       MSK_IINF_SIM_PRIMAL_DUAL_DEG_ITER,
       MSK_IINF_SIM_PRIMAL_DUAL_INF_ITER,
       MSK_IINF_SIM_PRIMAL_DUAL_ITER,
       MSK_IINF_SIM_PRIMAL_INF_ITER,
       MSK_IINF_SIM_PRIMAL_ITER
     };

     int nIter = 0;

     for(int i=0; i<n; i++){
       MSKint32t count;

       MSK_getintinf (
           task,
           iters[i],
           &count
           );
       nIter += count;
     }
     return nIter;

   };


   virtual int getNumRows(){
     int numcon;
     MSK_getnumcon (
         task,
         &numcon);
     return numcon;
   };

   virtual int getNumCols(){
     int numvar;
     MSK_getnumvar (
         task,
         &numvar);
     return numvar;
   };
   
   virtual void setupStdBasis(){
     for(int i= 0; i< getNumCols(); i++){
       colStatus[i] = MSK_SK_LOW;
     }
     for(int i= 0; i< getNumRows(); i++){
       rowStatus[i] = MSK_SK_BAS;
     }
   };



  private:

   Status convertFromMOSEK(MSKstakeye s){
     switch(s){
       case MSK_SK_BAS:
         return LPSolver::BASIC;
       case MSK_SK_UPR:
         return LPSolver::UPPER;
       case MSK_SK_LOW:
         return LPSolver::LOWER;
       case MSK_SK_UNK:
         return LPSolver::UNKNOWN;
       case MSK_SK_SUPBAS:
         return LPSolver::SUPERBASIC;
       case MSK_SK_FIX:
         return LPSolver::FIXED;
       case MSK_SK_INF:
         return LPSolver::INF;
       case MSK_SK_END:
         return LPSolver::END; 
     }
     return LPSolver::END; 
   };


   MSKstakeye convertToMOSEK(Status s){
     switch(s){
       case LPSolver::BASIC:
         return MSK_SK_BAS;
       case LPSolver::UPPER:
         return MSK_SK_UPR;
       case LPSolver::LOWER:
         return MSK_SK_LOW;
       case LPSolver::UNKNOWN:
       case LPSolver::FREE:
         return MSK_SK_UNK;
       case LPSolver::SUPERBASIC:
         return MSK_SK_SUPBAS;
       case LPSolver::FIXED:
         return MSK_SK_FIX;
       case LPSolver::INF:
         return MSK_SK_INF;
       case LPSolver::END:
         return MSK_SK_END;

     }
     return MSK_SK_END;
   };



};


#endif
