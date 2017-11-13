

#ifndef LEMONSOLVER_H
#define LEMONSOLVER_H


#include <lemon/smart_graph.h>
#include <lemon/network_simplex.h>


class LemonSolver : public LPSolver{

  private:
    typedef LPSolver::Status Status;
     DIGRAPH_TYPEDEFS(SmartDigraph);
     lemon::SmartGraph graph;
     leom::NetworkSimplex *solver;



  public:

    LemonSolver() {
      graph = NULL;
      solver = NULL;
    };

    ~LemonSolver(){
      deleteLP();
      if(solver != NULL){
        delete solver;
      }
    };




  

   virtual void createLP(int nSource, int nTarget){
     deleteLP();
     graph->reserveNode(nSource + nTarget);
   };



   virtual void solveLP(){
     solver = new NetworkSimplex( &graph );
   };



   virtual void deleteLP(){
     graph.clear();
   };



   virtual void addCols(int n){
     graph.reserveEdge( n );
   };



   virtual void addRows(int n){
     for(int i=0;i<n; i++){
       g.addNode();
     }
   };
 


   virtual double getRowDual(int row){
     return solver->potential( graph.nodeFromId(row) );
   };



   virtual double getColPrimal(int col){
     return solver->flow( graph.edgeFromId(col) );
   };
   


   virtual void setRowBounds(int i, double mass){
     
   };

   virtual double getRowBounds(int i){
     
   };

   virtual void setColBounds(int i){
   
   };
   
   virtual void setCoefficent(int i, double cost){
   };
   
   virtual void setColConstraints(int col, int s, int t){
   };

   virtual int getColConstraints(int col, int *ind, double *val){
   };
   
   virtual Status getColStatus(int col){
   };
   
   virtual Status getRowStatus(int row){
   };
   
   virtual void setColStatus(int col, Status s){
   };
   
   virtual void setRowStatus(int row, Status s){
   };

   virtual double getObjectiveValue(){
   };


   virtual int getIterationCount(){
   };


   virtual int getNumRows(){

   };

   virtual int getNumCols(){

   };
   
   virtual void setupStdBasis(){

   };



  private:

   Status convertFromLemon( s ){
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


   int convertToLemon(Status s){
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
