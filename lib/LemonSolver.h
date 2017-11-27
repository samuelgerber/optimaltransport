

#ifndef LEMONSOLVER_H
#define LEMONSOLVER_H

#include "LPSolver.h"

#include <lemon/smart_graph.h>
#include <lemon/network_simplex.h>


class LemonSolver : public LPSolver{

  private:
    typedef LPSolver::Status Status;

    std::vector<int> sInd;
    std::vector<int> tInd;

    std::vector<double> coeff;
    std::vector<double> mass;
    std::vector<double> primal;
    std::vector<double> dual;
    std::vector<double> colLB;
    std::vector<double> colUB;

    std::vector<Status> colStatus;
    std::vector<Status> rowStatus;

    int solstat;
    int iCount;
    double objValue;
    double lambda;
    int ns;
    int nt;
    bool success;




  public:

    LemonSolver(double l=-1) {
      lambda = l;
      ns = 0;
      nt = 0;
      success = false;
    };

    ~LemonSolver(){
      deleteLP();
    };





   virtual void createLP(int nSource, int nTarget){
     deleteLP();
     ns=nSource;
     nt=nTarget;

     if(lambda > 0){
       rowStatus.resize( 1, LPSolver::BASIC );
       dual.resize( 1, 0);

       primal.resize( nt, 0);
       colStatus.resize( nt, LPSolver::LOWER );
     }

   };



   virtual void solveLP(){

     using namespace lemon;

     DIGRAPH_TYPEDEFS(SmartDigraph);
     SmartDigraph graph;
     int status;

     //Lemon allwos integer only
     static const int maxVal = 10000;
     double costScaling = 0;
     for(int i=0; i<coeff.size(); i++){
       if(coeff[i] > costScaling){
         costScaling = coeff[i];
       }
     }
     if( costScaling == 0){
       costScaling = maxVal;
     }
     else{
       costScaling = ((double) maxVal ) / costScaling;
     }

     double capacityScaling = 0;
     for(int i=0; i<mass.size(); i++){
        double tmp = fabs(mass[i]);
        if( tmp > capacityScaling){
          capacityScaling = tmp;
        }
     }
     capacityScaling = ((double) maxVal ) / capacityScaling;

     std::cout << "capacityScaling: "<< capacityScaling << std::endl;
     std::cout << "costScaling: "<< costScaling << std::endl;

     //Setup nodes
     int massPositive = 0;
     int massNegative = 0;
     graph.reserveNode( mass.size() +1 );
     for(int i=0; i <= mass.size(); i++){
        graph.addNode();
     }
     IntNodeMap supply(graph);
     for(int i=0; i<mass.size(); i++){
        int m = (int) ( mass[i] * capacityScaling );
        supply[ graph.nodeFromId(i) ] = m;
        if( m > 0 ){
          massPositive += m;
        }
        else{
          massNegative += m;
        }
     }

     int massImbalance = massPositive + massNegative;
     Node terminal = graph.nodeFromId( mass.size() );
     if( lambda > 0){
       for(int i=ns; i<mass.size(); i++){
         supply[ graph.nodeFromId(i) ] = 0;
       }
       supply[ terminal ] = -massPositive;
     }
     else{
       std::cout << "Mass positive: " << massPositive << std::endl;
       std::cout << "Mass negative: " << massNegative << std::endl;
       std::cout << "Mass imbalance: " << massImbalance << std::endl;
       supply[ terminal ] = -massImbalance;
     }




     //Setup arcs
     if( lambda > 0 || massImbalance < 0 ){
       graph.reserveArc( sInd.size() + nt );
     }
     else{
       graph.reserveArc( sInd.size() + ns );
     }
     //Add regular node
     for( int i=0; i < sInd.size(); i++){
       graph.addArc( graph.nodeFromId(sInd[i]), graph.nodeFromId(tInd[i]) );
     }

     //Add arcs to terminal node to deal with fuzzy match or imbalances
     if(lambda > 0 || massImbalance < 0){
       for(int i=0; i < nt; i++){
         graph.addArc( graph.nodeFromId(i+ns), terminal );
       }
     }
     else{
       for(int i=0; i < ns; i++){
         graph.addArc( graph.nodeFromId(i), terminal );
       }
     }


     IntArcMap capacity(graph), lower(graph), cost(graph);
     for( int i=0; i < coeff.size(); i++){
       Arc a = graph.arcFromId(i);
       capacity[a] = (int) (capacityScaling * colUB[i] ) + 1;
       lower[a] = (int) std::max(0, (int) (capacityScaling * colLB[i] ) -1 );
       cost[a] = (int)( costScaling * coeff[i]);
     }

     if( lambda > 0 ){
       for(int i=coeff.size(); i<coeff.size()+nt; i++){
         Arc a = graph.arcFromId(i);
         double m = fabs(mass[i-coeff.size()+ns]);
         double delta = m*lambda;
         capacity[a] = (int) (capacityScaling * m+delta ) + 1;
         lower[a] = (int) std::max(0, (int) (capacityScaling * m-delta ) - 1);
         cost[a] = 0;
       }
     }
     else{
       if( massImbalance < 0 ){
         for(int i=coeff.size(); i<coeff.size()+nt; i++){
           Arc a = graph.arcFromId(i);
           capacity[a] = capacityScaling;
           lower[a] = 0;
           cost[a] = 0;
         }
       }
       else{
         for(int i=coeff.size(); i<coeff.size()+ns; i++){
           Arc a = graph.arcFromId(i);
           capacity[a] = capacityScaling;
           lower[a] = 0;
           cost[a] = 0;
         }
       }
     }


     //Solve
     NetworkSimplex<SmartDigraph> ns(graph);
     ns.upperMap(capacity);
     ns.lowerMap(lower);
     ns.costMap(cost);
     ns.supplyMap(supply);


     NetworkSimplex<SmartDigraph>::ProblemType res = ns.run();
     success = res == NetworkSimplex<SmartDigraph>::OPTIMAL;
     objValue  = ns.totalCost<double>();
     objValue /= capacityScaling;
     objValue /= costScaling;
     iCount = 1;

     std::cout << "success: " << success << std::endl;
     std::cout << "objective: " << objValue << std::endl;


     std::cout << "primal.size(): " << primal.size() << std::endl;
     for(int i=0; i< primal.size(); i++){
       primal[i] = ( (double) ns.flow( graph.arcFromId( i ) ) ) / capacityScaling;
     }
     std::cout << "dual.size(): " << dual.size() << std::endl;
     for(int i=0; i<dual.size(); i++){
       dual[i] = ns.potential( graph.nodeFromId(i) );
     }
     std::cout << "stored "<< std::endl;


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
     return success;
   };

   virtual void addCols(int n){
     sInd.resize( sInd.size() + n, -1  );
     tInd.resize( tInd.size() + n, -1 );
     coeff.resize( coeff.size() + n, 0 );
     primal.resize( primal.size() + n, 0 );
     colStatus.resize( colStatus.size() + n, LPSolver::LOWER );
     colLB.resize( colLB.size() + n, 0 );
     colUB.resize( colUB.size() + n, 1 );
   };



   virtual void addRows(int n){
     mass.resize( mass.size() + n , 0);
     rowStatus.resize( rowStatus.size() + n, LPSolver::BASIC );
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
     setColBounds(i, 0, 1);
   };

   virtual void setColBounds(int col, double lb, double ub){
     colLB[col] = lb;
     colUB[col] = ub;
   };

   virtual Status getColStatus(int col){
     return  colStatus[col];

   };

   virtual Status getRowStatus(int row){
     return rowStatus[row];
   };

   virtual void setColStatus(int col, Status s){
     colStatus[col] = s;
   };


   virtual void setRowStatus(int row, Status s){
     rowStatus[row] = s;
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
       colStatus[i] = LPSolver::LOWER;
     }
     for(int i= 0; i< getNumRows(); i++){
       rowStatus[i] =  LPSolver::BASIC;
     }
   };

};


#endif
