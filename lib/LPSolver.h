#ifndef LPSOLVER_H
#define LPSOLVER_H




class LPSolver {

  
  public:

    
   enum Status {UPPER, LOWER, BASIC, UNKNOWN, SUPERBASIC, FIXED, INF, END, FREE};

   virtual ~LPSolver(){};

   virtual void createLP(int nSource, int nTarget) = 0;
   virtual void solveLP() = 0;
   virtual bool isOptimal() = 0;
   
   virtual int getNumRows() = 0;
   virtual int getNumCols() = 0;
   virtual double getObjectiveValue() = 0;
   virtual int getIterationCount() = 0;
   virtual void addCols(int n) = 0;
   virtual void addRows(int n) = 0;
   virtual double getRowDual(int row) = 0;
   virtual double getColPrimal(int col) = 0;
   virtual void setRowBounds(int i, double mass) = 0;
   virtual double getRowBounds(int i) = 0;
   virtual void setColBounds(int i) = 0;
   virtual void setColBounds(int col, double lb, double ub) = 0;
   virtual void setCoefficent(int i, double cost) = 0;
   virtual void setColConstraints(int col, int s, int t) = 0;
   virtual int getColConstraints(int col, int *ind, double *val) = 0;
   virtual Status getColStatus(int col) = 0;
   virtual Status getRowStatus(int row) = 0;
   virtual void setColStatus(int col, Status s) = 0;
   virtual void setRowStatus(int row, Status s) = 0;
   virtual void setupStdBasis() = 0;



};


#endif


