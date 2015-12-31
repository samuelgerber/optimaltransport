#ifndef GLPKSOLVER_H
#define GLPKSOLVER_H

#include "LPSolver.h"

#include <glpk.h>





typedef struct DMP DMP;
typedef struct GLPROW GLPROW;
typedef struct GLPCOL GLPCOL;
typedef struct AVL AVL;
typedef struct BFD BFD;

struct glp_prob
{     /* LP/MIP problem object */
      unsigned magic;
      /* magic value used for debugging */
      DMP *pool;
      /* memory pool to store problem object components */
      glp_tree *tree;
      /* pointer to the search tree; set by the MIP solver when this
         object is used in the tree as a core MIP object */
      void *parms;
      /* reserved for backward compatibility */
      /*--------------------------------------------------------------*/
      /* LP/MIP data */
      char *name;
      /* problem name (1 to 255 chars); NULL means no name is assigned
         to the problem */
      char *obj;
      /* objective function name (1 to 255 chars); NULL means no name
         is assigned to the objective function */
      int dir;
      /* optimization direction flag (objective "sense"):
         GLP_MIN - minimization
         GLP_MAX - maximization */
      double c0;
      /* constant term of the objective function ("shift") */
      int m_max;
      /* length of the array of rows (enlarged automatically) */
      int n_max;
      /* length of the array of columns (enlarged automatically) */
      int m;
      /* number of rows, 0 <= m <= m_max */
      int n;
      /* number of columns, 0 <= n <= n_max */
      int nnz;
      /* number of non-zero constraint coefficients, nnz >= 0 */
      GLPROW **row; /* GLPROW *row[1+m_max]; */
      /* row[i], 1 <= i <= m, is a pointer to i-th row */
      GLPCOL **col; /* GLPCOL *col[1+n_max]; */
      /* col[j], 1 <= j <= n, is a pointer to j-th column */
      AVL *r_tree;
      /* row index to find rows by their names; NULL means this index
         does not exist */
      AVL *c_tree;
      /* column index to find columns by their names; NULL means this
         index does not exist */
      /*--------------------------------------------------------------*/
      /* basis factorization (LP) */
      int valid;
      /* the factorization is valid only if this flag is set */
      int *head; /* int head[1+m_max]; */
      /* basis header (valid only if the factorization is valid);
         head[i] = k is the ordinal number of auxiliary (1 <= k <= m)
         or structural (m+1 <= k <= m+n) variable which corresponds to
         i-th basic variable xB[i], 1 <= i <= m */
      glp_bfcp *bfcp;
      /* basis factorization control parameters; may be NULL */
      BFD *bfd; /* BFD bfd[1:m,1:m]; */
      /* basis factorization driver; may be NULL */
      /*--------------------------------------------------------------*/
      /* basic solution (LP) */
      int pbs_stat;
      /* primal basic solution status:
         GLP_UNDEF  - primal solution is undefined
         GLP_FEAS   - primal solution is feasible
         GLP_INFEAS - primal solution is infeasible
         GLP_NOFEAS - no primal feasible solution exists */
      int dbs_stat;
      /* dual basic solution status:
         GLP_UNDEF  - dual solution is undefined
         GLP_FEAS   - dual solution is feasible
         GLP_INFEAS - dual solution is infeasible
         GLP_NOFEAS - no dual feasible solution exists */
      double obj_val;
      /* objective function value */
      int it_cnt;
      /* simplex method iteration count; increased by one on performing
         one simplex iteration */
      int some;
      /* ordinal number of some auxiliary or structural variable having
         certain property, 0 <= some <= m+n */
      /*--------------------------------------------------------------*/
      /* interior-point solution (LP) */
      int ipt_stat;
      /* interior-point solution status:
         GLP_UNDEF  - interior solution is undefined
         GLP_OPT    - interior solution is optimal
         GLP_INFEAS - interior solution is infeasible
         GLP_NOFEAS - no feasible solution exists */
      double ipt_obj;
      /* objective function value */
      /*--------------------------------------------------------------*/
      /* integer solution (MIP) */
      int mip_stat;
      /* integer solution status:
         GLP_UNDEF  - integer solution is undefined
         GLP_OPT    - integer solution is optimal
         GLP_FEAS   - integer solution is feasible
         GLP_NOFEAS - no integer solution exists */
      double mip_obj;
      /* objective function value */
};

int get_it_cnt(glp_prob *P) { return P->it_cnt; };


class GLPKSolver : public LPSolver{

  private:
    typedef LPSolver::Status Status;
    glp_prob *lp;

    
  public:
 
    GLPKSolver(){
      lp = NULL;
    };

   ~GLPKSolver(){
      deleteLP();
   };

   virtual void createLP(int nSource, int nTarget){
      deleteLP();
      lp = glp_create_prob();
      glp_set_obj_dir(lp, GLP_MIN);
   };

   virtual void solveLP(){
      deleteLP()p;
      glp_smcp parm;
      glp_init_smcp(&parm);
      int ret = glp_simplex(lp, &parm);
   };

   virtual void deleteLP(){
     if(lp !=NULL){
       glp_delete_prob(lp);
       lp = NULL;
     }
   };

   virtual void addCols(int n){
      int start = glp_add_cols(lp, n);
      std::cout << "start: " << start << std::endl;
   };

   virtual void addRows(int n){
      glp_add_rows(lp, n);
   };
   
   virtual double getRowDual(int row){
     return glp_get_row_dual(lp, row+1);
   };

   virtual double getColPrimal(int col){
     return glp_get_col_prim(lp, col+1);
   };
   
   virtual void setRowBounds(int i, double mass){
     glp_set_row_bnds(lp, i+1, GLP_FX, mass, mass); 
   };

   virtual double getRowBounds(int i){
     return glp_get_row_ub(lp, i+1); 
   };

   virtual void setColBounds(int i){
     glp_set_col_bnds(lp, i+1, GLP_LO, 0, 0);
   };
   
   virtual void setCoefficent(int i, double cost){
     glp_set_obj_coef(lp, i+1, cost);
   };
   
   virtual void setColConstraints(int col, int s, int t){
     int ind[3] = {0,s+1,t+1};
     double val[3] = {1,1,-1};
     glp_set_mat_col(lp, col+1, 2, ind, val); 
   };

   virtual int getColConstraints(int col, int *ind, double *val){
     int i[3];
     double v[3];
     int n = glp_get_mat_col(lp, col+1, i, v); 
     ind[0] = i[1]-1;
     ind[1] = i[2]-1;
     val[0] = v[1];
     val[1] = v[2];
     return n;

   };
   
   virtual Status getColStatus(int col){
     int s = glp_get_col_stat(lp, col+1);
     return convertFromGPLK(s);
   };
   
   virtual Status getRowStatus(int row){
     int s = glp_get_row_stat(lp, row+1);
     return convertFromGPLK(s);
   };
   
   virtual void setColStatus(int col, Status s){
     int st = convertToGPLK(s);
     glp_set_col_stat(lp, col+1, st);
   };
   
   virtual void setRowStatus(int row, Status s){
     int st = convertToGPLK(s);
     glp_set_row_stat(lp, row+1, st);
   };

   virtual double getObjectiveValue(){
     return glp_get_obj_val(lp);
   };

   virtual int getIterationCount(){
     return get_it_cnt( lp );
   };

   virtual void initLP(){
   };

   virtual int getNumRows(){
     return glp_get_num_rows(lp);
   };

   virtual int getNumCols(){
     return glp_get_num_cols(lp);
   };
   
   virtual void setupStdBasis(){
     glp_std_basis(lp);
   };



  private:

   Status convertFromGPLK(int s){
     switch(s){
       case GLP_BS:
         return LPSolver::BASIC;
       case GLP_NU:
         return LPSolver::UPPER;
       case GLP_NL:
         return LPSolver::LOWER;
       case GLP_NS:
         return LPSolver::FIXED;
       case GLP_NF:
         return LPSolver::FREE;

     }
     return LPSolver::END; 
   };


   int convertToGPLK(Status s){
     switch(s){
       case LPSolver::BASIC:
         return GLP_BS;
       case LPSolver::UPPER:
         return GLP_NU;
       case LPSolver::LOWER:
         return GLP_NL;
       case LPSolver::FREE:
         return GLP_NF;
       case LPSolver::FIXED:
         return GLP_NS;

     }
     return GLP_NF;
   };



};


#endif
