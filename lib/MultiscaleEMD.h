#ifndef MULTISCALEEMD_H
#define MULTISCALEEMD_H

#include "Linalg.h"
#include "LinalgIO.h"
#include "DenseVector.h"
#include "DenseMatrix.h"

#include "Eigen/Sparse"

#include <glpk.h>

#include <set>
#include <map>
#include <list>
#include <set>
#include <vector>

typedef std::map< std::pair<int,int>, double> PMap;
typedef std::map< std::pair<int, int>, std::set< std::pair<int, int> > > VMap;

//typedef std::pair< VMAP::iterator, bool > VMapInsert;
 
typedef std::map< std::pair<int, int>, PMap > DMap;

typedef std::pair< std::pair<int, int>, std::pair<int, int> > Path; 

typedef std::map< int, Path > IMap;
typedef std::map< Path, int > PathMap;

typedef std::map<std::pair<int, int>, int> RMap;




template <typename TPrecision>
class ScaleHistogram{
public:


  ScaleHistogram(int m, int n){          
    p = DenseMatrix<TPrecision>(m, n);
    Linalg<TPrecision>::Zero(p);
    mI = DenseMatrix<TPrecision>(m, n);
    Linalg<TPrecision>::Zero(mI);
    mJ = DenseMatrix<TPrecision>(m, n);
    Linalg<TPrecision>::Zero(mJ);
  };


  DenseMatrix<TPrecision> p;
  DenseMatrix<TPrecision> mI;
  DenseMatrix<TPrecision> mJ;

  int radius;
};





template <typename TPrecision>
class EMD{
public:
 
  ScaleHistogram<TPrecision> *sourceHist;
  ScaleHistogram<TPrecision> *targetHist;

  //Mass sent from source to which targets
  VMap source;
  //Mass received from sources at  target
  VMap target;

  //Distances
  DMap D;
  //Weights (i.e. transport map)
  DMap W;

  //
  VMap neighbors;

  //Dual solution
  DMap delta;
  PMap phi;
  PMap psi;

  double cost;
  double optimizationStatus;
  int nVariables;
};




template <typename TPrecision>
class MultiscaleEMD{

  public:

    enum Neighborhood {BALL, DUAL, POLY};

    static std::list< EMD<TPrecision> *> solve(DenseMatrix<TPrecision> A,
        DenseMatrix<TPrecision> B, Neighborhood neighborhood = BALL, int p = 1,
        int radius = 1, int nScales = -1){

      //compute distances at each tree level
      typename std::list< EMD<TPrecision> *> solutions; 

      int mA = A.M();
      int nA = A.N();
      
      int mB = B.M();
      int nB = B.N();
      int splits = 2;


      std::list< ScaleHistogram<TPrecision> * > pyramidA = buildPyramid(A);
      std::list< ScaleHistogram<TPrecision> * > pyramidB = buildPyramid(B);

      int nScalesMax = std::max(pyramidA.size(), pyramidB.size());
      if(nScales == -1 || nScales >= nScalesMax){
        nScales = nScalesMax-1;
      }

      typename std::list< ScaleHistogram<TPrecision> * >::iterator itA =
        pyramidA.begin();
      typename std::list< ScaleHistogram<TPrecision> * >::iterator itB =
        pyramidB.begin();


      ScaleHistogram<TPrecision> *As = *itA;
      ScaleHistogram<TPrecision> *Bs = *itB;
      int scaleStart = nScalesMax-1-nScales;

      std::cout << scaleStart << std::endl;

      for(int i = 0; i < scaleStart; i++){
        if(itA != pyramidA.end()){
          ++itA;
          As = *itA;
        }
        if(itB != pyramidB.end()){
          ++itB;
          Bs = *itB;
        }
      }

      /*
      DMap pW;
      DMap pDelta;
      for(int ai=0; ai < As->mI.N()/2; ai++){
        for(int aj=0; aj < As->mJ.M()/2; aj++){
          std::pair<int, int> pA( ai, aj );
          for(int bi=0; bi < Bs->mI.N()/2; bi++){
            for(int bj=0; bj < Bs->mJ.M()/2; bj++){
              std::pair<int, int> pB(bi, bj );
              pW[pA][pB] = 1;
              pDelta[pA][pB] = 0;
            }
          }
        }
      }
*/

      EMD<TPrecision> *prevEMD = NULL;
      for(int n = scaleStart; n < nScalesMax; n++){
	      if(itA != pyramidA.end()){
          As = *itA;
          ++itA;
        }
        if(itB != pyramidB.end()){
          Bs = *itB;
          ++itB;
        }


        //Solve transport linaer program with hierarchical constraints and
	      //variable reductions
        EMD<TPrecision> *emd = solveLP(As, Bs, prevEMD, neighborhood, p, radius);
        
        /*
        if(prevEMD == NULL){
	        emd = solveLP(As, Bs, pW, pDelta, p, radius);
        }
        else{
	        emd = solveLP(As, Bs, prevEMD->W, prevEMD->delta, p, radius);
        }
        */


        solutions.push_back(emd);
        
        prevEMD = emd;

        


      }
      
      return solutions;
    };




  private:


    static EMD<TPrecision> *solveLP(ScaleHistogram<TPrecision> *source,
        ScaleHistogram<TPrecision> *target, EMD<TPrecision> *prevSol, Neighborhood
        neighborhood, int p=1, int r = 1){
  
        EMD<TPrecision> *sol = new EMD<TPrecision>();


        //Number of variables in LP. I.e. number of combinations of indices in A
        //to indices in B
        int nCombinations = 0;

        if(prevSol == NULL){
          //Add all variables
          for(int ai=0; ai < source->mI.M(); ai++){
            for(int aj=0; aj < source->mI.N(); aj++){
              std::pair<int, int> pA(ai, aj);
              if(source->p(ai, aj) == 0) continue;

              for(int bi=0; bi < target->mI.M(); bi++){
                for(int bj=0; bj < target->mI.N(); bj++){
                  std::pair<int, int> pB(bi, bj);
                
                  if(target->p(bi, bj) == 0) continue;
              
                  double t1 = source->mI(ai, aj) - target->mI(bi, bj);
                  double t2 = source->mJ(ai, aj) - target->mJ(bi, bj);
                  sol->D[pA][pB] = sqrt(t1*t1 + t2*t2);
                  sol->source[pA].insert(pB);
                  sol->target[pB].insert(pA);
                  nCombinations++;
                }
              }
            }
          }
        }
        else if( neighborhood == BALL ){
          nCombinations = BallNeighborhood(source, target, sol, prevSol, r, p);
        }
        else if(neighborhood == DUAL){
          //Add based on previous scales solution
          nCombinations = DualNeighborhood(source, target, sol, prevSol, r, p);
        }
        else{
         nCombinations = PolyNeighborhood(source, target, sol, prevSol, r, p);

        }

        sol->nVariables = nCombinations;


        //setup linear programm
        glp_prob *lp = glp_create_prob();
        glp_set_obj_dir(lp, GLP_MIN);
        
        //rows - constraints:
        // M of type sum of weights going from node A tree1 to nodes in tree2 =
        // weight of node A in tree1
        // N of type sum of weights from nodes in tree1 going to node A in tree2
        // = weight of node A in tree2
        int nConstraints = sol->source.size() + sol->target.size(); 
        glp_add_rows(lp, nConstraints);
        std::cout << "Sanity check: " << sol->source.size() << " x " << sol->target.size() << std::endl; 
        //columns - coefficents for each combination of nodes from tree1 to
        //tree2 = M*N
        glp_add_cols(lp, nCombinations );

        //set coefficients equal to cost of moving mass from node i in tree1 to
        //node j in tree2
        int index = 1;
        IMap imap;
        PathMap pmap;
        for(VMap::iterator it = sol->source.begin(); it != sol->source.end(); ++it){
          const std::pair<int, int> &from = it->first;
          std::set< std::pair<int, int> > &toList = it->second;
          for(std::set<std::pair<int, int> >::iterator toIt = toList.begin(); toIt!=toList.end(); ++toIt){
            glp_set_col_bnds(lp, index, GLP_LO, 0, 0);
            const std::pair<int, int> &to = *toIt;
            glp_set_obj_coef(lp, index, sol->D[from][to] ); 
            Path path(from, to);
            imap[index] = path;
            pmap[path] = index;
            ++index;
          }
        }

        //set constraints bounds

        int row=1;
        //The weights of each nodes in tree1
        RMap fromMap;
        for(VMap::iterator it = sol->source.begin(); it != sol->source.end(); ++it, ++row){
          const std::pair<int, int> &p = it->first;
          double w = source->p(p.first, p.second);
          glp_set_row_bnds(lp, row, GLP_FX, w, w);
          fromMap[p] = row;
        }


        //the weight of each node in tree2
        RMap toMap;
        for(VMap::iterator it = sol->target.begin(); it != sol->target.end(); ++it, ++row){
          const std::pair<int, int> &p = it->first;
          double w = target->p(p.first, p.second);
          glp_set_row_bnds(lp, row, GLP_FX, w, w);
          toMap[p] = row;
        }

        //Setup constraint matrix based on the following set of constraints:
        //1. the sum of the N weights leaving from A has to equal the weight in the
        //receiver in B 
        //2. the sum of the M weights ending up in node i of tree2 has to equal the
        //weight of node i in tree2
        //The constraint bounds above are the receing and sending amount of
        //mass. Each path needs to add a 1 in the corresponding rows (from and
        //to)
        int ind[3] = {0,0,0};
        double val[3] = {1,1,1};
        for(int i=0; i < nCombinations; i++){
          const Path &path = imap[i+1]; 
          const std::pair<int, int> &f1= path.first;
          const std::pair<int, int> &t1 = path.second;


          ind[1] = fromMap[f1];
          ind[2] = toMap[t1];

          fromMap[f1] = ind[1];
          toMap[t1] = ind[2];

          glp_set_mat_col(lp, i+1, 2, ind, val);
        }




        //
        //Set inital basic solution based on the solution at the previous
        //scale. Results in reducing the number of iterations required toi solve
        //the lienar program at this scale
        //
        if(prevSol != NULL){
        //if(false){
          glp_std_basis(lp);
          int nBaseCol = 0;
          
          DenseMatrix<TPrecision> targetP = Linalg<TPrecision>::Copy( target->p );
          DenseMatrix<TPrecision> sourceP = Linalg<TPrecision>::Copy( source->p );
          
          for(VMap::iterator it = prevSol->source.begin(); it !=
              prevSol->source.end(); ++it){
            const std::pair<int, int> &from = it->first;

            int fi = from.first*2;
            int fj = from.second*2;

            //number of variables within the partition sending mass

            int nFrom = 0;
            double sumFrom = 0;
            for(int i=0; i<2; i++){
              for(int j=0; j<2; j++){
                double wFrom = sourceP(fi+i, fj+j); 
                if(wFrom > 0){
                  nFrom ++;
                  sumFrom += wFrom;
                }
              }
            }

            std::set< std::pair<int, int> > &toList = it->second;

            //number of variables receiving mass and respecting coarser
            //partition structure 
            int nTo = 0;
            int nToPrev = 0;
            double sumTo = 0;

            for(std::set< std::pair<int,int> >::iterator toIt = toList.begin(); toIt!=toList.end(); ++toIt){
              const std::pair<int, int> &to = *toIt;
              double W = prevSol->W[from][to];
            

              if(W != 0){
                nToPrev++;
                int ti = to.first*2;
                int tj = to.second*2;
                for(int i=0; i<2; i++){
                  for(int j=0; j<2; j++){
                    double wTo = targetP(ti+i, tj+j);
                    if(wTo > 0){
                      nTo++;
                      sumTo += wTo;
                    }
                  }
                }
              }
            }

            /*
            double diff = 0;
            if(sumTo < sumFrom){
              double diff = sumFrom - sumTo;
              for(int i=0; i<2; i++){
                for(int j=0; j<2; j++){
                  double wFrom = sourceP(fi+i, fj+j); 
                  if(wFrom > diff){
                    sourceP(fi+i, fj+j) -= diff;
                    std::cout << "adjusted outgoing mass by: " << diff <<
                      std::endl;
                    i=10;
                    break;
                  }
                }
              }
              sumFrom = sumTo;
            }
            else if(sumFrom > sumTo){
              std::cout << "moo" << std::endl;
            }
*/




            //setup initializaing lp program
            glp_prob *lpInit = glp_create_prob();
            glp_set_obj_dir(lpInit, GLP_MIN);

            glp_add_rows(lpInit, nFrom + nTo + nToPrev);
            glp_add_cols(lpInit, nFrom * nTo );



            //set up coefficents (i.e. variable bounds and coefficients)
            int col= 1;
            for(int i=0; i<2; i++){
              for(int j=0; j<2; j++){
                double wFrom   = sourceP(fi+i, fj+j);

                
                if(wFrom > 0){
                  
                  std::pair<int, int> f(fi+i, fj+j);
                  
                  for(std::set< std::pair<int,int> >::iterator toIt = toList.begin(); toIt!=toList.end(); ++toIt){
                    const std::pair<int, int> &to = *toIt;
                    double W = prevSol->W[from][to];


                    if( W > 0 ){

                      int ti = to.first*2;             
                      int tj = to.second*2;

                      for(int ii=0; ii<2; ii++){
                        for(int jj=0; jj<2; jj++){
                          double wTo = targetP(ti+ii, tj+jj);

                          if(wTo > 0){
                            std::pair<int, int> t(ti+ii, tj+jj);


                            //set column coefficients and bounds
                            glp_set_col_bnds(lpInit, col, GLP_LO, 0, 0);
                            glp_set_obj_coef(lpInit, col, distance(source, target, f, t, p) );
                            ++col;

                          }

                        }
                      }
                    }

                  }
                }
              }
            }

            //Setup rows (i.e. constraints)

            //source constraints
            int row = 1;
            std::vector<int> indTo(nTo+1, 0);
            std::vector<double> valTo(nTo+1, 1);
            for(int i=0; i<2; i++){
              for(int j=0; j<2; j++){
                double wFrom = sourceP(fi+i, fj+j);
                if(wFrom > 0){
                  //wFrom *= fromFactor;

                  int off = (row-1) * nTo; 
                  for(int n=1; n<=nTo; n++){
                    indTo[n] = off + n;
                  }

                  glp_set_row_bnds(lpInit, row, GLP_FX, wFrom, wFrom);
                  glp_set_mat_row(lpInit, row, nTo, indTo.data(), valTo.data());
                  ++row;
                }
              }
            }



            //target constraints
            int targetInd = 1;
            std::vector<int> indFrom(nFrom+1, 0);
            std::vector<double> valFrom(nFrom+1, 1);
            for(std::set< std::pair<int,int> >::iterator toIt = toList.begin(); toIt!=toList.end(); ++toIt){
              const std::pair<int, int> &to = *toIt;
              double W = prevSol->W[from][to];

              int nToP = 0;
              if( W > 0 ){

                int ti = to.first*2;             
                int tj = to.second*2;
                
                for(int ii=0; ii<2; ii++){
                  for(int jj=0; jj<2; jj++){
                    double wTo = targetP(ti+ii, tj+jj);
                    if(wTo > 0){
                      for(int n=0; n<nFrom; n++){
                        indFrom[n+1] = targetInd + n * nTo;
                      }
                      glp_set_row_bnds(lpInit, row, GLP_UP, wTo, wTo);
                      glp_set_mat_row(lpInit, row, nFrom, indFrom.data(),
                          valFrom.data());
                      row++;
                      targetInd++;
                      nToP++;
                    }
                  }
                }

                //std::cout << "W-sumw: " << W - sumw << std::endl;
              

                if(nToP != 0){ 
                  //Add constraint from above partition            
                  std::vector<int> indP(nFrom*nToP +1, 0);
                  std::vector<double> valP(nFrom*nToP+1,1);
                  for(int n=0; n < nFrom; n++){
                    int off = n * nTo;
                    for(int i=0; i < nToP; i++){
                      indP[n*nToP + i+1] = off + targetInd - 1 - i;
                    }
                  }
                  glp_set_row_bnds(lpInit, row, GLP_FX, W, W);
                  glp_set_mat_row(lpInit, row, nToP*nFrom, indP.data(), valP.data());
                  row++;
                }
                else{
                  std::cout << "Partitioned zeroed out" << std::endl;
                  glp_set_row_bnds(lpInit, row, GLP_LO, 0, 0);
                  row++;
                }
              }
            }

            glp_smcp parm;
            glp_init_smcp(&parm);
            parm.msg_lev = GLP_MSG_ERR;

            int ret = glp_simplex(lpInit, &parm);
         
            //read basis rows from solution

            col= 1;
            PMap targetUpdates;
            for(int i=0; i<2; i++){
              for(int j=0; j<2; j++){
                double wFrom   = sourceP(fi+i, fj+j);
                
                if(wFrom > 0){      
                  std::pair<int, int> f(fi+i, fj+j);
                  
                  for(std::set< std::pair<int,int> >::iterator toIt = toList.begin(); toIt!=toList.end(); ++toIt){
                    
                    const std::pair<int, int> &to = *toIt;
                    double W = prevSol->W[from][to];


                    if( W > 0 ){

                      int ti = to.first*2;             
                      int tj = to.second*2;


                      for(int ii=0; ii<2; ii++){
                        for(int jj=0; jj<2; jj++){
                          double wTo = targetP(ti+ii, tj+jj);

                          if(wTo > 0){
                            std::pair<int, int> t(ti+ii, tj+jj);

                            
                            int status = glp_get_col_stat(lpInit, col);


                            if(status == GLP_BS){
                              if(nBaseCol < nConstraints - 1){
                                Path p(f, t);
                                int colIndex = pmap[p];
                                glp_set_col_stat(lp, colIndex, status);
                              }
                              nBaseCol++;
                            }
                            
                            double w = glp_get_col_prim(lpInit, col);

                            if(w > 0 ){
                              targetUpdates[t] = targetUpdates[t] + w;
                            }

                            col++;

                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            for(PMap::iterator it = targetUpdates.begin(); it !=
                targetUpdates.end(); ++it){
              const std::pair<int, int> &t = it->first;     
              double w = targetP(t.first, t.second) - it->second;
              if(w < parm.tol_bnd){
                w = 0;
              }
              targetP(t.first, t.second) = w;
            }


            glp_delete_prob(lpInit);

          }
         


          //Now also set all constraints as active
          for(int i=0; i< nConstraints; i++){
            glp_set_row_stat(lp, i+1, GLP_NU);
          }
/*
          int nBaseRow = 0;
          for(int i=0; i < targetP.M(); i++){
            for(int j=0; j < targetP.N(); j++){
              if(targetP(i, j) < 0 ){
                nBaseRow++;
                std::cout << targetP(i, j) << std::endl;
              }
            }
          }
  */


          //one constraint is linearly dependant on the rest
          glp_set_row_stat(lp, 1, GLP_BS);
          //for(int i=0; i < nConstraints - 1 - nBaseCol; i++){
          //  glp_set_row_stat(lp, nConstraints-i, GLP_BS);
          //}


          //failed to setup intial solution
          if(nConstraints - nBaseCol != 1){
            glp_std_basis(lp);
            //glp_adv_basis(lp, 0);
          }
          
          std::cout << "#basic columns: " << nBaseCol << std::endl; 
          std::cout << "#constraints: " << nConstraints << std::endl; 
          std::cout << "#variables: " << nCombinations << std::endl; 

          targetP.deallocate();
          sourceP.deallocate();

        }
        //
        //end initialization from previous scale
        //







        //
        //Solve current scale lp
        //

        std::cout << "Solving LP of size: " << source->p.M() * source->p.N() << " x " <<
          target->p.M() * target->p.N() << std::endl; 
        std::cout << "Reduced to " << nCombinations << " variables" << std::endl; 
         
        // set lp parameters
        glp_smcp parm;
        glp_init_smcp(&parm);
        //parm.meth = GLP_DUAL;

        //glp_adv_basis(lp, 0);

        //solve lp with simplex algorithm
        int ret = glp_simplex(lp, &parm);
        //int ret = glp_interior(lp, NULL);

        //sol.optimizationStatus = ret;
        sol->cost = glp_get_obj_val(lp);

        //resolve while increasing number of variables by including neighborhood
        //transport paths
        double prevCost = 0;
        while(prevCost != sol->cost){
          int added = addNeighborhodColums(lp, imap, fromMap, toMap, source,
              target, sol, p);
          prevCost = sol->cost;
          if(added != 0){
            nCombinations += added;
            std::cout << "#variables: " << nCombinations << std::endl; 
            int ret = glp_simplex(lp, &parm);
            sol->cost = glp_get_obj_val(lp);
          }
            
        }



        //store solution
        index = 1;
        for(VMap::iterator it = sol->source.begin(); it != sol->source.end(); ++it){
          const std::pair<int, int> &from = it->first;
          sol->psi[from] = glp_get_row_dual(lp, index);
          ++index;
        }


        for(VMap::iterator it = sol->target.begin(); it != sol->target.end(); ++it){
          const std::pair<int, int> &to = it->first;
          sol->phi[to] = glp_get_row_dual(lp, index);
          index++;
        }

        
        int ncols = glp_get_num_cols(lp);
        for(int i=1; i<= ncols; i++){
          const Path &path = imap[i];
          sol->W[path.first][path.second] = glp_get_col_prim(lp, i);
          sol->delta[path.first][path.second] = glp_get_col_dual(lp, i);
          //std::cout << sol->phi[to] << ", " << sol->psi[from] << ", " << sol->delta[from][to] << " : " << sol->D[from][to] << std::endl;
        } 
       





        //Create list of neighboring solutions to optimal one
        if(neighborhood == POLY){

          typedef Eigen::Triplet<double> Triplet;

          std::vector<Triplet> basic;
          basic.reserve( ( nConstraints - 1) * nConstraints );
          std::vector<Triplet> nonbasic;
          nonbasic.reserve( ( nCombinations - nConstraints +1 ) * nConstraints );

          int *ind = new int[nConstraints];
          double *val =  new double[nConstraints];
          int basicCol = 0;
          int nonbasicCol = 0;

          DenseVector<int> basicIndices(nConstraints);
          DenseVector<int> nonBasicIndices(nCombinations - nConstraints + 1);
          for(int i=0; i<nCombinations; i++){
            int status = glp_get_col_stat(lp, i+1);
            int l = glp_get_mat_col(lp, i+1, ind, val);
            double v = 1.0/sqrtf(l);
            if(status == GLP_BS){
              for(int j=1; j<=l; j++){
                Triplet t(ind[j]-1, basicCol, v);
                basic.push_back(t);
              }
              basicIndices(basicCol) = basicCol+nonbasicCol;
              basicCol++;
            }
            else{
              for(int j=1; j<=l; j++){
                //build transpose of non-basic
                Triplet t(nonbasicCol, ind[j] - 1, v);
                nonbasic.push_back(t);
              }
              nonBasicIndices(nonbasicCol) = basicCol+nonbasicCol;
              nonbasicCol++;
            }
          }
          delete ind;
          delete val;

          Eigen::SparseMatrix<double> Basic(nConstraints, nConstraints-1);
          Basic.setFromTriplets(basic.begin(), basic.end() );

          Eigen::SparseMatrix<double> NonBasic(nCombinations-nConstraints+1,
              nConstraints);
          NonBasic.setFromTriplets(nonbasic.begin(), nonbasic.end() );

          Eigen::SparseMatrix<double> P = (NonBasic * Basic).pruned();

          std::cout << Basic.rows() << "x" << Basic.cols() << std::endl;
          std::cout << NonBasic.rows() << "x" << NonBasic.cols() << std::endl;
          std::cout << P.rows() << "x" << P.cols() << std::endl;

          DenseVector<double> length( P.rows() );
          Linalg<double>::Zero(length);


          for (int k=0; k < P.outerSize(); ++k){
            for (Eigen::SparseMatrix<double>::InnerIterator it(P,k); it; ++it){
              double v = it.value();
              int r = it.row(); // row index
              int c = it.col(); // col index (here it is equal to k)
              length(r) += v*v;
            }
          }

          int nn = 0;
          for(int i=0; i < length.N(); i++){
            if(length(i) > 0.99){
              Path &path = imap[nonBasicIndices(i)+1];
              sol->neighbors[path.first].insert(path.second);
              ++nn;
            }
          }

          basicIndices.deallocate();
          nonBasicIndices.deallocate();


          std::cout << std::endl << "nn: " << nn << std::endl;
          std::cout << std::endl << "all-nn: " << length.N() - nn << std::endl;
          std::cout << "basicCol: " << basicCol << std::endl;
          std::cout << "nConstraints: " << nConstraints << std::endl;


        }




        //Clean up and store source and target distribution at this scale in
        //solution object
        glp_delete_prob(lp);

        sol->sourceHist = source;
        sol->targetHist = target;


        return sol;
    };




    static int addNeighborhodColums(glp_prob *lp, IMap &imap, RMap &fromMap,
        RMap &toMap, ScaleHistogram<TPrecision> *source,
        ScaleHistogram<TPrecision> *target, EMD<TPrecision> *sol, int p){

        int ncols = glp_get_num_cols(lp);
        int nAdded = 0;
        for(int i=1; i <= ncols; i++){
          int status = glp_get_col_stat(lp, i);

          if(status == GLP_BS){
            const Path &path = imap[i]; 
            const std::pair<int, int> &f1= path.first;
            const std::pair<int, int> &t1 = path.second;
            
            std::pair<int, int> f2;
            std::pair<int, int> t2;
            for(int i= -1; i <= 1; i++){
              f2.first = f1.first+i;
              for(int j= -1; j <= 1; j++){
                f2.second = f1.second + j;
                if( isInBounds(f2, source)  && source->p(f2.first, f2.second) != 0){

                  for(int ii= -1; ii <= 1; ii++){
                    t2.first = t1.first+ii;
                    for(int jj= -1; jj <= 1; jj++){
                      t2.second = t1.second + jj;

                      if( isInBounds(t2, target)   && target->p(t2.first, t2.second) != 0 && !entryExists(f2, t2, sol->D, -1) ){

                        //add ccolumn to lp
                        int col = glp_add_cols(lp, 1);
                       
                        //set coefficient 
                        double d = distance(source, target, f2, t2, p);
                        glp_set_obj_coef(lp, col, d);
                        glp_set_col_bnds(lp, col, GLP_DB, 0, 1);

                        //set constraints
                        int ind[3] = {0,0,0};
                        double val[3] = {1,1,1};
                        ind[1] = fromMap[f2];
                        ind[2] = toMap[t2];
                        glp_set_mat_col(lp, col, 2, ind, val);
                        
                        //update maps
                        fromMap[f2] = ind[1];
                        toMap[t2] = ind[2];
                        
                        imap[col] = Path(f2, t2);
                      
                        sol->source[f2].insert(t2);
                        sol->target[t2].insert(f2);
                        
                        sol->D[f2][t2] = d;

                        
                        nAdded++;
                      }

                    }
                  }

                }

              } 
            } 
          
          }
        }

        return nAdded;

    };






    static int BallNeighborhood(ScaleHistogram<TPrecision> *source,
        ScaleHistogram<TPrecision> *target, EMD<TPrecision> *sol, EMD<TPrecision>
        *prevSol, int r, int p){
      //Number of variables in LP. I.e. number of combinations of indices in A
      //to indices in B with nonzero weight at previous scale
      int nCombinations = 0;


      for(DMap::iterator it = prevSol->W.begin(); it != prevSol->W.end(); ++it){
        std::pair<int, int> f1  = it->first;
        f1.first *= 2;
        f1.second *= 2;

        std::map< std::pair<int, int>, double > &toList = it->second;
        for(std::map< std::pair<int, int>, double >::iterator toIt = toList.begin(); toIt != toList.end(); ++toIt){

          //No mass moved at previous solution 
          if(toIt->second == 0){
            continue;
          }

          std::pair<int, int> t1 = toIt->first;
          t1.first *= 2;
          t1.second *= 2;

          std::pair<int, int> f2;
          std::pair<int, int> t2;

          for(int i= -r; i <= r+1; i++){
            f2.first = f1.first+i;
            for(int j= -r; j <= r+1; j++){
              f2.second = f1.second + j;
              if( isInBounds(f2, source)  && source->p(f2.first, f2.second) != 0){

                for(int ii= -r; ii <= r+1; ii++){
                  t2.first = t1.first+ii;
                  for(int jj= -r; jj <= r+1; jj++){
                    t2.second = t1.second + jj;

                    if( isInBounds(t2, target)   && target->p(t2.first, t2.second) != 0 && !entryExists(f2, t2, sol->D, -1) ){

                      sol->D[f2][t2] = distance(source, target, f2, t2, p);

                      sol->source[f2].insert(t2);
                      sol->target[t2].insert(f2);
                      nCombinations++;
                    }

                  }
                }

              }

            } 
          } 

        }
      }

      return nCombinations;
    };



    static int PolyNeighborhood(ScaleHistogram<TPrecision> *source,
        ScaleHistogram<TPrecision> *target, EMD<TPrecision> *sol, EMD<TPrecision>
        *prevSol, int r, int p){
      //Number of variables in LP. I.e. number of combinations of indices in A
      //to indices in B with nonzero weight at previous scale
      int nCombinations = 0;

       
        for(DMap::iterator it = prevSol->W.begin(); it != prevSol->W.end(); ++it){
          std::pair<int, int> f1  = it->first;
          f1.first *= 2;
          f1.second *= 2;
          
          std::map< std::pair<int, int>, double > &toList = it->second;
          for(std::map< std::pair<int, int>, double >::iterator toIt = toList.begin(); toIt != toList.end(); ++toIt){
           
            //No mass moved at previous solution 
            if(toIt->second == 0){
              continue;
            }

            std::pair<int, int> t1 = toIt->first;
            t1.first *= 2;
            t1.second *= 2;

            std::pair<int, int> f2;
            std::pair<int, int> t2;
 
            for(int i= -r; i <= r+1; i++){
              f2.first = f1.first+i;
              for(int j= -r; j <= r+1; j++){
                f2.second = f1.second + j;
                if( isInBounds(f2, source)  && source->p(f2.first, f2.second) != 0){
                  
                  for(int ii=-r; ii <= r+1; ii++){
                    t2.first = t1.first+ii;
                    for(int jj=-r; jj <= r+1; jj++){
                      t2.second = t1.second + jj;
                    
                       if( isInBounds(t2, target)   && target->p(t2.first, t2.second) != 0 && !entryExists(f2, t2, sol->D, -1) ){

                         sol->D[f2][t2] = distance(source, target, f2, t2, p);

                         sol->source[f2].insert(t2);
                         sol->target[t2].insert(f2);
                         nCombinations++;
                      }

                    }
                  }

                }

              } 
            } 

          }
        }
       
        for(VMap::iterator it = prevSol->neighbors.begin(); it !=
            prevSol->neighbors.end(); ++it){
          std::pair<int, int> f1  = it->first;
          f1.first *= 2;
          f1.second *= 2;
          
          std::set< std::pair<int, int> > &toList = it->second;
          for(std::set< std::pair<int, int>, double >::iterator toIt = toList.begin(); toIt != toList.end(); ++toIt){
           

            std::pair<int, int> t1 = *toIt;
            t1.first *= 2;
            t1.second *= 2;

            std::pair<int, int> f2;
            std::pair<int, int> t2;
 
            for(int i=0; i <= 1; i++){
              f2.first = f1.first+i;
              for(int j= 0; j <= 1; j++){
                f2.second = f1.second + j;
                if( isInBounds(f2, source)  && source->p(f2.first, f2.second) != 0){
                  
                  for(int ii= 0; ii <= 1; ii++){
                    t2.first = t1.first+ii;
                    for(int jj=0; jj <= 1; jj++){
                      t2.second = t1.second + jj;
                    
                       if( isInBounds(t2, target)   && target->p(t2.first, t2.second) != 0 && !entryExists(f2, t2, sol->D, -1) ){

                         sol->D[f2][t2] = distance(source, target, f2, t2, p);

                         sol->source[f2].insert(t2);
                         sol->target[t2].insert(f2);
                         nCombinations++;
                      }

                    }
                  }

                }

              } 
            } 

          }
        }

        return nCombinations;
    };






    static int DualNeighborhood(ScaleHistogram<TPrecision> *source,
        ScaleHistogram<TPrecision> *target, EMD<TPrecision> *sol, EMD<TPrecision>
        *prevSol, int rFactor, int p){
        
        int nCombinations = 0;

        double r = (source->radius + target->radius) * 2 * rFactor;

        //sources to targets 
        for(VMap::iterator sIt = prevSol->source.begin(); sIt !=
            prevSol->source.end(); sIt++){
          
          const std::pair<int, int> &from = sIt->first;
         
          std::set< std::pair<int, int> > &targets = sIt->second;
          
          
          for(int index = 0; index < targets.size() - 1; index++){
             std::set< std::pair<int, int> >::iterator tIt  = targets.begin();
             for(int i=0; i<index; i++){
               ++tIt;
             }
             const std::pair<int, int> &to = *tIt;
             //TPrecision dist = prevSol->D[from][to];
             TPrecision delta = prevSol->delta[from][to];
             TPrecision pTo = prevSol->phi[to];
             
             if(delta != 0){ continue; }
             for(;tIt != targets.end(); ++tIt){
               const std::pair<int, int> &to2 = *tIt;
               
               //TPrecision dist2 = prevSol->D[from][to2];
               //TPrecision delta2 = prevSol->delta[from][to2];
               TPrecision pTo2 = prevSol->phi[to2];
               //if(dist2 - delta2 + r > dist + delta ){
                 //Add all pairs of the childrens partitions of from and to2
                 nCombinations += addSubPartitionsSource(source, target, from,
                     to, to2, sol, pTo, pTo2, p);
               //}
             } 
          } 
        }

        //targets to sources
        for(VMap::iterator tIt = prevSol->target.begin(); tIt !=
            prevSol->target.end(); tIt++){
          
          const std::pair<int, int> &to = tIt->first;
          
          std::set< std::pair<int, int> > &sources = tIt->second;
          
          
          for(int index= 0; index < sources.size() - 1; index++){
             std::set< std::pair<int, int> >::iterator sIt  = sources.begin();
             for(int i=0; i<index; i++){
               ++sIt;
             }
             const std::pair<int, int> &from = *sIt;
             
             //TPrecision dist = prevSol->D[from][to];
             TPrecision delta = prevSol->delta[from][to];
             TPrecision pFrom = prevSol->psi[from];
             
             if(delta != 0 ){ continue; }
             for(;sIt != sources.end(); ++sIt){
               const std::pair<int, int> &from2 = *sIt;
               
               //TPrecision dist2 = prevSol->D[from2][to];
               //TPrecision delta2 = prevSol->delta[from2][to];
               TPrecision pFrom2 = prevSol->psi[from2];
               //if(dist2 - delta2 + r > dist + delta ){
                 //Add all pairs of the childrens partitions of from2 and to
                 nCombinations += addSubPartitionsTarget(source, target, from, from2, to,
                     sol, -pFrom, -pFrom2, p);
               //}
             } 
          } 
        }

        return nCombinations;
    };





    static int addSubPartitionsSource(ScaleHistogram<TPrecision> *source,
        ScaleHistogram<TPrecision> *target, const std::pair<int, int> &from,
        const std::pair<int, int> &to,  const std::pair<int, int> to2,
        EMD<TPrecision> *sol, double pTo1, double pTo2, int p){

      int fi = from.first*2;
      int fj = from.second*2;
      
      int ti = to.first*2;
      int tj = to.second*2;
      
      int t2i = to2.first*2;
      int t2j = to2.second*2;
      
      int nAdded = 0;

      for(int i=0;i<2;i++){
        for(int j=0; j<2;j++){
          std::pair<int, int> f(fi+i, fj+j);
          if( source->p(f.first, f.second) == 0 ) continue;

          for(int ii=0;ii<2;ii++){
            for(int jj=0; jj<2;jj++){
              std::pair<int, int> t(ti+ii, tj+jj);
              if( target->p(t.first, t.second) == 0 ) continue;
                
              TPrecision dist = distance(source, target, f, t, p);
              
              for(int iii=0; iii<2; iii++){
                for(int jjj=0; jjj<2; jjj++){
                  std::pair<int, int> t2(t2i+iii, t2j+jjj);
                  if( target->p(t2.first, t2.second) == 0 ) continue;
          
                  if( !entryExists(f, t2, sol->D, -1) ){

                    TPrecision dist2 = distance(source, target, f, t2, p);

                    if(dist2 + pTo2 <= dist + pTo1 ){
                      sol->D[f][t2] = dist;
                      sol->source[f].insert(t2);
                      sol->target[t2].insert(f);
                      nAdded++;
                    }
                  }
                }

              }
            }
          }
        }
      }
      return nAdded;
    };






    static int addSubPartitionsTarget(ScaleHistogram<TPrecision> *source,
        ScaleHistogram<TPrecision> *target, const std::pair<int, int> &from,
        const std::pair<int, int> &from2,  const std::pair<int, int> to,
        EMD<TPrecision> *sol, double pFrom, double pFrom2, int p){

      int fi = from.first*2;
      int fj = from.second*2;

      int ti = to.first*2;
      int tj = to.second*2;
      
      int f2i = from2.first*2;
      int f2j = from2.second*2;

      int nAdded = 0;
      for(int i=0;i<2;i++){
        for(int j=0; j<2;j++){
          std::pair<int, int> f(fi+i, fj+j);
          if( source->p(f.first, f.second) == 0 ) continue;

          for(int ii=0;ii<2;ii++){
            for(int jj=0; jj<2;jj++){
              std::pair<int, int> t(ti+ii, tj+jj);
              if( target->p(t.first, t.second) == 0 ) continue;
                
              TPrecision dist = distance(source, target, f, t, p);
              
              for(int iii=0; iii<2; iii++){
                for(int jjj=0; jjj<2; jjj++){
                  std::pair<int, int> f2(f2i+iii, f2j+jjj);
                  if( source->p(f2.first, f2.second) == 0 ) continue;
          
                  if( !entryExists(f2, t, sol->D, -1) ){

                    TPrecision dist2 = distance(source, target, f2, t, p);

                    if(dist2 - pFrom2 <= dist - pFrom ){
                      sol->D[f2][t] = dist;
                      sol->source[f2].insert(t);
                      sol->target[t].insert(f2);
                      nAdded++;
                    }
                  }
                }

              }
            }
          }
        }
      }
      return nAdded;
    };







    static std::list< ScaleHistogram<TPrecision> * > buildPyramid(DenseMatrix<TPrecision> X){
       ScaleHistogram<TPrecision> *H = new ScaleHistogram<TPrecision>(X.M(), X.N());
       
       for(int i=0; i<H->p.M(); i++){
         for(int j=0; j<H->p.N(); j++){
            H->p(i, j) = X(i, j);
            H->mI(i, j) = i;
            H->mJ(i, j) = j;
          }
        }
        H->radius = 0.5;
  
        std::list< ScaleHistogram<TPrecision> * > pyramid;
        pyramid.push_front(H); 
        while( std::min(H->p.M(), H->p.N()) > 2 ){
          
          ScaleHistogram<TPrecision> *H2 = new
            ScaleHistogram<TPrecision>(H->p.M()/2, H->p.N()/2);
          H2->radius = H->radius*2;
          
          
          for(int i=0; i<H2->p.M(); i++){
            for(int j=0; j<H2->p.N(); j++){
              TPrecision p1 = H->p(2*i, 2*j);
              TPrecision p2 = H->p(2*i+1, 2*j);
              TPrecision p3 = H->p(2*i, 2*j+1);
              TPrecision p4 = H->p(2*i+1, 2*j+1);
              TPrecision ps = p1+p2+p3+p4;
              H2->p(i, j) = ps;
              if(ps!=0){
                H2->mI(i, j) = (p1*H->mI(2*i, 2*j) +  p2*H->mI(2*i+1, 2*j) +  p3*H->mI(2*i, 2*j+1) +  p4*H->mI(2*i+1, 2*j+1) ) / ps;
                H2->mJ(i, j) = (p1*H->mJ(2*i, 2*j) +  p2*H->mJ(2*i+1, 2*j) +  p3*H->mJ(2*i, 2*j+1) +  p4*H->mJ(2*i+1, 2*j+1) ) / ps;
              }
              else{
                H2->mI(i, j) = (H->mI(2*i, 2*j) +  H->mI(2*i+1, 2*j) +  H->mI(2*i, 2*j+1) +  H->mI(2*i+1, 2*j+1) ) / 4.0;
                H2->mJ(i, j) = (H->mJ(2*i, 2*j) +  H->mJ(2*i+1, 2*j) +  H->mJ(2*i, 2*j+1) +  H->mJ(2*i+1, 2*j+1) ) / 4.0;
              }
            }
          }


          //TODO: add left overs
          pyramid.push_front(H2);
          H = H2;
        }
 
       return pyramid; 
    };





    

    static bool entryExists(std::pair<int, int> &from, std::pair<int, int> &to, DMap &pW, double tol){
      DMap::iterator it = pW.find(from);
      if(it == pW.end() ){
        return false;
      }
      std::map< std::pair<int, int>, double> &m2 = it->second;
      std::map< std::pair<int, int>, double>::iterator it2 = m2.find(to);
      if(it2 == m2.end() ){
        return false;
      }     
      double w = it2->second;
      return w > tol;

    };





    static bool isInBounds(std::pair<int, int> &index,
        ScaleHistogram<TPrecision> *X){
      if(index.first < 0){
        return false;
      }
      if(index.first >= X->p.M() ){
        return false;
      }
      if(index.second < 0){
        return false;
      }
      if(index.second >= X->p.N() ){
        return false;
      } 
      return true;
    };





    static TPrecision distance(ScaleHistogram<TPrecision> *source,
        ScaleHistogram<TPrecision> *target, const std::pair<int, int> &from,
        const std::pair<int, int> &to, int p){
      double x = source->mI(from.first, from.second) - target->mI(to.first, to.second);
      double y = source->mJ(from.first, from.second) - target->mJ(to.first, to.second);
      double dist = sqrt(x*x + y*y);

      return pow(dist, p);
    };


};



         
          /*
          int nBaseCol = 0;
          for(VMap::iterator it = prevSol->source.begin(); it !=
              prevSol->source.end(); ++it){
            const std::pair<int, int> &from = it->first;

            std::set< std::pair<int, int> > &toList = it->second;
            for(std::set< std::pair<int,int> >::iterator toIt = toList.begin(); toIt!=toList.end(); ++toIt){
              const std::pair<int, int> &to = *toIt;
              double W = prevSol->W[from][to];


              if( W != 0 ){
                //solve lp to figure out intializiation with this partition
                int fi = from.first*2; 
                int fj = from.second*2; 

                int ti = to.first*2; 
                int tj = to.second*2;
                
                for(int i=0; i<2; i++){
                  for(int j=0; j<2; j++){
                    if(target->p(ti+i, tj+j) != 0){
                      ti+=i;
                      tj+=j;
                      break;
                    }
                  }
                  if(target->p(ti, tj) != 0) break;
                }
                std::pair<int, int> t(ti, tj);

                for(int i=0; i<2; i++){
                  for(int j=0; j<2; j++){
                    if(source->p(fi+i, fj+j) != 0){
                      std::pair<int, int> f(fi+i, fj+j);
                      int cIndex = varmap[f][t];
                      glp_set_col_stat(lp, cIndex, GLP_BS);
                      nBaseCol++;
                    }
                  }
                }
              }
            }
          }
         
          for(VMap::iterator it = sol->target.begin(); it !=
             sol->target.end(); ++it){
            const std::pair<int, int> &to = it->first;
            std::set< std::pair<int, int> > &fromList = it->second;
            const std::pair<int, int> &from = *fromList.rbegin();
            int col = varmap[from][to];
            glp_set_col_stat(lp, col, GLP_BS);
            nBaseCol++;
          }
          
          for(int i=0; i< sol->source.size(); i++){
            glp_set_row_stat(lp, i+1, GLP_BS);
          }
          for(int i=0; i< sol->target.size(); i++){
            glp_set_row_stat(lp, sol->source.size() + i + 1, GLP_NU);
          }
          std::cout << "#basic columns: " << nBaseCol << std::endl; 
          std::cout << "#constraints: " << sol->source.size() << std::endl;
          std::cout << "#constraints: " << sol->target.size() << std::endl;
           */
           

   /*       
          DenseMatrix<TPrecision> targetP = Linalg<TPrecision>::Copy( target->p );
          DenseMatrix<TPrecision> sourceP = Linalg<TPrecision>::Copy( source->p );
          
          int nBaseCol = 0;
          for(VMap::iterator it = sol->source.begin(); it != sol->source.end(); ++it){
          
            const std::pair<int, int> &from = it->first;

            std::set< std::pair<int, int> > &toList = it->second;
            double sumTo = 0;
            for(std::set< std::pair<int,int> >::iterator toIt = toList.begin();
                toIt!=toList.end(); ++toIt, ++row){
              const std::pair<int, int> &to = *toIt;
              sumTo += targetP(to.first, to.second);
            }

            glp_prob *lpInit = glp_create_prob();
            glp_set_obj_dir(lpInit, GLP_MIN);

            glp_add_rows(lpInit, toList.size() + 1 );
            glp_add_cols(lpInit, toList.size()  );

            double wFrom = source->p(from.first, from.second);
            
            bool massExcess = wFrom > sumTo;

            if(massExcess){
              glp_set_row_bnds(lpInit, 1, GLP_DB, 0, wFrom);
            }
            else{
              glp_set_row_bnds(lpInit, 1, GLP_FX, wFrom, wFrom);
            }
            int toInd[toList.size()+1];
            double toVal[toList.size()+1];
            for(int i=0; i< toList.size(); i++){
              toInd[i+1] = i+1;
              toVal[i+1] = 1;
            }
            glp_set_mat_row(lpInit, 1, toList.size(), toInd, toVal);
            
            int row=2;
            double val[] = {1,1};
            for(std::set< std::pair<int,int> >::iterator toIt = toList.begin();
                toIt!=toList.end(); ++toIt, ++row){
              
              const std::pair<int, int> &to = *toIt;
            
              double wTo = targetP(to.first, to.second);
              if(wTo != 0 && !massExcess){
                glp_set_row_bnds(lpInit, row, GLP_DB, 0, wTo);
              }
              else{
                glp_set_row_bnds(lpInit, row, GLP_FX, wTo, wTo);
              }

              int fromInd[] = {0, row-1};
              glp_set_mat_row(lpInit, row, 1, fromInd, val);
                        
              glp_set_col_bnds(lpInit, row-1, GLP_DB, 0, 1);
              glp_set_obj_coef(lpInit, row-1, distance(source, target, from, to, p) );
            
            }
                      
            glp_smcp parm;
            glp_init_smcp(&parm);
            parm.msg_lev = GLP_MSG_ERR;

            int ret = glp_simplex(lpInit, &parm);
            
            int col=1;
            for(std::set< std::pair<int,int> >::iterator toIt = toList.begin();
                toIt!=toList.end(); ++toIt, ++row){

              const std::pair<int, int> &to = *toIt;
              double W = glp_get_col_prim(lpInit, col);

              if(W != 0){ 
                targetP(to.first, to.second) -= W;
                sourceP(from.first, from.second) -= W;

                int colstat = glp_get_col_stat(lpInit, col);
                //std::cout << "W: " << W << std::endl;
                //std::cout << "colstat: " << colstat << std::endl;

                int cIndex = varmap[from][to];
                glp_set_col_stat(lp, cIndex, colstat);
                if(colstat == GLP_BS){
                  nBaseCol++;
                }
              }
              ++col;
            }
                
            
            glp_delete_prob(lpInit);
          }

          for(int i=0; i< sol->source.size() + sol->target.size(); i++){
            glp_set_row_stat(lp, i+1, GLP_NU);
          }

          //check for violated constraints
          std::list< std::pair<double, int> > violations; 
          int rIndex = 1;
          for(VMap::iterator it = sol->source.begin(); it != sol->source.end(); ++it){
            const std::pair<int, int> &from = it->first;
            double wFrom = sourceP(from.first, from.second);
            if(wFrom != 0){ 
              std::pair<double, int> v(wFrom, rIndex);
              violations.push_back(v);
            }
            ++rIndex;
          }
          
          for(VMap::iterator it = sol->target.begin(); it != sol->target.end(); ++it){
            const std::pair<int, int> &to = it->first;
            double wTo = targetP(to.first, to.second);
            if(wTo != 0){ 
              std::pair<double, int> v(wTo, rIndex);
              violations.push_back(v);
            }
            ++rIndex;
          }
          
          violations.sort();
          
          int nConstraintBasis = 0;
          for(int i=0; i<nConstraints - nBaseCol; i++){
            const std::pair<double, int> &v = violations.back();
            violations.pop_back();

            glp_set_row_stat(lp, v.second, GLP_BS);
            nConstraintBasis++;
          }
          
          if(nConstraintBasis + nBaseCol != nConstraints){
            glp_set_row_stat(lp, 1, GLP_BS);
          }

          targetP.deallocate();
          sourceP.deallocate();


          std::cout << "#basic columns: " << nBaseCol << std::endl; 
          std::cout << "#basic rows: " << nConstraintBasis << std::endl; 
          std::cout << "#constraints: " << nConstraints << std::endl; 
          std::cout << "#variables: " << nCombinations << std::endl; 
*/




          /*
          int nBasePrev = 0;
          int nBaseCol = 0;
          for(VMap::iterator it = prevSol->source.begin(); it !=
              prevSol->source.end(); ++it){
            const std::pair<int, int> &from = it->first;

            int fi = from.first*2;
            int fj = from.second*2;

            std::set< std::pair<int, int> > &toList = it->second;

            for(std::set< std::pair<int,int> >::iterator toIt = toList.begin(); toIt!=toList.end(); ++toIt){
              const std::pair<int, int> &to = *toIt;
              double W = prevSol->W[from][to];


              if( W != 0 ){
                nBasePrev++;

                int ti = to.first*2;
                int tj = to.second*2;
              
                double pFrom = 0;
                double pTo = 0;
                for(int i=0; i<2; i++){
                  for(int j=0; j<2; j++){
                    pFrom += sourceP(fi+i, fj+j);
                    pTo   += targetP(ti+i, tj+j);
                  }
                }
                //solve lp to figure out intializiation within this partition

                glp_prob *lpInit = glp_create_prob();
                glp_set_obj_dir(lpInit, GLP_MIN);

                glp_add_rows(lpInit, 8 );
                glp_add_cols(lpInit, 16 );

                //add constraints and row indices
                int row = 1;
                double val[] = {0,1,1,1,1};
                for(int i=0; i<2; i++){
                  for(int j=0; j<2; j++){
                    double wFrom = sourceP(fi+i, fj+j);
                    double wTo = targetP(ti+i, tj+j);
                    if(pTo > pFrom){
                      glp_set_row_bnds(lpInit, row, GLP_FX, wFrom, wFrom);
                      if(wTo == 0){
                        glp_set_row_bnds(lpInit, 4+row, GLP_FX, 0, wTo);
                      }
                      else{
                        glp_set_row_bnds(lpInit, 4+row, GLP_DB, 0, wTo);
                      }
                    }
                    else{
                      if(wFrom == 0){
                        glp_set_row_bnds(lpInit, row, GLP_FX, 0, wFrom);
                      }
                      else{
                        glp_set_row_bnds(lpInit, row, GLP_DB, 0, wFrom);
                      }
                      glp_set_row_bnds(lpInit, 4+row, GLP_FX, wTo, wTo);
                    }

                    //glp_set_row_bnds(lpInit, row, GLP_FX, wFrom, wFrom);
                    //glp_set_row_bnds(lpInit, 4+row, GLP_FX, wTo, wTo);
                    
                    int c = (row-1)*4;
                    int toInd[] = {0, c+1, c+2, c+3, c+4};
                    int fromInd[] = {0, row, row+4, row+8, row+12};
                    glp_set_mat_row(lpInit, row, 4, toInd, val);
                    glp_set_mat_row(lpInit, 4+row, 4, fromInd, val);
                    ++row;
                  }
                }

                //add variable coefficients
                int col = 1;
                for(int i=0; i<2; i++){
                  for(int j=0; j<2; j++){
                    std::pair<int, int> f(fi+i, fj+j);
                    for(int ii=0; ii<2; ii++){
                      for(int jj=0; jj<2; jj++){
                        std::pair<int, int> t(ti+ii, tj+jj);
                        glp_set_col_bnds(lpInit, col, GLP_DB, 0, 1);
                        glp_set_obj_coef(lpInit, col, distance(source, target, f, t, p) );
                        ++col;
                      }
                    }
                  }
                } 

                int ret = glp_simplex(lpInit, NULL);

                //use above solution to set up intial solution for finer scale
                col = 1;
                double pFromP = prevSol->sourceHist->p(from.first, from.second);
                double pToP = prevSol->targetHist->p(to.first, to.second);
                
                double factor = 1;//W / std::min(pToP, pFromP);
                //factor = std::min(1.0, factor);

                for(int i=0; i<2; i++){
                  for(int j=0; j<2; j++){
                    std::pair<int, int> f(fi+i, fj+j);
                    for(int ii=0; ii<2; ii++){
                      for(int jj=0; jj<2; jj++){
                        
                        std::pair<int, int> t(ti+ii, tj+jj);
                        double W = glp_get_col_prim(lpInit, col);

                        if(W != 0){ 
                          targetP(t.first, t.second) -= W*factor;
                          sourceP(f.first, f.second) -= W*factor;

                          int colstat = glp_get_col_stat(lpInit, col);
                          //std::cout << "W: " << W << std::endl;
                          //std::cout << "colstat: " << colstat << std::endl;

                          int cIndex = varmap[f][t];
                          glp_set_col_stat(lp, cIndex, colstat);
                          if(colstat == GLP_BS){
                            nBaseCol++;
                          }
                        }
                        ++col;
                       
                      }
                    }
                  }
                }


                glp_delete_prob(lpInit);

              }
            }
          }
          sourceP.deallocate();
          targetP.deallocate();

          std::cout << "#basic columns previous: " << nBasePrev << std::endl; 
          std::cout << "#basic columns: " << nBaseCol << std::endl; 
          std::cout << "#constraints: " << nConstraints << std::endl; 
          std::cout << "#variables: " << nCombinations << std::endl; 
          //Now also set all constraints as active
          for(int i=0; i< sol->source.size() + sol->target.size(); i++){
            glp_set_row_stat(lp, i+1, GLP_NU);
          }
          for(int i=0; i < nConstraints - nBaseCol; i++){
            glp_set_row_stat(lp, i+1, GLP_BS);
          }
          */
           







#endif
