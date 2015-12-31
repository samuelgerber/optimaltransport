#include "Precision.h"

#include <iostream>

#include "LinalgIO.h"
#include "IPCATree.h"
#include "DenseMatrix.h"
#include <tclap/CmdLine.h>
#include "Wasserstein.h"


DenseMatrix<Precision> createMatrix(std::list<int> ind, DenseMatrix<Precision>
    X){ 
  DenseMatrix<Precision> M(X.M(), ind.size());
  int index = 0;
  for(std::list<int>::iterator it = ind.begin(); it != ind.end(); ++it){
    int i = *it;
    Linalg<Precision>::SetColumn(M, index, X, i);
    ++index;
  }
  return M;
};




int main(int argc, char **argv){

  //Command line parsing
  TCLAP::CmdLine cmd("GWT Wasserstein distances", ' ', "1");

  TCLAP::ValueArg<std::string> tArg("t","tree", "Tree data file",  true, "",
      "tree data file");
  cmd.add(tArg);

  TCLAP::ValueArg<std::string> xArg("x","data1", "Data set 1",  true, "",
      "data matrix header file");
  cmd.add(xArg);

  TCLAP::ValueArg<std::string> outArg("o","out", "Data file",  true, "", "file to store results in");
  cmd.add(outArg);

  try{
    cmd.parse( argc, argv );
  } 
  catch (TCLAP::ArgException &e){ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    return -1;
  }


  std::string outFileName  = outArg.getValue();
  std::string treeFileName = tArg.getValue();

  std::string xFileName = xArg.getValue();

  DenseMatrix<Precision> X = LinalgIO<Precision>::readMatrix(xFileName);

  IPCATree<Precision> ipcaTree;

  std::ifstream treeFile(treeFileName.c_str());
  ipcaTree.unflatten(treeFile);
  treeFile.close();

  //Compare densities in each node and store wasserstein distance and indices
  //into patch matrix

  std::list< GMRANode * > nodes;
  nodes.push_back( ipcaTree.getRoot() );

  std::list< std::list<int> > indices;

  std::list< int > start;
  for(int i=0; i<X.N(); i++){
    start.push_back(i);
  }
  indices.push_back(start);

  std::list<Precision> wdists;

  DenseVector<int> leafID(X.N());

  DenseMatrix< Precision > sumD(2, X.N()); 
  Linalg<Precision>::Zero(sumD);

  std::list<int> nodeID;
  nodeID.push_back(0);

  int IDcount = 0;
  int ID = 0; 
  while( !nodes.empty() ){
    IPCANode<Precision> *node = (IPCANode<Precision> *) nodes.front();
    nodes.pop_front();

    std::list<int> ind = indices.front();
    indices.pop_front();

    int curID = nodeID.front();
    nodeID.pop_front();

    /*
       DenseVector<int> treeInd = Linalg<int>::ToVector(node->getPoints());
       std::stringstream ss1;
       ss1 << outFileName << << "-node-" << nodeID << "-tree.data";
       LinalgIO<int>::writeVector(ss1.str(), treeInd);

       DenseVector<int> queryInd = Linalg<int>::ToVector(ind);
       std::stringstream ss2;
       ss2 << outFileName << << "-node-" << queryInd << "-query.data";
       LinalgIO<int>::writeVector(ss1.str(), treeInd);

       std::stringstream ss3;
       ss3 << outFileName << "-lid.data";
       LinalgIO<int>::writeVector(ss2.str(), leafID);
     */


    if(ind.size() > 0){
      DenseMatrix<Precision> A = createMatrix(ind, X);


      DenseMatrix<Precision> U;
      DenseVector<Precision> S;
      DenseVector<Precision> center;
     /* if(A.N() > 1000){
        RandomSVD<Precision> svd(A, node->phi.N()+6, 2, true);
        U = svd.U;
        S = svd.S;
        center = svd.c;
      }
      else{*/
        SVD<Precision> svd(A, true);
        U = svd.U;
        S = svd.S;
        center = svd.c;
        svd.Vt.deallocate();
      //}
      A.deallocate();

      if(U.N() > node->phi.N()){
        DenseMatrix<Precision> tmp = Linalg<Precision>::ExtractColumns(U, 0,
            node->phi.N());
        U.deallocate();
        U = tmp;
        S.shorten(node->phi.N()); 
      }
      else if( U.N() < node->phi.N() ){
        DenseMatrix<Precision> tmp(U.M(), node->phi.N());
        Linalg<Precision>::Zero(tmp);
        Linalg<Precision>::SetColumns(tmp, 0, U.N(), U, 0);
        U.deallocate();
        U = tmp;

        DenseVector<Precision> s(node->phi.N());
        Linalg<Precision>::Zero(s);
        for(int i=0; i<S.N(); i++){
          s(i) = S(i);
        }
        S.deallocate();
        S = s;
      }
      Linalg<Precision>::Scale(S, 1.0/sqrt(A.N()), S);
      for(int i=0; i<S.N(); i++){
        S(i) *= S(i);
      }

      Precision dist = Wasserstein<Precision>::distance2(U, S, center, node->phi,
          node->sigma2, node->center);
      U.deallocate();
      S.deallocate();
      center.deallocate();



      wdists.push_back(dist);
      for(std::list<int>::iterator it=ind.begin(); it != ind.end(); ++it){
        int i= *it;
        sumD(0, i) += dist; 
      }      

     std::vector<int> &nodeInd = node->getPoints(); 
     for(std::vector<int>::iterator it=nodeInd.begin(); it != nodeInd.end(); ++it){
        int i= *it;
        sumD(1, i) += dist; 
      }

      if(node->getChildren().size() != 0){
        std::list<int> i1;
        std::list<int> i2;
        DenseVector<Precision> tmp(X.M());
        for(std::list<int>::iterator it=ind.begin(); it != ind.end(); ++it){
          int i= *it;
          Linalg<Precision>::ExtractColumn(X, i, tmp);
          if(node->split(tmp) > 0){ 
            i1.push_back(i);
          }
          else{
            i2.push_back(i);
          }
        }
        tmp.deallocate();

        indices.push_back(i1);
        indices.push_back(i2);
      }
    }
    else{
      wdists.push_back(-1);
      if(node->getChildren().size() != 0){
         std::list<int> empty;
         indices.push_back(empty);
         indices.push_back(empty);
      }
    }

    std::list< GMRANode* > children = node->getChildren();
    for(std::list< GMRANode * >::iterator it =
        children.begin(); it != children.end(); ++it){
      IDcount++;
      nodeID.push_back(IDcount);
      nodes.push_back(*it);
    }


    if(children.size() == 0){
      std::vector<int> pts = node->getPoints();
      for(int i=0; i<pts.size(); i++){
        leafID(pts[i]) = curID;
      }
    }
  }


  DenseVector<Precision> wd = Linalg<Precision>::ToVector(wdists);
  std::stringstream ss1;
  ss1 << outFileName << "-wdist.data";
  LinalgIO<Precision>::writeVector(ss1.str(), wd);

  std::stringstream ss2;
  ss2 << outFileName << "-lid.data";
  LinalgIO<int>::writeVector(ss2.str(), leafID);

  std::stringstream ss3;
  ss3 << outFileName << "-sum-wdist.data";
  LinalgIO<Precision>::writeMatrix(ss3.str(), sumD);

  return 0;

}
