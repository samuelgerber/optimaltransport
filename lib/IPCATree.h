#ifndef IPCATREE_H
#define IPCATREE_H

#include "GMRATree.h"
#include "IPCANode.h"
#include "IPCANodeFactory.h"

#include <map>
#include <iostream>
#include <fstream>

#include <set>


template <typename TPrecision>
class IPCATree : public GMRATree<TPrecision>{

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
    



  private:

    IPCANodeFactory<TPrecision> *nodeFactory;

    

    TPrecision rootVariance;
    TPrecision rootRadius;
    TPrecision rootMSE;
    int nPoints;


    void buildTreeRecursive(IPCANode<TPrecision> *node){

      //Stop tree building?
      if( nodeFactory->isStop(node, nPoints, rootRadius, rootMSE, rootVariance) ){
        return;
      }



      //Create sub partitions
      int size =  node->getMaxKids() ;
      std::vector< std::vector<int> > children(size);
      

      std::vector<int> &nodePts = node->getPoints();
      for(std::vector<int>::iterator it = nodePts.begin();
          it!=nodePts.end(); ++it){
        int childIndex = node->getChildIndex( this->data->getPoint(*it) );
        children[childIndex].push_back(*it);
      }

      for(int i=0; i< children.size(); i++){
        if(children[i].size() > 0){

          IPCANode<TPrecision> *n = 
            dynamic_cast<IPCANode<TPrecision>*>( nodeFactory->createNode(children[i]) );

          node->addChild(n, i);
          buildTreeRecursive(n);
        }
      }

    };





  public:





    IPCATree(GMRADataObject<TPrecision> *D, IPCANodeFactory<TPrecision> *nf) :
      GMRATree<TPrecision>(D),  nodeFactory(nf) { 
    };



    ~IPCATree(){
      delete nodeFactory;
    };




    //TODO: implment tree updating for adding points when root is not NULL
    void addPoints(std::vector<int> &pts){


      IPCANode<TPrecision> *root = dynamic_cast< IPCANode<TPrecision>* >( nodeFactory->createNode(pts) );
      rootVariance = root->getTotalVariance();
      rootMSE = rootVariance - root->getSigma().array().square().sum(); 
      rootRadius = root->getL2Radius();
      nPoints = pts.size();

      buildTreeRecursive( (IPCANode<TPrecision>*) root ); 

      this->setRoot(root);
      this->setupParents();


    };





    std::vector<GMRANode<TPrecision> *> getLeafPath(const VectorXp &x ) {
   
      GMRANode<TPrecision> *node = this->getRoot();

      std::vector<GMRANode<TPrecision> *> path;
      while(node != NULL ){
        path.push_back( node );
        node = node->findDescendant( x );
      }

      return path;
    }; 






};


#endif
