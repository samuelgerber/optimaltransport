#ifndef IKMTREE_H
#define IKMTREE_H

#include "GMRATree.h"
#include "IKMNode.h"
#include "IKMKmeansData.h"
#include "Kmeans.h"


#include <map>
#include <iostream>
#include <fstream>

#include <set>



template <typename TPrecision>
class IKMTree : public GMRATree<TPrecision>{

  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
    

    enum StoppingCriterium {R2, MSE, RADIUS, RELATIVE_RADIUS};
    enum SplitCriterium {FIXED, ADAPTIVE, ADAPTIVE_FIXED};

    IKMKmeansDataFactory<TPrecision> *dataFactory;
    SplitCriterium split;
    StoppingCriterium stop;
    TPrecision epsilon;
    int minPoints;
    int nKids;
    TPrecision threshold;
    int maxIter;




  private:
    
    TPrecision rootRadius;
    TPrecision rootMSE;
    int nPoints;




    KmeansData<TPrecision> *getKmeansData(std::vector<int> &pts){
      return dataFactory->createKmeansData(this->getDataObject(), pts);
    };




    std::vector< KmeansCenter<TPrecision> > runKmeans(std::vector<int> &pts,
        std::vector<VectorXp> &means, TPrecision radius){

      KmeansData<TPrecision> *tmpData = getKmeansData(pts);
      Kmeans<TPrecision> kmeans(maxIter, threshold);
      std::vector< KmeansCenter<TPrecision> > centers;
      if(split != FIXED){
        centers =  kmeans.run( radius/2 , nKids, means, *tmpData);
      }
      else{
        centers =  kmeans.run( means, *tmpData, minPoints);
      }
      delete tmpData;

      return centers;
    };




    std::vector< KmeansCenter<TPrecision> > runKmeans(std::vector<int> &pts,
        TPrecision radius){

      KmeansData<TPrecision> *tmpData = getKmeansData(pts);
      Kmeans<TPrecision> kmeans(maxIter, threshold);
      std::vector< KmeansCenter<TPrecision> > centers;
      if(split != FIXED){
        centers =  kmeans.run( radius/2, nKids, *tmpData);
      }
      else{
        centers =  kmeans.run( std::min((int)pts.size(), nKids), *tmpData, minPoints);
      }
      delete tmpData;

      return centers;
    };


    void buildTreeRecursive(IKMNode<TPrecision> *node, int scale){
#ifdef VERBOSE
      std::cout << "Node MSE : " << node->getMSE() << std::endl;
      std::cout << "Node size : " << node->getPoints().size() << std::endl;
      std::cout << "Node R^2 : " << node->getMSE()/rootMSE << std::endl;
      std::cout << "Node Kmeans radius : " << node->getKMRadius() << std::endl;
      std::cout << "Node relative Kmeans radius : " << node->getKMRadius()/rootRadius << std::endl;
#endif


      
      //Stop tree building?
      if(node->getPoints().size() <= std::max(1, minPoints) ){
        return;
      }
      if(stop == R2 && (node->getMSE() / rootMSE) < epsilon){
        return;
      }
      if(stop == MSE && node->getMSE()  < epsilon){
        return;
      }
      if(stop == RADIUS && node->getKMRadius()  < epsilon){
        return;
      }
      if(stop == RELATIVE_RADIUS && (node->getKMRadius()/rootRadius)  < epsilon){
        return;
      }


      //Create sub partitions
      std::vector<int> &nodePts = node->getPoints();
      std::vector< GMRANode<TPrecision> *> kids = node->getChildren();
      
      if(kids.size() == 0){
       
        TPrecision r = node->getKMRadius();
        if(split == ADAPTIVE_FIXED){
          r = rootRadius / pow(2.0, scale );
        }

        std::vector< KmeansCenter<TPrecision> > centers = runKmeans(nodePts, r);
       
        if(centers.size() < 2){
         return;
        } 
        for(int i=0; i< centers.size(); i++){
         
          if(centers[i].points.size() > 0 ){ 
            std::vector<int> kidPts(centers[i].points);
            for(int j=0; j<kidPts.size(); j++){
              kidPts[j] = nodePts[centers[i].points[j]];
            }

            IKMNode<TPrecision> *n = new IKMNode<TPrecision>(centers[i].center,
                kidPts, centers[i].radius, centers[i].mse);
            node->addChild(n);

            buildTreeRecursive(n, scale+1);
          }
        }
      
      }
      else{
        
        std::vector<VectorXp> means(kids.size());
        for(int i=0; i<means.size(); i++){
          means[i] = kids[i]->getCenter();
        }
        
        TPrecision r = node->getKMRadius();
        if(split == ADAPTIVE_FIXED){
          r = rootRadius / pow(2.0, scale );
        }
        std::vector< KmeansCenter<TPrecision> > centers = runKmeans(nodePts, means, r );
        if(centers.size() < 2){
         return;
        } 
        
        for(int i=0; i<kids.size(); i++){
          IKMNode<TPrecision> *node = dynamic_cast<IKMNode<TPrecision>*>( kids[i] );
          
          std::vector<int> kidPts(centers[i].points);
          for(int j=0; j<kidPts.size(); j++){
            kidPts[j] = nodePts[centers[i].points[j]];
          }

          node->setCenter(centers[i].center);
          node->setPoints(kidPts);
          node->setKMRadius(centers[i].radius);
          node->setMSE(centers[i].mse);
          
          buildTreeRecursive(node, scale+1);

        }

      }



    };







  public:





    IKMTree(GMRADataObject<TPrecision> *D) : GMRATree<TPrecision>(D) {
      stop = RELATIVE_RADIUS;
      split = ADAPTIVE;
      epsilon = 0.05;
      nKids = 16;
      threshold= 0.01;
      maxIter = 100;
      minPoints = 1;
    };




    virtual ~IKMTree(){
      delete dataFactory;
    };





    void setSplitCriterium(int s){
      split = ADAPTIVE;
      switch(s){
        case 1:
          split = FIXED;
          break;
        case 2:
          split = ADAPTIVE;
          break;
        case 3:
          split = ADAPTIVE_FIXED;
          break;
      }
    };



    void setStoppingCriterium(int s){
      stop = RADIUS;
      switch(s){
        case 1:
          stop = R2;
          break;
        case 2:
          stop = MSE;
          break;
        case 3:
          stop = RADIUS;
          break;
        case 4:
          stop = RELATIVE_RADIUS;
          break;
      }
    };







    void addPoints(std::vector<int> &pts){


      GMRANode<TPrecision> *root = this->getRoot();
      if(root != NULL){
        std::vector<int> &ppts = root->getPoints();
        pts.insert(pts.end(), ppts.begin(), ppts.end() );
      }

      KmeansData<TPrecision> *tmpData = getKmeansData(pts);
      VectorXp mean = tmpData->getMean( pts );
      
      rootRadius = 0;
      rootMSE = 0;
      for(int i = 0; i < pts.size(); i++){
        TPrecision tmp = tmpData->getSquaredSimilarity(i, mean);
        rootMSE += tmp;
        if(tmp > rootRadius){
          rootRadius = tmp;
        }
      
      }
      rootMSE /= pts.size();
      rootRadius = sqrt(rootRadius);

      delete tmpData;



      if(root == NULL){
        root = new IKMNode<TPrecision>(mean, pts, rootRadius, rootMSE);
      }
      else{
        IKMNode<TPrecision> *node = dynamic_cast<IKMNode<TPrecision>*>( root );
        node->setCenter(mean);
        node->setPoints(pts);
        node->setKMRadius(rootRadius);
        node->setMSE(rootMSE);
      };

      IKMNode<TPrecision> *node = dynamic_cast<IKMNode<TPrecision>*>( root );
      buildTreeRecursive( node , 0); 

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
