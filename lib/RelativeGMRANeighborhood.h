//author: Samuel Gerber

#ifndef RELATIVEGMRANEIGHBORHOOD_H
#define RELATIVEGMRANEIGHBORHOOD_H


#include "GMRATree.h"
#include "GMRANeighborhood.h"





//Neighborhood search using only the GMRA tree structure
template <typename TPrecision>
class RelativeGMRANeighborhood : public GMRANeighborhood<TPrecision>{
  
  private:
    typedef typename GMRANeighborhood<TPrecision>::Neighbor Neighbor;
    typedef typename GMRANeighborhood<TPrecision>::NeighborList NeighborList;
    typedef typename NeighborList::iterator NeighborListIterator;
    
    typedef typename GMRANode<TPrecision>::NodeVector  NodeVector;
    typedef typename NodeVector::iterator NodeVectorIterator;
    

  public:

    int neighbors(GMRANode<TPrecision> *x, TPrecision eps,
        NeighborList &collected, int stopScale) const{
      
      NeighborList nodes;
      GMRANode<TPrecision> *root = this->tree->getRoot();
      TPrecision d = this->dist->distance(x, root);
      nodes.push_back( Neighbor(d, root) );
      
      TPrecision xr = x->getRadius();

      int nCollected = 0;
      while( !nodes.empty()  ){
        
        Neighbor &n = nodes.front();

        NodeVector kids = n.second->getChildren();
        
        if( n.second->isStop() || kids.size() == 0 || n.second->getScale() == stopScale ){
          if(n.first <= eps + n.second->getRadius() + xr){
            collected.push_back(n);
            nCollected++;
          }
        }
        else if(n.second->getScale() < stopScale) {
          //For each kid check if nearest neighbors within epsilon are possible.
          for(NodeVectorIterator it = kids.begin(); it != kids.end(); ++it){
            GMRANode<TPrecision> *kid = *it;
            TPrecision d = this->dist->distance(kid, x);
            if(d <= kid->getRadius() + xr  + eps){
              nodes.push_back( Neighbor(d, kid) );
            }
          }
        }
        
        nodes.pop_front();
      }


      return nCollected;
    };




    
    RelativeGMRANeighborhood(GMRATree<TPrecision> *t, NodeDistance<TPrecision>
        *d) : GMRANeighborhood<TPrecision>(t, d){
    };

    ~RelativeGMRANeighborhood(){
    };
    


};



#endif 
