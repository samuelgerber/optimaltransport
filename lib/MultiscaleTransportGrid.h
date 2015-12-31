#ifndef MULTISCALETRANSPORTGRID_H
#define MULTISCALETRANSPORTGRID_H


#include "MultiscaleTransport.h"

#include <map>


template <typename TPrecision>
class Index{
  public:
    virtual TPrecision distance(Index<TPrecision> &other, int p) const = 0;
    virtual bool operator == (const Index<TPrecision> &other) const = 0;
    virtual bool operator <  (const Index<TPrecision> &other) const = 0;
    virtual bool operator >  (const Index<TPrecision> &other) const = 0;

    std::list<Index> getNeighborhood(){

};



template <typename TPrecision>
class Index2d : public Index<TPrecision> {
  public:
    int i, j;
    Index2d():i(0), j(0){};
    Index2d(int ii, int jj):i(ii), j(jj){};
    
    virtual TPrecision distance(Index<TPrecision> &other, int p) const{
      int x = other.i  - i;  
      int y = other.i  - i; 
      return powf(x*x +y*y, p/2.0); 
    };

    virtual bool operator == (const Index<TPrecision> &other) const{
      return i==other.i && j==other.j;
    };

    virtual bool operator <  (const Index<TPrecision> &other) const{
      if(i < other.i){
        return true;
      }
      else if( i = =other.i){
        return j < other.j;
      }
      return false;
    };

    virtual bool operator >  (const Index<TPrecision> &other) const{
      return other < *this;
    };

   
}



template <typename TPrecision, typename TIndex>
class HistogramEntry : public TransportNode<TPrecision>{
  
private:
  TIndex i;
  Histogram2d<TPrecision, TIndex> *pHist;
  Histogram2d<TPrecision, TIndex> *hist;
  Histogram2d<TPrecision, TIndex> *cHist;

public:

    HistogramEntry(TIndex i, Histogram2d<TPrecision, TIndex> *current,
        Histogram2d<TPrecision, TIndex> *parent) : hist(current), pHist(parent){
      i = index;

      children = child->getNeighborhood(i, 
    };

    virtual TPrecision getMass() const { 
      return hist->p(index);
    };
    
    virtual TPrecision getTransportCost(TransportNode<TPrecision> &other, int p) const {
       HistogramEntry<TPrecision, TIndex> &o = (HistogramEntry<TPrecision, Index> other;
       return i.distance(o.i, p);  
    };

    virtual TPrecision getRadius() const{ 
      return hist->getRadius(); 
    };
    
    virtual TransportNodeVector &getChildren(){ 
      return children; 
    };

    void setChildren(TransportNodeList &kids){ children= kids 

    };
    
    virtual TransportNode<TPrecision> *getParent() const{ return parent; };
    void setParent(TransportNode<TPrecision> *p){ parent = p; };
    
    virtual bool operator == (const TransportNode<TPrecision> &other) const{
      return i == other.i;
    };

    virtual bool operator <  (const TransportNode<TPrecision> &other) const{
      return i < other.i;
    };

    virtual bool operator >  (const TransportNode<TPrecision> &other) const{
      return i > other.i;
    };
  

};


template <typename TPrecision>
class Histogram2d : public MultiscaleTransportLevel{
  public:
    typedef typename DenseMatrix< HistogramEntry<TPrecision, Index2d>* > Hist;
    Hist hist;
    
    typedef typename TransportNode<TPrecision>::TransportNodeVector
      TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;





    Histogram2d(int m, int n){
      hist = Hist(m, n);     
    };

    virtual ~Histogram2d(){
      delete hist;
    };

    HistogramEntry<TPrecision, Index2d> *getEntry(Index2d i){
      return hist(i.i, i.j);
    };

    void insertEntry(Index2d i, TPrecision mass, TPrecision radius){
      hist(i.i, i.j) = new HistogramEntry<TPrecision, TIndex>(i, mass, radius);
    };

    virtual TransportNodeVector &getNodes() const{
      TransportNodeVector nodes( hist.M()*hist.N(), hist.data() )
      return nodes;
    };


    virtual TransportNodeVector getNeighborhood(TransportNode<tPrecision> *node, TPrecision eps) const{
      TPrecision r = node->getRadius();
      int n = floor(2*r/eps);
      TransportNodeVector nodes;
      HistogramEntry<TPrecision> *entry = (HistogramEntry<TPrecision,
          Index2d> *) node;

      for(int i= entry->i.i -n; i<= entry.i.i+n; i++){
        if(i < 0 || i >= hist.M() ){ continue; }

        for(int j= entry->i.j -n; i<= entry.i.j+n; j++){
          if(j < 0 || j >= hist.N() ){ continue; }
          
          nodes.push_back( hist(i, j) );
        }
      }

      return nodes;
    };


    static std::vector< Histogram2d<TPrecision> * >
      buildMultiscaleHistogram2d(DenseMatrix<TPrecision> A){
        return buildPyramid(A);
    };


  private:


    Hist hist;

    static std::list< Histogram2d<TPrecision> * > buildPyramid(DenseMatrix<TPrecision> X){
      
      Histogram2d<TPrecision> *H = new Histogram2d<TPrecision>();

      for(int i=0; i < X.M(); i++){
        for(int j=0; j < X.N(); j++){
          if(X(i, j)  != 0){

          }
        }
      }
      H->radius = 0.5;

      std::list< Histogram2D<TPrecision> * > pyramid;
      pyramid.push_front(H); 
      while( std::min(H->p.M(), H->p.N()) > 2 ){

        Histogram2d<TPrecision> *H2 = new Histogram2d<TPrecision>(H->p.M()/2, H->p.N()/2);
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

};





