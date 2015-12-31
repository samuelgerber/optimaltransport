#ifndef IPCAGWT_H
#define IPCAGWT_H


#include "GWT.h"
#include "IPCATree.h"



template <typename TPrecision> 
class IPCAGWTNode : public GWTNode<TPrecision>{
  public:  
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
  
  private:
    IPCANode<TPrecision> *node;

  public:

    IPCAGWTNode<TPrecision>(int id, GWTNode<TPrecision> *p, IPCANode<TPrecision>
        *n ) : GWTNode<TPrecision>(n->getPhi(), n->getCenter(), id, p, n), node(n) {
    };

    virtual ~IPCAGWTNode(){};


    virtual MatrixXp &getPhi(){
      return node->getPhi();
    };

    virtual VectorXp &getSigma(){
      return node->getSigma();
    };

    virtual VectorXp &getCenter(){
      return node->getCenter();
    };

};


template <typename TPrecision> 
class IPCAGWTDecorator : public Decorator<TPrecision>{
  private:
    int id;
  public:

    IPCAGWTDecorator(){
      id=0;
    };

    GMRANode<TPrecision> *decorate(GMRANode<TPrecision> *node, GMRANode<TPrecision> *parent){
      GWTNode<TPrecision> *pgwt = dynamic_cast<GWTNode<TPrecision> *>(parent);
      IPCANode<TPrecision> *ipcanode = dynamic_cast<IPCANode<TPrecision> *>( node );

      GWTNode<TPrecision> *gwt = new IPCAGWTNode<TPrecision>(id, pgwt, ipcanode);
      ++id;
      return gwt;
    };

};


template <typename TPrecision> 
class IPCAGWT : public GWT<TPrecision>{

  protected:

    void setupGWT(GMRATree<TPrecision> *tree){
      IPCAGWTDecorator<TPrecision> decorator;
      tree->decorate(decorator);
    };  


};


#endif
